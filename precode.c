#include "precode.h"
#include "heap.h"

static int rowcmp(const rowstat *a, const rowstat *b) {
  if (a->nz == b->nz)
    return a->idx > b->idx;
  return a->nz > b->nz;
}

HEAP_GENERATE(rh, rowstat, rowcmp)

static void precode_matrix_permute(octmat *D, int P[], int n) {
  for (int i = 0; i < n; i++) {
    int at = i, mark = -1;
    while (P[at] >= 0) {
      oswaprow(om_P(*D), i, P[at], D->cols);
      int tmp = P[at];
      P[at] = mark;
      at = tmp;
    }
  }
}

static void precode_matrix_apply_sched(octmat *D, schedule *S) {
  for (int i = 0; i < kv_size(S->ops); i++) {
    sched_op op = kv_A(S->ops, i);
    if (op.beta)
      oaxpy(om_P(*D), om_P(*D), op.i, op.j, D->cols, op.beta);
    else
      oscal(om_P(*D), op.i, D->cols, op.j);
  }
}

static void precode_matrix_make_identity(spmat *A, int dim, int m, int n) {
  for (int diag = 0; diag < dim; diag++)
    spmat_push(A, m + diag, n + diag);
}

static void precode_matrix_make_LDPC1(spmat *A, int S, int B) {
  for (int col = 0; col < B; col++) {
    int submtx = col / S;
    int b1 = (col % S);
    int b2 = (col + submtx + 1) % S;
    int b3 = (col + 2 * (submtx + 1)) % S;
    spmat_push(A, b1, col);
    spmat_push(A, b2, col);
    spmat_push(A, b3, col);
  }
}

static void precode_matrix_make_LDPC2(spmat *A, int W, int S, int P) {
  for (int idx = 0; idx < S; idx++) {
    int b1 = idx % P;
    int b2 = (idx + 1) % P;
    spmat_push(A, idx, W + b1);
    spmat_push(A, idx, W + b2);
  }
}

static octmat precode_matrix_make_HDPC(params *P) {
  int m = P->H;
  int n = P->Kprime + P->S;

  octmat HDPC = OM_INITIAL;
  assert(m > 0 && n > 0);
  om_resize(&HDPC, m, n);

  for (int row = 0; row < m; row++)
    om_A(HDPC, row, n - 1) = OCT_EXP[row];

  for (int col = n - 2; col >= 0; col--) {
    for (int row = 0; row < m; row++)
      om_A(HDPC, row, col) = OCT_EXP[OCT_LOG[om_A(HDPC, row, col + 1)] + 1];
    int b1 = rnd_get(col + 1, 6, m);
    int b2 = (b1 + rnd_get(col + 1, 7, m - 1) + 1) % m;
    om_A(HDPC, b1, col) = om_A(HDPC, b1, col) ^ 1;
    om_A(HDPC, b2, col) = om_A(HDPC, b2, col) ^ 1;
  }
  return HDPC;
}

static void precode_matrix_make_G_ENC(spmat *A, params *P) {
  for (int row = P->S + P->H; row < P->L; row++) {
    uint32_t isi = (row - P->S) - P->H;
    int_vec idxs = params_get_idxs(isi, P);
    for (int idx = 0; idx < kv_size(idxs); idx++)
      spmat_push(A, row, kv_A(idxs, idx));
    kv_destroy(idxs);
  }
}

spmat *precode_matrix_gen(params *P, int overhead) {
  spmat *A = spmat_new(P->L + overhead, P->L);
  precode_matrix_make_LDPC1(A, P->S, P->B);
  precode_matrix_make_identity(A, P->S, 0, P->B);
  precode_matrix_make_LDPC2(A, P->W, P->S, P->P);
  precode_matrix_make_G_ENC(A, P);
  return A;
}

static void precode_matrix_sort(params *P, spmat *A, schedule *S) {
  rowstat sorted[A->rows];
  for (int row = 0; row < A->rows; row++) {
    sorted[row].idx = row;
    sorted[row].nz = spmat_nnz(A, row, 0, A->cols - P->P);
  }
  for (int row = P->S; row < P->S + P->H; row++)
    sorted[row].nz = A->cols - P->P;

  HEAPIFY(rh, sorted, A->rows);
  size_t n = A->rows;
  size_t row = 0;
  while (n > 0) {
    S->d[row] = sorted[0].idx;
    S->nz[S->d[row]] = sorted[0].nz;
    HEAP_POP(rh, sorted, n);
    n--;
    row++;
  }
  for (int row = 0; row < A->rows; row++)
    S->di[S->d[row]] = row;
}

static int precode_matrix_choose(int V0, int Vrows, int Srows, int Vcols,
                                 schedule *S) {
  int chosen = Vrows, r = Vcols + 1;
  for (int row = V0; row < Srows; row++) {
    int nz = S->nz[S->d[row]];
    if (nz > 0 && nz < r) {
      chosen = row;
      r = nz;
      if (nz < 2)
        break;
    }
  }
  return chosen;
}

static int precode_row_nz_fwd(spmat *A, int row, int s, int e, schedule *S) {
  int min = e;
  int_vec rs = A->idxs[S->d[row]];
  for (int it = 0; it < kv_size(rs); it++) {
    int col = S->ci[kv_A(rs, it)];
    if (col >= s && col < e && col < min) {
      min = col;
    }
  }
  return min;
}

static int precode_row_nz_rev(spmat *A, int row, int s, int e, schedule *S) {
  int tailsz = 10, tail[10]; // keep track of last N and find first zero col
  for (int t = 0; t < tailsz; t++)
    tail[t] = 0;
  int_vec rs = A->idxs[S->d[row]];
  for (int it = 0; it < kv_size(rs); it++) {
    int col = S->ci[kv_A(rs, it)];
    if (col > (e - tailsz) && col < e) {
      tail[col - (e - tailsz)] = col;
    }
  }
  for (int j = (tailsz - 1); j >= 0; j--) {
    if (tail[j] == 0) {
      return j + e - tailsz;
    }
  }
  return s;
}

static void precode_matrix_swap_first_col(spmat *A, int V0, int Vcols,
                                          schedule *S) {
  int *c = S->c, *ci = S->ci;
  int first = precode_row_nz_fwd(A, V0, V0, V0 + Vcols, S);
  if (first != V0) {
    TMPSWAP(int, c[V0], c[first]);
    TMPSWAP(int, ci[c[V0]], ci[c[first]]);
  }
}

static void precode_matrix_swap_tail_cols(spmat *A, int row, int s, int e,
                                          schedule *S) {
  while (s < e) {
    int swap1 = precode_row_nz_fwd(A, row, s, e, S);
    if (swap1 == e)
      return;
    int swap2 = precode_row_nz_rev(A, row, s, e, S);
    if (swap2 == s)
      return;
    if (swap1 >= swap2)
      return;
    TMPSWAP(int, S->c[swap1], S->c[swap2]);
    TMPSWAP(int, S->ci[S->c[swap1]], S->ci[S->c[swap2]]);
    s = swap1;
    e = swap2;
  }
}

static void precode_matrix_update_nnz(spmat *AT, int V0, int Vcols, int r,
                                      schedule *S) {
  int_vec cs = AT->idxs[S->c[V0]];
  for (int it = 0; it < kv_size(cs); it++) {
    int row = kv_A(cs, it);
    --S->nz[row];
  }
  for (int col = 0; col < r - 1; col++) {
    cs = AT->idxs[S->c[V0 + Vcols - col - 1]];
    for (int it = 0; it < kv_size(cs); it++) {
      int row = kv_A(cs, it);
      --S->nz[row];
    }
  }
}

static bool precode_matrix_precond(params *P, spmat *A, spmat *AT,
                                   schedule *S) {
  int i = 0, u = P->P, rows = A->rows, Srows = A->rows - P->H, cols = A->cols;
  int *d = S->d, *di = S->di;

  while (i + u < P->L) {
    int Vrows = rows - i, Vcols = cols - i - u, V0 = i;
    int chosen = precode_matrix_choose(V0, Vrows, Srows, Vcols, S);
    if (chosen >= Srows)
      return false;
    if (V0 != chosen) {
      TMPSWAP(int, d[V0], d[chosen]);
      TMPSWAP(int, di[d[V0]], di[d[chosen]]);
    }
    int r = S->nz[d[V0]];
    precode_matrix_swap_first_col(A, V0, Vcols, S);
    precode_matrix_swap_tail_cols(A, V0, V0 + 1, V0 + Vcols, S);
    precode_matrix_update_nnz(AT, V0, Vcols, r, S);
    i++;
    u += r - 1;
  }
  S->i = i;
  S->u = u;
  return true;
}

static void precode_matrix_unwind_X(params *P, wrkmat *U, schedule *S) {
  for (int i = S->Xe; i >= S->Xs; i--) {
    sched_op op = kv_A(S->ops, i);
    wrkmat_axpy(U, op.i, op.j, op.beta);
    sched_push(S, op.i, op.j, op.beta);
  }
}

static void precode_matrix_rewind_X(params *P, wrkmat *U, schedule *S) {
  for (int i = S->Xs; i <= S->Xe; i++) {
    sched_op op = kv_A(S->ops, i);
    sched_push(S, op.i, op.j, op.beta);
  }
}

static void precode_matrix_fwd_GE(wrkmat *U, schedule *S, spmat *AT, int s,
                                  int e) {
  int *c = S->c, *d = S->d, *di = S->di;
  for (int row = 0; row < S->i; row++) {
    int mv = s < row ? row : s;
    int_vec cs = AT->idxs[c[row]];
    for (int it = 0; it < kv_size(cs); it++) {
      int tmp = kv_A(cs, it);
      int h = di[tmp];
      if (h > mv && h < e) {
        wrkmat_axpy(U, d[h], d[row], 1);
        sched_push(S, d[h], d[row], 1);
      }
    }
  }
}

static void precode_matrix_fill_U_upper(wrkmat *U, spmat *A, spmat *AT,
                                        schedule *S) {
  for (int i = 0; i < A->rows; i++) {
    int_vec rs = A->idxs[i];
    for (int it = 0; it < kv_size(rs); it++) {
      int col = S->ci[kv_A(rs, it)];
      if (col >= S->i)
        wrkmat_set(U, i, col - S->i, 1);
    }
  }
}

static void precode_matrix_fill_U_lower(wrkmat *U, octmat *HDPC, params *P,
                                        schedule *S) {
  octmat UL = OM_INITIAL;
  om_resize(&UL, P->H, S->u);
  for (int row = 0; row < UL.rows; row++) {
    for (int col = 0; col < UL.cols - P->H; col++)
      om_A(UL, row, col) =
          om_A(*HDPC, row, S->c[HDPC->cols - (S->u - P->H) + col]);
    om_A(UL, row, row + (UL.cols - P->H)) = 1;
  }
  wrkmat_assign_block(U, &UL, P->S, 0, P->H, S->u);
}

static wrkmat *precode_matrix_make_U(params *P, spmat *A, spmat *AT,
                                     schedule *S) {
  int rows = A->rows, row_start = S->i;
  int *c = S->c, *d = S->d;

  wrkmat *U = wrkmat_new(rows, S->u);
  precode_matrix_fill_U_upper(U, A, AT, S);

  // forward GE on U upper
  precode_matrix_fwd_GE(U, S, AT, 0, S->i);
  S->Xs = 0;
  S->Xe = kv_size(S->ops) - 1;
  precode_matrix_fwd_GE(U, S, AT, S->i - 1, rows - P->H);

  // build U lower from rightmost cols of HDPC and I_H
  octmat HDPC = precode_matrix_make_HDPC(P);
  precode_matrix_fill_U_lower(U, &HDPC, P, S);

  // forward GE on HDPC rows
  for (int row = 0; row < row_start; row++) {
    for (int h = 0; h < P->H; h++) {
      uint8_t beta = om_A(HDPC, h, c[row]);
      if (beta) {
        wrkmat_axpy(U, d[rows - P->H + h], d[row], beta);
        sched_push(S, d[rows - P->H + h], d[row], beta);
      }
    }
  }
  om_destroy(&HDPC);

  return U;
}

static bool precode_matrix_solve(params *P, wrkmat *U, schedule *S) {
  int rows = U->rows, row_start = S->i;
  int *d = S->d, *di = S->di;
  uint8_t beta = 0;

  for (int row = row_start; row < P->L; row++) {
    int nzrow = row;
    int urow = row - S->i;

    for (; nzrow < rows; nzrow++) {
      beta = wrkmat_at(U, d[nzrow], urow);
      if (beta != 0)
        break;
    }

    if (nzrow == rows)
      return false;

    if (row != nzrow) {
      TMPSWAP(int, d[row], d[nzrow]);
      TMPSWAP(int, di[d[row]], di[d[nzrow]]);
    }

    if (beta > 1) {
      wrkmat_scal(U, d[row], OCT_INV[beta]);
      sched_push(S, d[row], OCT_INV[beta], 0);
    }

    for (int del_row = row + 1; del_row < rows; del_row++) {
      beta = wrkmat_at(U, d[del_row], urow);
      if (beta == 0)
        continue;
      wrkmat_axpy(U, d[del_row], d[row], beta);
      sched_push(S, d[del_row], d[row], beta);
    }
  }
  return true;
}

static void precode_matrix_backsolve(params *P, spmat *AT, wrkmat *U,
                                     schedule *S) {
  int *c = S->c, *d = S->d, row_start = S->i;

  for (int row = P->L - 1; row >= row_start; row--) {
    int_vec cs = AT->idxs[c[row]];
    for (int it = 0; it < kv_size(cs); it++) {
      int del_row = S->di[kv_A(cs, it)];
      if (del_row < S->i) {
        sched_push(S, d[del_row], d[row], 1);
      }
    }
    for (int del_row = S->i; del_row < row; del_row++) {
      uint8_t beta = wrkmat_at(U, d[del_row], row - S->i);
      if (beta == 0)
        continue;
      sched_push(S, d[del_row], d[row], beta);
    }
  }
}

static void *precode_matrix_cleanup(spmat *A, spmat *AT, schedule *S,
                                    wrkmat *U) {
  spmat_free(A);
  spmat_free(AT);
  if (S)
    sched_free(S);
  if (U)
    wrkmat_free(U);
  return NULL;
}

schedule *precode_matrix_invert(params *P, spmat *A) {
  int rows = A->rows, cols = A->cols;
  schedule *S = sched_new(rows, cols, P->L);
  wrkmat *U = NULL;

  precode_matrix_sort(P, A, S);
  spmat *AT = spmat_transpose(A);

  if (!precode_matrix_precond(P, A, AT, S))
    return precode_matrix_cleanup(A, AT, S, U);

  U = precode_matrix_make_U(P, A, AT, S);
  if (U == NULL)
    return precode_matrix_cleanup(A, AT, S, U);

  if (!precode_matrix_solve(P, U, S))
    return precode_matrix_cleanup(A, AT, S, U);

  precode_matrix_unwind_X(P, U, S);
  precode_matrix_backsolve(P, AT, U, S);
  precode_matrix_rewind_X(P, U, S);

  precode_matrix_cleanup(A, AT, NULL, U);

  return S;
}

bool precode_matrix_intermediate(params *P, spmat *A, octmat *D) {
  int rows = A->rows, cols = A->cols;
  if (rows <= 0)
    return false;

  schedule *S = precode_matrix_invert(P, A);
  if (S == NULL)
    return false;

  precode_matrix_apply_sched(D, S);
  precode_matrix_permute(D, S->di, rows);
  precode_matrix_permute(D, S->c, cols);

  sched_free(S);
  return true;
}
