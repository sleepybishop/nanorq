#include "precode.h"

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

static void precode_matrix_apply_op(octmat *D, schedule *S, int i) {
  sched_op op = kv_A(S->ops, i);
  if (op.beta)
    oaxpy(om_P(*D), om_P(*D), op.i, op.j, D->cols, op.beta);
  else
    oscal(om_P(*D), op.i, D->cols, op.j);
}

static void precode_matrix_apply_sched(octmat *D, schedule *S) {
  for (int i = 0; i < S->marks[1]; i++)
    precode_matrix_apply_op(D, S, i);
  for (int i = S->marks[0]; i >= 0; i--)
    precode_matrix_apply_op(D, S, i);
  for (int i = S->marks[1]; i < kv_size(S->ops); i++)
    precode_matrix_apply_op(D, S, i);
  for (int i = 0; i <= S->marks[0]; i++)
    precode_matrix_apply_op(D, S, i);
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

  assert(m > 0 && n > 0);
  octmat HDPC = OM_INITIAL;
  om_resize(&HDPC, m, n);

  for (int row = 0; row < m; row++)
    om_A(HDPC, row, n - 1) = OCT_EXP[row];

  for (int col = n - 2; col >= 0; col--) {
    for (int row = 0; row < m; row++)
      om_A(HDPC, row, col) =
          (om_A(HDPC, row, col + 1) == 0)
              ? 0
              : OCT_EXP[OCT_LOG[om_A(HDPC, row, col + 1)] + 1];
    int b1 = rnd_get(col + 1, 6, m);
    int b2 = (b1 + rnd_get(col + 1, 7, m - 1) + 1) % m;
    om_A(HDPC, b1, col) ^= 1;
    om_A(HDPC, b2, col) ^= 1;
  }
  return HDPC;
}

static void precode_matrix_make_G_ENC(spmat *A, params *P) {
  for (int row = P->S + P->H; row < P->L; row++)
    params_set_idxs(row - P->S - P->H, P, &A->idxs[row]);
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
  for (int row = 0; row < A->rows; row++)
    S->d[row] = (row + P->S + P->H) % A->rows; // move HDPC to bottom
  for (int i = 0; i < A->rows; i++)
    S->di[S->d[i]] = i;
  for (int row = 0; row < A->rows; row++) {
    S->nz[S->d[row]] = spmat_nnz(A, S->d[row], 0, A->cols - P->P);
    if (S->nz[S->d[row]] == 0)
      S->nz[S->d[row]] = A->cols;
  }
}

/* shortcuts are taken here
 *  - component / original degree tracking is skipped for speed
 *  - it might result in decoding failures until more symbols are added
 */
static int precode_matrix_choose(int V0, int Vrows, int Srows, int Vcols,
                                 schedule *S, spmat *NZT) {
  int chosen = Vrows;
  for (int b = 1; b < 3; b++) {
    while (kv_size(NZT->idxs[b]) > 0) {
      chosen = kv_pop(NZT->idxs[b]);
      if (S->di[chosen] >= V0 && S->nz[chosen] == b)
        return S->di[chosen];
    }
  }
  return Srows;
}

int precode_row_nz_at(spmat *A, int row, int s, int e, schedule *S, int *at) {
  int r = 0;
  at[0] = at[1] = e;
  uint_vec rs = A->idxs[S->d[row]];
  for (int it = 0; it < kv_size(rs) && r < S->nz[S->d[row]]; it++) {
    int col = S->ci[kv_A(rs, it)];
    if (col >= s && col < e)
      at[r++] = col;
  }
  if (at[0] > at[1])
    TMPSWAP(int, at[0], at[1]);
  return r;
}

static int precode_matrix_swap_cols(spmat *A, int V0, int Vcols, schedule *S) {
  int *c = S->c, *ci = S->ci, ones[2], Vlast = V0 + Vcols - 1;
  int r = precode_row_nz_at(A, V0, V0, V0 + Vcols, S, ones);
  if (ones[0] != V0) {
    TMPSWAP(int, c[V0], c[ones[0]]);
    TMPSWAP(int, ci[c[V0]], ci[c[ones[0]]]);
  }
  if (r == 2 && ones[1] != Vlast) {
    TMPSWAP(int, c[Vlast], c[ones[1]]);
    TMPSWAP(int, ci[c[Vlast]], ci[c[ones[1]]]);
  }
  return r;
}

static void precode_matrix_update_nnz(spmat *AT, int V0, int Vcols, int r,
                                      schedule *S, spmat *NZT) {
  uint_vec cs = AT->idxs[S->c[V0]];
  for (int it = 0; it < kv_size(cs); it++) {
    int row = kv_A(cs, it);
    int nz = --S->nz[row];
    if (nz > 0 && nz < 3)
      spmat_push(NZT, nz, row);
  }
  for (int col = 0; col < r - 1; col++) {
    cs = AT->idxs[S->c[V0 + Vcols - col - 1]];
    for (int it = 0; it < kv_size(cs); it++) {
      int row = kv_A(cs, it);
      int nz = --S->nz[row];
      if (nz > 0 && nz < 3)
        spmat_push(NZT, nz, row);
    }
  }
}

static void precode_matrix_precond(params *P, spmat *A, spmat *AT,
                                   schedule *S) {
  int i = 0, u = P->P, rows = A->rows, Srows = A->rows - P->H, cols = A->cols;
  int *d = S->d, *di = S->di;

  spmat *NZT = spmat_new(3, rows);
  for (int row = 0; row < Srows; row++) {
    if (S->nz[S->d[row]] < 3)
      spmat_push(NZT, S->nz[S->d[row]], S->d[row]);
  }
  while (i + u < P->L) {
    int Vrows = rows - i, Vcols = cols - i - u, V0 = i;
    int chosen = precode_matrix_choose(V0, Vrows, Srows, Vcols, S, NZT);
    if (chosen >= Srows)
      break;
    if (V0 != chosen) {
      TMPSWAP(int, d[V0], d[chosen]);
      TMPSWAP(int, di[d[V0]], di[d[chosen]]);
    }
    int r = precode_matrix_swap_cols(A, V0, Vcols, S);
    precode_matrix_update_nnz(AT, V0, Vcols, r, S, NZT);
    i++;
    u += r - 1;
  }
  spmat_free(NZT);
  S->i = i;
  S->u = P->L - i;
}

static void precode_matrix_fwd_GE(wrkmat *U, schedule *S, spmat *AT, int s,
                                  int e) {
  int *c = S->c, *d = S->d, *di = S->di;
  for (int row = 0; row < S->i; row++) {
    int mv = s < row ? row : s;
    uint_vec cs = AT->idxs[c[row]];
    for (int it = 0; it < kv_size(cs); it++) {
      int tmp = kv_A(cs, it), h = di[tmp];
      if (h > mv && h < e) {
        wrkmat_axpy(U, tmp, d[row], 1);
        sched_push(S, tmp, d[row], 1);
      }
    }
  }
}

static void precode_matrix_fill_U(wrkmat *U, spmat *A, spmat *AT, schedule *S) {
  for (int i = 0; i < A->rows; i++) {
    uint_vec rs = A->idxs[i];
    for (int it = 0; it < kv_size(rs); it++) {
      int col = S->ci[kv_A(rs, it)];
      if (col >= S->i)
        wrkmat_set(U, i, col - S->i, 1);
    }
  }
}

static void precode_matrix_fill_HDPC(params *P, wrkmat *U, schedule *S) {
  octmat UL = OM_INITIAL, HDPC = precode_matrix_make_HDPC(P);
  om_resize(&UL, 2 * P->H, S->u);
  for (int row = 0; row < P->H; row++) {
    for (int col = 0; col < UL.cols - P->H; col++)
      om_A(UL, row, col) =
          om_A(HDPC, row, S->c[HDPC.cols - (S->u - P->H) + col]);
    om_A(UL, row, row + (UL.cols - P->H)) = 1; // I_H
  }
  wrkmat_assign_block(U, &UL, P->S, 0, P->H, S->u);
  for (int row = 0; row < S->i; row++) {
    for (int h = 0; h < P->H; h++) {
      uint8_t beta = om_A(HDPC, h, S->c[row]);
      if (beta) {
        wrkmat_axpy(U, S->d[U->rows - P->H + h], S->d[row], beta);
        sched_push(S, S->d[U->rows - P->H + h], S->d[row], beta);
      }
    }
  }
  om_destroy(&HDPC);
}

static wrkmat *precode_matrix_make_U(params *P, spmat *A, spmat *AT,
                                     schedule *S) {
  wrkmat *U = wrkmat_new(A->rows, S->u);
  precode_matrix_fill_U(U, A, AT, S);
  precode_matrix_fwd_GE(U, S, AT, 0, S->i);
  S->marks[0] = kv_size(S->ops) - 1;
  precode_matrix_fwd_GE(U, S, AT, S->i - 1, A->rows - P->H);
  return U;
}

static int precode_matrix_solve_gf2(params *P, wrkmat *U, schedule *S) {
  int *d = S->d, *di = S->di, row, nzrow, rows = U->rows - P->H;
  for (row = S->i; row < P->L; row++) {
    int col = row - S->i;
    for (nzrow = row; nzrow < rows; nzrow++)
      if (wrkmat_at(U, d[nzrow], col))
        break;
    if (nzrow == rows)
      break;
    if (row != nzrow) {
      TMPSWAP(int, d[row], d[nzrow]);
      TMPSWAP(int, di[d[row]], di[d[nzrow]]);
    }
    for (int del_row = row + 1; del_row < rows; del_row++) {
      if (wrkmat_at(U, d[del_row], col) == 0)
        continue;
      wrkmat_axpy(U, d[del_row], d[row], 1);
      sched_push(S, d[del_row], d[row], 1);
    }
  }
  return row;
}

static int precode_matrix_solve_gf256(params *P, wrkmat *U, schedule *S) {
  int *d = S->d, *di = S->di, row, nzrow, rows = U->rows;
  for (row = S->i; row < P->L; row++) {
    int col = row - S->i, beta = 0;
    for (nzrow = row; nzrow < rows; nzrow++) {
      beta = wrkmat_at(U, d[nzrow], col);
      if (beta != 0)
        break;
    }
    if (nzrow == rows)
      break;
    if (row != nzrow) {
      TMPSWAP(int, d[row], d[nzrow]);
      TMPSWAP(int, di[d[row]], di[d[nzrow]]);
    }
    if (beta > 1) {
      wrkmat_scal(U, d[row], OCT_INV[beta]);
      sched_push(S, d[row], OCT_INV[beta], 0);
    }
    for (int del_row = row + 1; del_row < rows; del_row++) {
      beta = wrkmat_at(U, d[del_row], col);
      if (beta == 0)
        continue;
      wrkmat_axpy(U, d[del_row], d[row], beta);
      sched_push(S, d[del_row], d[row], beta);
    }
  }
  return row;
}

static void precode_matrix_backsolve(params *P, spmat *AT, wrkmat *U,
                                     schedule *S) {
  int *c = S->c, *d = S->d;
  for (int row = P->L - 1; row >= S->i; row--) {
    uint_vec cs = AT->idxs[c[row]];
    for (int it = 0; it < kv_size(cs); it++) {
      int del_row = S->di[kv_A(cs, it)];
      if (del_row < S->i)
        sched_push(S, d[del_row], d[row], 1);
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
  schedule *S = sched_new(rows, cols, 3 * P->L);
  wrkmat *U = NULL;

  precode_matrix_sort(P, A, S);
  spmat *AT = spmat_transpose(A);

  precode_matrix_precond(P, A, AT, S);

  U = precode_matrix_make_U(P, A, AT, S);
  if (U == NULL)
    return precode_matrix_cleanup(A, AT, S, U);

  int rank = 0;
  if ((A->rows - P->H) >= P->L)
    rank = precode_matrix_solve_gf2(P, U, S);

  if (rank < P->L) {
    precode_matrix_fill_HDPC(P, U, S);
    rank = precode_matrix_solve_gf256(P, U, S);
    if (rank < P->L) {
      return precode_matrix_cleanup(A, AT, S, U);
    }
  }
  S->marks[1] = kv_size(S->ops) - 1;
  precode_matrix_backsolve(P, AT, U, S);
  precode_matrix_cleanup(A, AT, NULL, U);

  return S;
}

void precode_matrix_intermediate(params *P, octmat *D, schedule *S) {
  precode_matrix_apply_sched(D, S);
  int *rm = calloc(sizeof(int), S->rows);
  int *cm = calloc(sizeof(int), S->cols);
  memcpy(rm, S->di, sizeof(int) * S->rows);
  memcpy(cm, S->c, sizeof(int) * S->cols);
  precode_matrix_permute(D, rm, S->rows);
  precode_matrix_permute(D, cm, S->cols);
  free(rm);
  free(cm);
}
