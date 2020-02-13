#include <stdint.h>
#include <stdlib.h>
#include <string.h>

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

static void precode_matrix_apply_sched(octmat *D, schedule *S) {
  for (int i = 0; i < kv_size(S->ops); i++) {
    struct sch_op op = kv_A(S->ops, i);
    if (op.beta)
      oaxpy(om_P(*D), om_P(*D), op.i, op.j, D->cols, op.beta);
    else
      oscal(om_P(*D), op.i, D->cols, op.j);
  }
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

static void precode_matrix_make_identity(spmat *A, int dim, int m, int n) {
  for (int diag = 0; diag < dim; diag++) {
    spmat_push(A, m + diag, n + diag);
  }
}

static octmat precode_matrix_make_MT(int rows, int cols) {
  octmat MT = OM_INITIAL;
  om_resize(&MT, rows, cols);

  for (int col = 0; col < MT.cols - 1; col++) {
    int b1 = rnd_get(col + 1, 6, MT.rows);
    int b2 = (b1 + rnd_get(col + 1, 7, MT.rows - 1) + 1) % MT.rows;
    om_A(MT, b1, col) = 1;
    om_A(MT, b2, col) = 1;
  }
  for (int row = 0; row < MT.rows; row++) {
    om_A(MT, row, MT.cols - 1) = OCT_EXP[row];
  }
  return MT;
}

static octmat precode_matrix_make_HDPC(params *P) {
  int m = P->H;
  int n = P->Kprime + P->S;

  octmat MT = precode_matrix_make_MT(m, n);
  octmat MTxGAMMA = OM_INITIAL;
  om_resize(&MTxGAMMA, m, n);

  uint8_t gv[2 * n];
  memset(gv, 0, n);
  for (int col = n; col < 2 * n; col++) {
    gv[col] = OCT_EXP[col % OCT_EXP_SIZE];
  }

  uint8_t *ap, *cp = om_P(MTxGAMMA);
  for (int row = 0; row < m; row++, cp += MTxGAMMA.cols_al) {
    ap = om_P(MT) + (row * MT.cols_al);
    for (int idx = 0; idx < n; idx++) {
      oaxpy(cp, gv + idx, 0, 0, n, ap[idx]);
    }
  }

  om_destroy(&MT);
  return MTxGAMMA;
}

static void precode_matrix_make_G_ENC(spmat *A, params *P) {
  for (int row = P->S + P->H; row < P->L; row++) {
    uint32_t isi = (row - P->S) - P->H;
    int_vec idxs = params_get_idxs(isi, P);
    for (int idx = 0; idx < kv_size(idxs); idx++) {
      spmat_push(A, row, kv_A(idxs, idx));
    }
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

static void decode_patch(params *P, spmat *A, struct bitmask *mask,
                         repair_vec *repair_bin, size_t num_symbols) {
  size_t padding = P->Kprime - num_symbols;
  int num_gaps = bitmask_gaps(mask, num_symbols);
  int rep_idx = 0;
  for (int gap = 0; gap < P->L && num_gaps > 0; gap++) {
    if (bitmask_check(mask, gap))
      continue;
    int row = gap + P->H + P->S;
    spmat_clear_row(A, row);

    int_vec idxs =
        params_get_idxs(kv_A(*repair_bin, rep_idx++).esi + padding, P);
    for (int idx = 0; idx < kv_size(idxs); idx++) {
      spmat_push(A, row, kv_A(idxs, idx));
    }
    kv_destroy(idxs);
    num_gaps--;
  }

  for (int rep_row = P->L; rep_row < A->rows; rep_row++) {
    spmat_clear_row(A, rep_row);
    int_vec idxs =
        params_get_idxs(kv_A(*repair_bin, rep_idx++).esi + padding, P);
    for (int idx = 0; idx < kv_size(idxs); idx++) {
      spmat_push(A, rep_row, kv_A(idxs, idx));
    }
    kv_destroy(idxs);
  }
}

static int decode_nz_fwd(spmat *A, int row, int s, int e, schedule *S) {
  int *d = S->d, *ci = S->ci;
  int min = e;
  int_vec rs = A->idxs[d[row]];
  for (int it = 0; it < kv_size(rs); it++) {
    int col = ci[kv_A(rs, it)];
    if (col >= s && col < e && col < min) {
      min = col;
    }
  }
  return min;
}

static int decode_nz_rev(spmat *A, int row, int s, int e, schedule *S) {
  int *d = S->d, *ci = S->ci;
  int tailsz = 10, tail[10]; // keep track of last N and find first empty col

  for (int t = 0; t < tailsz; t++) {
    tail[t] = 0;
  }
  int_vec rs = A->idxs[d[row]];
  for (int it = 0; it < kv_size(rs); it++) {
    int col = ci[kv_A(rs, it)];
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

static void decode_rear_swap(spmat *A, int row, int s, int e, schedule *S) {
  while (s < e) {
    int swap1 = decode_nz_fwd(A, row, s, e, S);
    if (swap1 == e)
      return;
    int swap2 = decode_nz_rev(A, row, s, e, S);
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

static bool decode_amd(params *P, spmat *A, spmat *AT, schedule *S) {
  int i = 0, u = P->P, rows = A->rows, cols = A->cols;
  int *c = S->c, *d = S->d, *ci = S->ci, *di = S->di;

  int counts[rows];
  for (int row = 0; row < rows; row++) {
    counts[d[row]] = spmat_nnz(A, d[row], 0, cols - u);
  }

  while (i + u < P->L) {
    int Vrows = rows - i, Vcols = cols - i - u, V0 = i;
    int r = Vcols + 1, chosen = Vrows;

    for (int row = V0; row < rows - P->H; row++) {
      int nz = counts[d[row]];
      if (nz > 0 && nz < r) {
        chosen = row;
        r = nz;
        if (nz < 2)
          break;
      }
    }
    if (r == Vcols + 1) {
      return false;
    }
    if (V0 != chosen) {
      TMPSWAP(int, d[V0], d[chosen]);
      TMPSWAP(int, di[d[V0]], di[d[chosen]]);
    }
    // find first one
    int first = decode_nz_fwd(A, V0, V0, V0 + Vcols, S);
    if (first != V0) {
      TMPSWAP(int, c[V0], c[first]);
      TMPSWAP(int, ci[c[V0]], ci[c[first]]);
    }
    decode_rear_swap(A, V0, V0 + 1, V0 + Vcols, S);
    // decrement nz counts if row had nz at V0 or nz's at last r - 1 cols
    int_vec cs = AT->idxs[c[V0]];
    for (int it = 0; it < kv_size(cs); it++) {
      counts[kv_A(cs, it)]--;
    }
    for (int col = 0; col < (r - 1); col++) {
      cs = AT->idxs[c[V0 + Vcols - col - 1]];
      for (int it = 0; it < kv_size(cs); it++) {
        counts[kv_A(cs, it)]--;
      }
    }
    i++;
    u += r - 1;
  }
  S->i = i;
  S->u = u;
  return true;
}

static void decode_unwind_X(params *P, wrkmat *U, schedule *S) {
  for (int i = S->Xe; i >= S->Xs; i--) {
    struct sch_op op = kv_A(S->ops, i);
    wrkmat_axpy(U, op.i, op.j, op.beta);
    sched_push(S, op.i, op.j, op.beta);
  }
}

static void decode_rewind_X(params *P, wrkmat *U, schedule *S) {
  for (int i = S->Xs; i <= S->Xe; i++) {
    struct sch_op op = kv_A(S->ops, i);
    sched_push(S, op.i, op.j, op.beta);
  }
}

static void decode_fwd_GE(wrkmat *U, schedule *S, spmat *AT, int s, int e) {
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

static wrkmat *decode_make_U(params *P, spmat *A, spmat *AT, schedule *S) {
  int rows = A->rows, row_start = S->i;
  int *c = S->c, *d = S->d, *ci = S->ci;

  // build U upper from indexes
  wrkmat *U = wrkmat_new(rows, S->u);
  for (int i = 0; i < rows; i++) {
    int_vec rs = A->idxs[i];
    for (int it = 0; it < kv_size(rs); it++) {
      int col = ci[kv_A(rs, it)];
      if (col >= S->i) {
        wrkmat_set(U, i, col - S->i, 1);
      }
    }
  }

  // forward GE on U upper
  decode_fwd_GE(U, S, AT, 0, S->i);
  S->Xs = 0;
  S->Xe = kv_size(S->ops) - 1;
  decode_fwd_GE(U, S, AT, S->i - 1, rows - P->H);

  // build U lower from rightmost cols of HDPC and I_H
  octmat HDPC = precode_matrix_make_HDPC(P);
  octmat UL = OM_INITIAL;
  om_resize(&UL, P->H, S->u);
  for (int row = 0; row < UL.rows; row++) {
    for (int col = 0; col < UL.cols - P->H; col++) {
      om_A(UL, row, col) = om_A(HDPC, row, c[HDPC.cols - (S->u - P->H) + col]);
    }
    om_A(UL, row, row + (UL.cols - P->H)) = 1;
  }
  wrkmat_assign_block(U, &UL, P->S, 0, P->H, S->u);

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

static bool decode_solve(params *P, wrkmat *U, schedule *S) {
  int rows = U->rows, row_start = S->i;
  int *d = S->d, *di = S->di;
  uint8_t beta = 0;

  for (int row = row_start; row < P->L; row++) {
    int nzrow = row;
    int urow = row - S->i;

    for (; nzrow < rows; nzrow++) {
      beta = wrkmat_at(U, d[nzrow], urow);
      if (beta != 0) {
        break;
      }
    }

    if (nzrow == rows) {
      return false;
    }

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

static void decode_backsolve(params *P, spmat *AT, wrkmat *U, schedule *S) {
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

static void *inv_cleanup(spmat *A, spmat *AT, schedule *S, wrkmat *U) {
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

  // roll precode matrix to put G_ENC rows on top
  for (int row = 0; row < rows; row++) {
    S->d[row] = (row + P->S + P->H) % rows;
  }
  for (int i = 0; i < rows; i++) {
    S->di[S->d[i]] = i;
  }
  spmat *AT = spmat_transpose(A);

  if (!decode_amd(P, A, AT, S))
    return inv_cleanup(A, AT, S, U);

  U = decode_make_U(P, A, AT, S);
  if (U == NULL)
    return inv_cleanup(A, AT, S, U);

  if (!decode_solve(P, U, S))
    return inv_cleanup(A, AT, S, U);

  decode_unwind_X(P, U, S);
  decode_backsolve(P, AT, U, S);
  decode_rewind_X(P, U, S);

  inv_cleanup(A, AT, NULL, U);

  return S;
}

bool precode_matrix_intermediate1(params *P, spmat *A, octmat *D) {
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

void precode_matrix_fill_slot(params *P, octmat *D, uint32_t isi, uint8_t *ptr,
                              size_t len) {
  int_vec idxs = params_get_idxs(isi, P);
  for (int idx = 0; idx < kv_size(idxs); idx++) {
    oaddrow(ptr, om_P(*D), 0, kv_A(idxs, idx), len);
  }
  kv_destroy(idxs);
}

bool precode_matrix_intermediate2(params *P, spmat *A, octmat *D, octmat *M,
                                  repair_vec *repair_bin,
                                  struct bitmask *repair_mask,
                                  int num_symbols) {
  int num_gaps, gap = 0, row = 0;

  if (D->cols == 0) {
    return false;
  }

  decode_patch(P, A, repair_mask, repair_bin, num_symbols);

  if (!precode_matrix_intermediate1(P, A, D)) {
    return false;
  }

  num_gaps = bitmask_gaps(repair_mask, num_symbols);
  om_resize(M, num_gaps, D->cols);
  for (gap = 0; gap < num_symbols && num_gaps > 0; gap++) {
    if (bitmask_check(repair_mask, gap))
      continue;
    precode_matrix_fill_slot(P, D, gap, om_R(*M, row), M->cols);
    row++;
    num_gaps--;
  }
  return true;
}

bool precode_matrix_decode(params *P, octmat *D, octmat *M,
                           repair_vec *repair_bin,
                           struct bitmask *repair_mask) {
  int rep_idx, num_gaps, num_repair, overhead, skip = P->S + P->H;
  int num_symbols = P->K;

  num_repair = kv_size(*repair_bin);
  num_gaps = bitmask_gaps(repair_mask, num_symbols);

  if (num_gaps == 0)
    return true;

  if (num_repair < num_gaps)
    return false;

  overhead = num_repair - num_gaps;
  rep_idx = 0;

  if (D->rows < P->S + P->H + P->Kprime + overhead) {
    // overhead estimate was insuffucient, have to reallocate
    octmat X = OM_INITIAL;
    om_resize(&X, P->S + P->H + P->Kprime + overhead, D->cols);
    memcpy(X.data, D->data, D->rows * D->cols_al);
    uint8_t *tmp = D->data;
    D->data = X.data;
    X.data = tmp;
    om_destroy(&X);
  }
  D->rows = P->S + P->H + P->Kprime + overhead;

  for (int gap = 0; gap < num_symbols && rep_idx < num_repair; gap++) {
    if (bitmask_check(repair_mask, gap))
      continue;
    int row = skip + gap;
    repair_sym rs = kv_A(*repair_bin, rep_idx++);
    ocopy(om_P(*D), om_P(rs.row), row, 0, D->cols);
  }

  for (int row = skip + P->Kprime; rep_idx < num_repair; row++) {
    repair_sym rs = kv_A(*repair_bin, rep_idx++);
    ocopy(om_P(*D), om_P(rs.row), row, 0, D->cols);
  }

  spmat *A = precode_matrix_gen(P, overhead);
  bool precode_ok = precode_matrix_intermediate2(P, A, D, M, repair_bin,
                                                 repair_mask, num_symbols);
  return precode_ok;
}
