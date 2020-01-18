#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <oblas.h>

#include "bitmask.h"
#include "params.h"
#include "precode.h"
#include "rand.h"
#include "sched.h"

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

static void precode_matrix_make_LDPC1(wrkmat *A, int S, int B) {
  int row, col;
  for (row = 0; row < S; row++) {
    for (col = 0; col < B; col++) {
      uint16_t submtx = col / S;
      if ((row == (col % S)) || (row == (col + submtx + 1) % S) ||
          (row == (col + 2 * (submtx + 1)) % S)) {
        wrkmat_set(A, row, col, 1);
      }
    }
  }
}

static void precode_matrix_make_LDPC2(wrkmat *A, uint16_t skip, uint16_t rows,
                                      uint16_t cols) {
  for (int row = 0; row < rows; row++) {
    uint16_t start = row % cols;
    for (int col = 0; col < cols; col++) {
      uint8_t val = (col == start || col == (start + 1) % cols) > 0;
      wrkmat_set(A, row, skip + col, val);
    }
  }
}

static void precode_matrix_make_identity(wrkmat *A, int size, int skip_row,
                                         int skip_col) {
  for (int diag = 0; diag < size; diag++) {
    wrkmat_set(A, skip_row + diag, skip_col + diag, 1);
  }
}

static octmat precode_matrix_make_MT(int rows, int cols) {
  octmat MT = OM_INITIAL;
  om_resize(&MT, rows, cols);

  for (int row = 0; row < MT.rows; row++) {
    int col;
    for (col = 0; col < MT.cols - 1; col++) {
      uint32_t tmp = rnd_get(col + 1, 6, MT.rows);
      if ((row == tmp) ||
          (row == (tmp + rnd_get(col + 1, 7, MT.rows - 1) + 1) % MT.rows)) {
        om_A(MT, row, col) = 1;
      } else {
        om_A(MT, row, col) = 0;
      }
    }
    om_A(MT, row, col) = OCT_EXP[row];
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
  for (int col = 0; col < n; col++) {
    gv[col] = 0;
  }
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

static void precode_matrix_make_G_ENC(wrkmat *A, params *P) {
  for (int row = P->S + P->H; row < P->L; row++) {
    uint32_t isi = (row - P->S) - P->H;
    uint16_vec idxs = params_get_idxs(isi, P);
    for (int idx = 0; idx < kv_size(idxs); idx++) {
      wrkmat_set(A, row, kv_A(idxs, idx), 1);
    }
    kv_destroy(idxs);
  }
}

wrkmat *precode_matrix_gen(params *P, int overhead) {
  wrkmat *A = wrkmat_new(P->L + overhead, P->L);

  precode_matrix_make_LDPC1(A, P->S, P->B);
  precode_matrix_make_identity(A, P->S, 0, P->B);
  precode_matrix_make_LDPC2(A, P->W, P->S, P->P);
  precode_matrix_make_G_ENC(A, P);
  precode_matrix_make_identity(A, P->H, P->S, P->L - P->H);

  octmat HDPC = precode_matrix_make_HDPC(P);
  wrkmat_assign_block(A, &HDPC, P->S, 0, P->H, P->Kprime + P->S);

  return A;
}

static void decode_patch(params *P, wrkmat *A, struct bitmask *mask,
                         repair_vec *repair_bin, uint16_t num_symbols,
                         uint16_t overhead) {
  size_t padding = P->Kprime - num_symbols;
  int num_gaps = bitmask_gaps(mask, num_symbols);
  int rep_idx = 0;
  for (int gap = 0; gap < P->L && num_gaps > 0; gap++) {
    if (bitmask_check(mask, gap))
      continue;
    int row = gap + P->H + P->S;
    for (int col = 0; col < A->cols; col++) {
      wrkmat_set(A, row, col, 0);
    }

    uint16_vec idxs =
        params_get_idxs(kv_A(*repair_bin, rep_idx++).esi + padding, P);
    for (int idx = 0; idx < kv_size(idxs); idx++) {
      wrkmat_set(A, row, kv_A(idxs, idx), 1);
    }
    kv_destroy(idxs);
    num_gaps--;
  }

  int rep_row = P->L;
  for (; rep_row < A->rows; rep_row++) {
    for (int col = 0; col < A->cols; col++) {
      wrkmat_set(A, rep_row, col, 0);
    }
    uint16_vec idxs =
        params_get_idxs(kv_A(*repair_bin, rep_idx++).esi + padding, P);
    for (int idx = 0; idx < kv_size(idxs); idx++) {
      wrkmat_set(A, rep_row, kv_A(idxs, idx), 1);
    }
    kv_destroy(idxs);
  }
}

static int decode_nz_fwd(wrkmat *A, int row, int s, int e, int c[], int d[]) {
  for (int col = s; col < e; col++) {
    int tmp = wrkmat_chk(A, d[row], c[col]);
    if (tmp > 0)
      return col;
  }
  return e;
}

static int decode_nz_rev(wrkmat *A, int row, int s, int e, int c[], int d[]) {
  for (int col = e - 1; col >= s; col--) {
    int tmp = wrkmat_chk(A, d[row], c[col]);
    if (tmp == 0)
      return col;
  }
  return s;
}

static void decode_rear_swap(wrkmat *A, int row, int s, int e, int c[],
                             int d[]) {
  while (s < e) {
    int swap1 = decode_nz_fwd(A, row, s, e, c, d);
    if (swap1 == e)
      return;
    int swap2 = decode_nz_rev(A, row, s, e, c, d);
    if (swap2 == s)
      return;
    if (swap1 >= swap2)
      return;
    TMPSWAP(int, c[swap1], c[swap2]);
    s++;
    e--;
  }
}

static bool decode_amd(params *P, wrkmat *A, schedule *S, int *vp) {
  int i = 0, u = P->P, *c = S->c, *d = S->d;
  int rows = A->rows, cols = A->cols;

  int counts[rows];
  for (int row = 0; row < rows; row++) {
    counts[d[row]] = wrkmat_nnz(A, d[row], 0, cols - u);
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
    }
    // find first one
    if (wrkmat_chk(A, d[V0], c[V0]) == 0) {
      int first = decode_nz_fwd(A, V0, V0 + 1, V0 + Vcols, c, d);
      TMPSWAP(int, c[V0], c[first]);
    }
    decode_rear_swap(A, V0, V0 + 1, V0 + Vcols, c, d);
    // decrement nz counts if row had nz at V0 or nz's at last r - 1 cols
    for (int row = V0 + 1; row < rows - P->H; row++) {
      for (int col = 0; col < (r - 1); col++) {
        counts[d[row]] -= !!wrkmat_chk(A, d[row], c[V0 + Vcols - col - 1]);
      }
      uint8_t beta = wrkmat_chk(A, d[row], c[V0]);
      if (beta) {
        counts[d[row]]--;
        wrkmat_axpy(A, d[row], d[V0], beta);
        sched_push(S, d[row], d[V0], beta);
      }
    }
    i++;
    u += r - 1;
  }
  *vp = i;
  return true;
}

static bool decode_hdpc(params *P, wrkmat *A, schedule *S, int row_start) {
  int rows = A->rows;
  int *c = S->c, *d = S->d;

  for (int row = 0; row < row_start; row++) {
    for (int h = P->S + P->Kprime; h < rows; h++) {
      uint8_t beta = wrkmat_at(A, d[h], c[row]);
      if (beta) {
        wrkmat_axpy(A, d[h], d[row], beta);
        sched_push(S, d[h], d[row], beta);
      }
    }
  }

  return true;
}

static bool decode_solve(params *P, wrkmat *A, schedule *S, int row_start) {
  int rows = A->rows, cols = A->cols;
  int *c = S->c, *d = S->d;
  uint8_t beta;

  for (int row = row_start; row < rows; row++) {
    int nzrow = row;

    // past diagonal
    if (row >= cols) {
      break;
    }

    for (; nzrow < rows; nzrow++) {
      beta = wrkmat_at(A, d[nzrow], c[row]);
      if (beta != 0) {
        break;
      }
    }

    if (nzrow == rows) {
      return false;
    }

    if (row != nzrow) {
      TMPSWAP(int, d[row], d[nzrow]);
    }

    beta = wrkmat_at(A, d[row], c[row]);
    if (beta > 1) {
      wrkmat_scal(A, d[row], OCT_INV[beta]);
      sched_push(S, d[row], OCT_INV[beta], 0);
    }

    for (int del_row = row + 1; del_row < rows; del_row++) {
      beta = wrkmat_at(A, d[del_row], c[row]);
      if (beta == 0)
        continue;
      wrkmat_axpy(A, d[del_row], d[row], beta);
      sched_push(S, d[del_row], d[row], beta);
    }
  }

  for (int row = P->L - 1; row >= row_start; row--) {
    for (int del_row = 0; del_row < row; del_row++) {
      beta = wrkmat_at(A, d[del_row], c[row]);
      if (beta == 0)
        continue;
      // wrkmat_axpy(A, d[del_row], d[row], beta);
      sched_push(S, d[del_row], d[row], beta);
    }
  }

  return true;
}

schedule *precode_matrix_invert(params *P, wrkmat *A) {
  int rows = A->rows, cols = A->cols, vp = 0;
  schedule *S = sched_new(rows, cols, P->L);

  // roll precode matrix to put G_ENC rows on top
  for (int row = 0; row < rows; row++) {
    S->d[row] = (row + P->S + P->H) % rows;
  }

  if (!decode_amd(P, A, S, &vp)) {
    wrkmat_free(A);
    sched_free(S);
    return NULL;
  }

  decode_hdpc(P, A, S, vp);

  if (!decode_solve(P, A, S, vp)) {
    wrkmat_free(A);
    sched_free(S);
    return NULL;
  }

  wrkmat_free(A);
  return S;
}

bool precode_matrix_intermediate1(params *P, wrkmat *A, octmat *D) {
  int rows = A->rows, cols = A->cols;
  if (rows <= 0)
    return false;

  schedule *S = precode_matrix_invert(P, A);
  if (S == NULL)
    return false;

  int dinv[rows];
  for (int i = 0; i < rows; i++)
    dinv[S->d[i]] = i;

  precode_matrix_apply_sched(D, S);
  precode_matrix_permute(D, dinv, rows);
  precode_matrix_permute(D, S->c, cols);

  sched_free(S);
  return true;
}

void precode_matrix_fill_slot(params *P, octmat *D, uint32_t isi, uint8_t *ptr,
                              size_t len) {
  uint16_vec idxs = params_get_idxs(isi, P);
  for (int idx = 0; idx < kv_size(idxs); idx++) {
    oaddrow(ptr, om_P(*D), 0, kv_A(idxs, idx), len);
  }
  kv_destroy(idxs);
}

bool precode_matrix_intermediate2(params *P, wrkmat *A, octmat *D, octmat *M,
                                  repair_vec *repair_bin,
                                  struct bitmask *repair_mask, int num_symbols,
                                  int overhead) {
  int num_gaps, gap = 0, row = 0;

  if (D->cols == 0) {
    return false;
  }

  decode_patch(P, A, repair_mask, repair_bin, num_symbols, overhead);

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
    uint16_t row = skip + gap;
    struct repair_sym rs = kv_A(*repair_bin, rep_idx++);
    ocopy(om_P(*D), om_P(rs.row), row, 0, D->cols);
  }

  for (int row = skip + P->Kprime; rep_idx < num_repair; row++) {
    struct repair_sym rs = kv_A(*repair_bin, rep_idx++);
    ocopy(om_P(*D), om_P(rs.row), row, 0, D->cols);
  }

  wrkmat *A = precode_matrix_gen(P, overhead);
  bool precode_ok = precode_matrix_intermediate2(
      P, A, D, M, repair_bin, repair_mask, num_symbols, overhead);
  return precode_ok;
}
