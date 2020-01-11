#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <oblas.h>

#include "bitmask.h"
#include "params.h"
#include "precode.h"
#include "rand.h"

static void sched_push(schedule *S, int type, int i, int j, int beta) {
  struct sch_op op = {.type = type, .i = i, .j = j, .beta = beta};
  kv_push(struct sch_op, *S, op);
}

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
  for (int i = 0; i < kv_size(*S); i++) {
    struct sch_op op = kv_A(*S, i);
    switch (op.type) {
    case OP_SWAP:
      oswaprow(om_P(*D), op.i, op.j, D->cols);
      break;
    case OP_SCAL:
      oscal(om_P(*D), op.i, D->cols, op.beta);
      break;
    case OP_AXPY:
      oaxpy(om_P(*D), om_P(*D), op.i, op.j, D->cols, op.beta);
      break;
    }
  }
}

static void precode_matrix_init_LDPC1(octmat *A, uint16_t S, uint16_t B) {
  int row, col;
  for (row = 0; row < S; row++) {
    for (col = 0; col < B; col++) {
      uint16_t submtx = col / S;
      if ((row == (col % S)) || (row == (col + submtx + 1) % S) ||
          (row == (col + 2 * (submtx + 1)) % S)) {
        om_A((*A), row, col) = 1;
      }
    }
  }
}

static void precode_matrix_init_LDPC2(octmat *A, uint16_t skip, uint16_t rows,
                                      uint16_t cols) {
  for (int row = 0; row < rows; row++) {
    uint16_t start = row % cols;
    for (int col = 0; col < cols; col++) {
      uint8_t val = (col == start || col == (start + 1) % cols) > 0;
      om_A((*A), row, skip + col) = val;
    }
  }
}

static void precode_matrix_add_identity(octmat *A, uint16_t size,
                                        uint16_t skip_row, uint16_t skip_col) {
  for (int diag = 0; diag < size; diag++) {
    om_A((*A), skip_row + diag, skip_col + diag) = 1;
  }
}

static octmat precode_matrix_make_MT(uint16_t rows, uint16_t cols) {
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

static void precode_matrix_init_HDPC(params *P, octmat *A) {
  uint16_t m = P->H;
  uint16_t n = P->Kprime + P->S;

  if (m == 0 || n == 0)
    return;

  octmat MT = precode_matrix_make_MT(m, n);
  octmat MTxGAMMA = OM_INITIAL;
  om_resize(&MTxGAMMA, m, n);

  uint8_t gv[n];
  for (int col = 0; col < n; col++) {
    gv[col] = OCT_EXP[col % OCT_EXP_SIZE];
  }

  uint8_t *ap, *cp = om_P(MTxGAMMA);
  for (int row = 0; row < m; row++, cp += MT.cols_al) {
    ap = om_P(MT) + (row * MT.cols_al);
    for (int idx = 0; idx < n; idx++) {
      uint8_t tmp[n];
      for (int col = 0; col < n; col++)
        tmp[col] = (col > idx) ? 0 : gv[idx - col];
      oaxpy(cp, tmp, 0, 0, n, ap[idx]);
    }
  }

  int row, col;
  for (col = 0; col < n; col++) {
    for (row = 0; row < m; row++) {
      om_A(*A, P->S + row, col) = om_A(MTxGAMMA, row, col);
    }
  }
  om_destroy(&MT);
  om_destroy(&MTxGAMMA);
}

static void precode_matrix_add_G_ENC(params *P, octmat *A) {
  for (int row = P->S + P->H; row < P->L; row++) {
    uint32_t isi = (row - P->S) - P->H;
    uint16_vec idxs = params_get_idxs(isi, P);
    for (int idx = 0; idx < kv_size(idxs); idx++) {
      om_A(*A, row, kv_A(idxs, idx)) = 1;
    }
    kv_destroy(idxs);
  }
}

static void decode_patch(params *P, octmat *A, struct bitmask *mask,
                         repair_vec *repair_bin, uint16_t num_symbols,
                         uint16_t overhead) {

  size_t padding = P->Kprime - num_symbols;
  uint16_t num_gaps = bitmask_gaps(mask, num_symbols);
  uint16_t rep_idx = 0;
  for (int gap = 0; gap < P->L && num_gaps > 0; gap++) {
    if (bitmask_check(mask, gap))
      continue;
    uint16_t row = gap + P->H + P->S;
    for (int col = 0; col < A->cols; col++) {
      om_A(*A, row, col) = 0;
    }

    uint16_vec idxs =
        params_get_idxs(kv_A(*repair_bin, rep_idx++).esi + padding, P);
    for (int idx = 0; idx < kv_size(idxs); idx++) {
      om_A(*A, row, kv_A(idxs, idx)) = 1;
    }
    kv_destroy(idxs);
    num_gaps--;
  }

  int rep_row = (uint16_t)(A->rows - overhead);
  for (; rep_row < A->rows; rep_row++) {
    for (int col = 0; col < A->cols; col++) {
      om_A(*A, rep_row, col) = 0;
    }
    uint16_vec idxs =
        params_get_idxs(kv_A(*repair_bin, rep_idx++).esi + padding, P);
    for (int idx = 0; idx < kv_size(idxs); idx++) {
      om_A(*A, rep_row, kv_A(idxs, idx)) = 1;
    }
    kv_destroy(idxs);
  }
}

static bool decode_amd(params *P, octmat *A, octmat *D, int c[], schedule *S) {
  int i = 0;
  int u = P->P;

  // defer the hdpc rows
  for (int l = 0; l < P->H; l++) {
    oswaprow(om_P(*A), l + P->S, A->rows - P->H + l, A->cols);
    oswaprow(om_P(*D), l + P->S, A->rows - P->H + l, D->cols);
  }

  while (i + u < P->L) {
    int Vrows = A->rows - i, Vcols = A->cols - i - u, V0 = i;
    int nnz = Vcols + 1;
    int first_one = 0;
    // find min nnz
    for (int row = 0; row < Vrows; row++) {
      int vnnz = 0;
      for (int col = V0; col < V0 + Vcols; col++) {
        vnnz += (om_A(*A, V0 + row, c[col]) > 0);
      }
      if (vnnz == 0)
        continue;
      if (vnnz > nnz)
        break;
      nnz = vnnz;
    }

    if (nnz == Vcols + 1) {
      return false;
    }

    // skipping row swaps based on # of ones/graph etc.
    /*
    int chosen = 0;
    if (chosen != 0) {
      oswaprow(om_P(*A), V0, V0 + chosen, A->cols);
      sched_push(S, OP_SWAP, V0, V0 + chosen, 0);
    }
    */

    // find first 1
    for (int col = V0; col < V0 + Vcols; col++) {
      if (om_A(*A, V0, c[col]) > 0) {
        first_one = col - V0;
        break;
      }
    }

    // find the first column in V with a nonzero and move to V[0]
    if (first_one > 0) {
      int col = first_one;
      TMPSWAP(int, c[V0], c[V0 + col]);
    }

    // move all columns with non zeros to the back
    for (int col = Vcols - 1, swap = 1; col > Vcols - nnz; col--) {
      if (om_A(*A, V0, c[V0 + col]) != 0)
        continue;
      while (swap < col && om_A(*A, V0, c[V0 + swap]) == 0) {
        swap++;
      }

      if (swap >= col)
        break;

      TMPSWAP(int, c[V0 + col], c[V0 + swap]);
    }

    // cancel out non zeros in rows below V[0]
    for (int row = 1; row < Vrows; row++) {
      if (om_A(*A, V0 + row, c[V0]) != 0) {
        uint8_t beta = om_A(*A, V0 + row, c[V0]); // assuming V[0][0] is 1
        oaxpy(om_P(*A), om_P(*A), V0 + row, V0, A->cols, beta);
        sched_push(S, OP_AXPY, V0 + row, V0, beta);
      }
    }
    i++;
    u += nnz - 1;
  }

  return true;
}

static bool decode_solve(params *P, octmat *A, octmat *D, int c[],
                         schedule *S) {
  int row_start = P->S; // start after ldpc rows
  int rows = A->rows;
  uint8_t beta;

  for (int row = row_start; row < rows; row++) {
    int nzrow = row;

    // past diagonal
    if (row >= A->cols) {
      break;
    }

    for (; nzrow < rows; nzrow++) {
      beta = om_A(*A, nzrow, c[row]);
      if (beta != 0) {
        break;
      }
    }

    if (nzrow == rows) {
      return false;
    }

    if (row != nzrow) {
      oswaprow(om_P(*A), row, nzrow, A->cols);
      sched_push(S, OP_SWAP, row, nzrow, 0);
    }

    beta = om_A(*A, row, c[row]);
    if (beta > 1) {
      oscal(om_P(*A), row, A->cols, OCT_INV[beta]);
      sched_push(S, OP_SCAL, row, 0, OCT_INV[beta]);
    }

    for (int del_row = row + 1; del_row < rows; del_row++) {
      beta = om_A(*A, del_row, c[row]);
      if (beta == 0)
        continue;
      oaxpy(om_P(*A), om_P(*A), del_row, row, A->cols, beta);
      sched_push(S, OP_AXPY, del_row, row, beta);
    }
  }

  for (int row = P->L - 1; row >= 0; row--) {
    for (int del_row = 0; del_row < row; del_row++) {
      beta = om_A(*A, del_row, c[row]);
      if (beta == 0)
        continue;
      // oaxpy(om_P(*A), om_P(*A), del_row, row, A->cols, beta);
      sched_push(S, OP_AXPY, del_row, row, beta);
    }
  }

  return true;
}

void precode_matrix_gen(params *P, octmat *A, uint16_t overhead) {
  om_resize(A, P->L + overhead, P->L);

  precode_matrix_init_LDPC1(A, P->S, P->B);
  precode_matrix_add_identity(A, P->S, 0, P->B);
  precode_matrix_init_LDPC2(A, P->W, P->S, P->P);
  precode_matrix_init_HDPC(P, A);
  precode_matrix_add_identity(A, P->H, P->S, P->L - P->H);
  precode_matrix_add_G_ENC(P, A);
}

bool precode_matrix_intermediate1(params *P, octmat *A, octmat *D) {
  bool success;
  int c[P->L];
  schedule S;

  kv_init(S);
  kv_resize(struct sch_op, S, (P->L * P->L) / 2);
  for (int l = 0; l < P->L; l++) {
    c[l] = l;
  }

  if (P->L == 0 || A == NULL || A->rows == 0 || A->cols == 0) {
    om_destroy(A);
    return false;
  }

  success = decode_amd(P, A, D, c, &S);
  if (!success) {
    om_destroy(A);
    return false;
  }

  success = decode_solve(P, A, D, c, &S);
  if (!success) {
    om_destroy(A);
    return false;
  }
  om_destroy(A);

  precode_matrix_apply_sched(D, &S);
  precode_matrix_permute(D, c, P->L);

  kv_destroy(S);
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

bool precode_matrix_intermediate2(params *P, octmat *A, octmat *D, octmat *M,
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
  octmat A = OM_INITIAL;

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

  precode_matrix_gen(P, &A, overhead);
  bool precode_ok = precode_matrix_intermediate2(
      P, &A, D, M, repair_bin, repair_mask, num_symbols, overhead);
  return precode_ok;
}
