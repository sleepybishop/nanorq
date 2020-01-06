#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <oblas.h>

#include "bitmask.h"
#include "chooser.h"
#include "graph.h"
#include "params.h"
#include "precode.h"
#include "rand.h"

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

static octmat precode_matrix_make_GAMMA(uint16_t dim) {
  octmat GAMMA = OM_INITIAL;
  om_resize(&GAMMA, dim, dim);

  for (int row = 0; row < GAMMA.rows; row++) {
    int col;
    for (col = 0; col <= row; col++)
      om_A(GAMMA, row, col) = OCT_EXP[(row - col) % OCT_EXP_SIZE];
    for (; col < GAMMA.cols; col++) {
      om_A(GAMMA, row, col) = 0;
    }
  }
  return GAMMA;
}

static void precode_matrix_init_HDPC(params *P, octmat *A) {
  uint16_t m = P->H;
  uint16_t n = P->Kprime + P->S;

  if (m == 0 || n == 0)
    return;

  octmat MT = precode_matrix_make_MT(m, n);
  octmat GAMMA = precode_matrix_make_GAMMA(n);
  octmat MTxGAMMA = OM_INITIAL;

  om_resize(&MTxGAMMA, MT.rows, GAMMA.cols);

  ogemm(om_P(MT), om_P(GAMMA), om_P(MTxGAMMA), MT.rows, MT.cols, GAMMA.cols);

  int row, col;
  for (col = 0; col < GAMMA.cols; col++) {
    for (row = 0; row < MT.rows; row++) {
      om_A(*A, P->S + row, col) = om_A(MTxGAMMA, row, col);
    }
  }
  om_destroy(&MT);
  om_destroy(&GAMMA);
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

static bool decode_amd(params *P, octmat *A, octmat *X, octmat *D,
                          int c[], int *i_val, int *u_val) {
  uint16_t i = 0;
  uint16_t u = P->P;

  struct chooser ch = chooser_init(A->rows);

  for (int row = 0; row < A->rows; row++) {
    bool is_hdpc = (row >= P->S && row < (P->S + P->H));
    size_t row_degree = 0;
    for (int col = 0; col < A->cols - u; col++) {
      row_degree += (uint8_t)(om_A(*A, row, col));
    }
    chooser_add_tracking_pair(&ch, is_hdpc, row_degree);
  }

  while (i + u < P->L) {
    uint16_t sub_rows = A->rows - i;
    uint16_t sub_cols = A->cols - i - u;
    uint16_t chosen, non_zero;
    struct graph *G = graph_new(sub_cols);

    non_zero = chooser_non_zero(&ch, A, G, i, sub_rows, sub_cols);
    if (non_zero == sub_cols + 1) {
      chooser_clear(&ch);
      graph_free(G);

      *i_val = 0;
      *u_val = 0;
      return false;
    }
    chosen = chooser_pick(&ch, G, i, sub_rows, non_zero);
    if (chosen != 0) {
      oswaprow(om_P(*A), i, chosen + i, A->cols);
      // oswaprow(om_P(*X), i, chosen + i, X->cols);
      oswaprow(om_P(*D), i, chosen + i, D->cols);

      kv_swap(struct tracking_pair, ch.tracking, i, chosen + i);
    }

    if (om_A(*A, i, i) == 0) {
      int idx = 1;
      for (; idx < sub_cols; idx++) {
        if (om_A(*A, i, idx + i) != 0)
          break;
      }

      oswapcol(om_P(*A), i, i + idx, A->rows, A->cols);
      // oswapcol(om_P(*X), i, i + idx, X->rows, X->cols);

      TMPSWAP(int, c[i], c[i + idx]);
    }
    int col = sub_cols - 1;
    int swap = 1;
    for (; col > sub_cols - non_zero; col--) {
      if (om_A(*A, i, col + i) != 0)
        continue;
      while (swap < col && om_A(*A, i, swap + i) == 0) {
        swap++;
      }

      if (swap >= col)
        break;

      oswapcol(om_P(*A), col + i, swap + i, A->rows, A->cols);
      // oswapcol(om_P(*X), col + i, swap + i, X->rows, X->cols);

      TMPSWAP(int, c[col + i], c[swap + i]);
    }

    for (int row = 1; row < sub_rows; row++) {
      if (om_A(*A, row + i, i) != 0) {
        uint8_t mnum = om_A(*A, row + i, i);
        uint8_t mden = om_A(*A, i, i);
        uint8_t multiple = (mnum > 0 && mden > 0) ? OCTET_DIV(mnum, mden) : 0;
        if (multiple == 0)
          continue;
        oaxpy(om_P(*A), om_P(*A), row + i, i, A->cols, multiple);
        oaxpy(om_P(*D), om_P(*D), row + i, i, D->cols, multiple);
      }
    }
    i++;
    u += non_zero - 1;

    graph_free(G);
  }
  chooser_clear(&ch);

  *i_val = i;
  *u_val = u;
  return true;
}

static bool decode_solve(params *P, octmat *A, octmat *D, uint16_t i,
                          uint16_t u) {

  uint16_t row_start = P->S; // start after ldcp rows
  int rows = A->rows;
  uint8_t multiple;

  // defer the hdpc rows
  for (int l = 0; l < P->H; l++) {
    oswaprow(om_P(*A), l + P->S, A->rows - P->H + l, A->cols);
    oswaprow(om_P(*D), l + P->S, D->rows - P->H + l, D->cols);
  }

  for (int row = row_start; row < rows; row++) {
    int nzrow = row;

    // past diagonal
    if (row >= A->cols) {
      break;
    }

    for (; nzrow < rows; nzrow++) {
      multiple = om_A(*A, nzrow, row);
      if (multiple != 0) {
        break;
      }
    }

    if (nzrow == rows) {
      abort();
      return false;
    }

    if (row != nzrow) {
      oswaprow(om_P(*A), row, nzrow, A->cols);
      oswaprow(om_P(*D), row, nzrow, D->cols);
    }

    multiple = om_A(*A, row, row);
    if (multiple > 1) {
      oscal(om_P(*A), row, A->cols, OCT_INV[multiple]);
      oscal(om_P(*D), row, D->cols, OCT_INV[multiple]);
    }

    for (int del_row = row + 1; del_row < rows; del_row++) {
      multiple = om_A(*A, del_row, row);
      if (multiple == 0)
        continue;
      oaxpy(om_P(*A), om_P(*A), del_row, row, A->cols, multiple);
      oaxpy(om_P(*D), om_P(*D), del_row, row, D->cols, multiple);
    }
  }

  for (int del_row = P->L - 1; del_row >= row_start; del_row--) {
    for (int row = row_start; row < del_row; row++) {
      multiple = om_A(*A, row, del_row);
      oaxpy(om_P(*D), om_P(*D), row, del_row, D->cols, multiple);
    }
  }
  for (int row = row_start; row < P->L - 1; row++) {
    ozero(om_P(*A), row, A->cols);
    om_A(*A, row, row) = 1;
  }

  for (int row = 0; row < i; row++) {
    for (int col = 0; col < u; col++) {
      uint8_t multiple = om_A(*A, row, col + i);
      if (multiple == 0)
        continue;
      oaxpy(om_P(*D), om_P(*D), row, i + col, D->cols, 1);
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

void precode_matrix_permute(octmat *D, int P[], int n) {
  for (int i = 0; i < n; i++) { 
    int at = i, mark = -1 ;
    while (P[at] >= 0) { 
      oswaprow(om_P(*D), i, P[at], D->cols);
      int tmp = P[at]; 
      P[at] = mark;
      at = tmp; 
    } 
  }
}

bool precode_matrix_intermediate1(params *P, octmat *A, octmat *D) {
  bool success;
  int i, u, c[P->L];

  octmat X = OM_INITIAL;

  for (int l = 0; l < P->L; l++) {
    c[l] = l;
  }

  if (P->L == 0 || A == NULL || A->rows == 0 || A->cols == 0) {
    om_destroy(A);
    return false;
  }

  success = decode_amd(P, A, &X, D, c, &i, &u);
  if (!success) {
    om_destroy(A);
    return false;
  }

  success = decode_solve(P, A, D, i, u);

  if (!success) {
    om_destroy(A);
    return false;
  }

  precode_matrix_permute(D, c, P->L);
  om_destroy(A);

  return true;
}

bool precode_matrix_intermediate2(octmat *M, octmat *A, octmat *D, params *P,
                                  repair_vec *repair_bin, struct bitmask *mask,
                                  uint16_t num_symbols, uint16_t overhead) {

  int num_gaps, gap = 0, row = 0;

  if (D->cols == 0) {
    return false;
  }

  decode_patch(P, A, mask, repair_bin, num_symbols, overhead);

  if(!precode_matrix_intermediate1(P, A, D)) {
    return false;
  }

  num_gaps = bitmask_gaps(mask, num_symbols);
  om_resize(M, num_gaps, D->cols);
  for (gap = 0; gap < num_symbols && num_gaps > 0; gap++) {
    if (bitmask_check(mask, gap))
      continue;
    octmat ret = precode_matrix_encode(P, D, gap);
    ocopy(om_P(*M), om_P(ret), row, 0, M->cols);
    om_destroy(&ret);
    row++;
    num_gaps--;
  }
  return true;
}

octmat precode_matrix_encode(params *P, octmat *C, uint32_t isi) {
  octmat ret = OM_INITIAL;
  om_resize(&ret, 1, C->cols);

  uint16_vec idxs = params_get_idxs(isi, P);
  for (int idx = 0; idx < kv_size(idxs); idx++) {
    oaddrow(om_P(ret), om_P(*C), 0, kv_A(idxs, idx), ret.cols);
  }
  kv_destroy(idxs);

  return ret;
}

bool precode_matrix_decode(params *P, octmat *X, repair_vec *repair_bin,
                           struct bitmask *mask) {
  uint16_t num_symbols = X->rows, rep_idx, num_gaps, num_repair, overhead;

  octmat A = OM_INITIAL;
  octmat D = OM_INITIAL;
  octmat M = OM_INITIAL;

  num_repair = kv_size(*repair_bin);
  num_gaps = bitmask_gaps(mask, num_symbols);

  if (num_gaps == 0)
    return true;

  if (num_repair < num_gaps)
    return false;

  overhead = num_repair - num_gaps;
  rep_idx = 0;
  precode_matrix_gen(P, &A, overhead);

  om_resize(&D, P->S + P->H + P->Kprime + overhead, X->cols);

  int skip = P->S + P->H;
  for (int row = 0; row < X->rows; row++) {
    ocopy(om_P(D), om_P(*X), skip + row, row, D.cols);
  }

  for (int gap = 0; gap < num_symbols && rep_idx < num_repair; gap++) {
    if (bitmask_check(mask, gap))
      continue;
    uint16_t row = skip + gap;
    struct repair_sym rs = kv_A(*repair_bin, rep_idx++);
    ocopy(om_P(D), om_P(rs.row), row, 0, D.cols);
  }

  for (int row = skip + P->Kprime; rep_idx < num_repair; row++) {
    struct repair_sym rs = kv_A(*repair_bin, rep_idx++);
    ocopy(om_P(D), om_P(rs.row), row, 0, D.cols);
  }

  bool precode_ok = precode_matrix_intermediate2(&M, &A, &D, P, repair_bin,
                                                 mask, num_symbols, overhead);
  om_destroy(&D);

  if (!precode_ok)
    return false;

  int miss_row = 0;
  for (int row = 0; row < num_symbols && miss_row < M.rows; row++) {
    if (bitmask_check(mask, row))
      continue;
    miss_row++;

    ocopy(om_P(*X), om_P(M), row, miss_row - 1, X->cols);
    bitmask_set(mask, row);
  }
  om_destroy(&M);
  return true;
}
