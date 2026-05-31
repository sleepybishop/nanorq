#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "gf2.h"
#include "wrkmat.h"

wrkmat *wrkmat_new(int rows, int cols) {
  wrkmat *w = calloc(1, sizeof(wrkmat));
  w->rows = rows;
  w->cols = cols;

  w->GF2 = gf2mat_new(rows, cols);
  w->rowmap = calloc(sizeof(int), rows);
  w->type = calloc(sizeof(int), rows);

  return w;
}

void wrkmat_free(wrkmat *w) {
  if (!w)
    return;
  if (w->rowmap)
    free(w->rowmap);
  if (w->type)
    free(w->type);
  if (w->GF2)
    gf2mat_free(w->GF2);
  om_destroy(&w->GF256);
  free(w);
}

void wrkmat_assign_block(wrkmat *w, octmat *B, int i, int j, int m, int n) {
  w->GF256 = *B;

  // overlay GF256 block over GF2
  for (int row = i; row < i + m; row++) {
    w->type[row] = 1;
    w->rowmap[row] = row - i;
  }
  w->blkidx = m;
}

void wrkmat_print(wrkmat *w, FILE *stream) {
  int m = w->rows, n = w->cols;
  fprintf(stream, "wrk [%ux%u]\n", m, n);
  fprintf(stream, "|     ");
  for (int j = 0; j < n; j++) {
    fprintf(stream, "| %03d ", j);
  }
  fprintf(stream, "|\n");
  for (int i = 0; i < m; i++) {
    fprintf(stream, "| %03d | %3d ", i, wrkmat_get(w, i, 0));
    for (int j = 1; j < n; j++) {
      fprintf(stream, "| %3d ", wrkmat_get(w, i, j));
    }
    fprintf(stream, "|\n");
  }
}

uint8_t wrkmat_get(wrkmat *w, int i, int j) {
  return w->type[i] ? om_A(w->GF256, w->rowmap[i], j)
                    : gf2mat_get(w->GF2, i, j);
}

void wrkmat_set(wrkmat *w, int i, int j, uint8_t b) {
  if (w->type[i]) {
    om_A(w->GF256, w->rowmap[i], j) = b;
  } else if (b <= 1) {
    gf2mat_set(w->GF2, i, j, b);
  } else {
    assert(0 && "unhandled set");
  }
}

void wrkmat_axpy(wrkmat *w, int i, int j, int beta) {
  if (w->type[i] == w->type[j]) {
    if (w->type[i]) {
      oaxpy(om_P(w->GF256), om_P(w->GF256), w->rowmap[i], w->rowmap[j], w->cols,
            beta);
    } else {
      gf2mat_xor(w->GF2, w->GF2, i, j);
    }
  } else {
    // if target row is in gf256, axpy in place from gf2 row
    if (w->type[i]) {
      // uint8_t *tmp = om_R(w->GF256, w->rowmap[i]);
      // gf2mat_axpy(w->GF2, j, tmp, beta);
      uint32_t *tmp = w->GF2->bits + w->GF2->stride * j;
      oaxpy_b32(om_P(w->GF256), tmp, w->rowmap[i], w->cols, beta);
    } else {
      if (w->blkidx >= w->GF256.rows) {
        int new_rows = w->GF256.rows * 2;
        if (new_rows == 0)
          new_rows = 16;
        octmat new_mat;
        om_resize(&new_mat, new_rows, w->GF256.cols);
        memcpy(new_mat.data, w->GF256.data, w->GF256.rows * w->GF256.cols_al);
        om_destroy(&w->GF256);
        w->GF256 = new_mat;
      }
      uint8_t *tmp = om_R(w->GF256, w->blkidx);
      memcpy(tmp, om_R(w->GF256, w->rowmap[j]), w->GF256.cols_al);
      if (beta != 1) {
        oscal(om_P(w->GF256), w->blkidx, w->cols, beta);
      }
      uint32_t *bit_row = w->GF2->bits + w->GF2->stride * i;
      oaxpy_b32(om_P(w->GF256), bit_row, w->blkidx, w->cols, 1);

      w->type[i] = 1; // row i is now a gf256 row
      w->rowmap[i] = w->blkidx;
      w->blkidx++;
    }
  }
}

void wrkmat_scal(wrkmat *w, int i, int beta) {
  if (w->type[i]) {
    oscal(om_P(w->GF256), w->rowmap[i], w->cols, beta);
  } else {
    assert(0 && "unhandled scal");
  }
}
