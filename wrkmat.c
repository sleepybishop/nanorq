#include "wrkmat.h"
#include <stdlib.h>
#include <string.h>

wrkmat *wrkmat_new(int rows, int cols) {
  wrkmat *w = calloc(1, sizeof(wrkmat));
  w->rows = rows;
  w->cols = cols;

  om_resize(&w->A, rows, cols);

  return w;
}

void wrkmat_free(wrkmat *w) {
  if (!w)
    return;
  om_destroy(&w->A);
  free(w);
}

void wrkmat_assign_block(wrkmat *w, octmat *B, int i, int j, int m, int n) {
  for (int row = i; row < i + m; row++) {
    for (int col = j; col < j + n; col++) {
      om_A(w->A, row, col) = om_A(*B, row - i, col - j);
    }
  }
  om_destroy(B);
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

uint8_t wrkmat_get(wrkmat *w, int i, int j) { return om_A(w->A, i, j); }

void wrkmat_set(wrkmat *w, int i, int j, uint8_t b) { om_A(w->A, i, j) = b; }

void wrkmat_axpy(wrkmat *w, int i, int j, int beta) {
  oaxpy(om_P(w->A), om_P(w->A), i, j, w->cols, beta);
}

void wrkmat_scal(wrkmat *w, int i, int beta) {
  oscal(om_P(w->A), i, w->cols, beta);
}

int wrkmat_nnz(wrkmat *w, int i, int s, int e) {
  int nz = 0;
  for (int col = s; col < e; col++) {
    nz += (om_A(w->A, i, col) > 0);
  }
  return nz;
}
