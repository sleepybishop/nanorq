#include <stdlib.h>
#include <string.h>

#include "wrkmat.h"

wrkmat *wrkmat_new(unsigned rows, unsigned cols) {
  wrkmat *w = calloc(1, sizeof(wrkmat));
  w->rows = rows;
  w->cols = cols;

  w->GF2 = binmat_new(rows, cols);
  w->rowmap = calloc(sizeof(unsigned), rows);
  w->type = calloc(sizeof(unsigned), rows);

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
    binmat_free(w->GF2);
  if (w->GF256)
    octmat_free(w->GF256);
  free(w);
}

void wrkmat_assign_block(wrkmat *w, octmat *B, unsigned i, unsigned j,
                         unsigned m, unsigned n) {
  w->GF256 = B;
  for (int row = 0; row < m; row++) {
    for (int col = 0; col < n; col++)
      binmat_set(w->GF2, row + i, col + j, om_A(w->GF256, row, col) > 0);
    for (int col = n; col < w->cols; col++)
      om_A(w->GF256, row, col) = binmat_get(w->GF2, row + i, col);
  }

  // overlay GF256 block over GF2
  for (int row = i; row < i + m; row++) {
    w->type[row] = 1;
    w->rowmap[row] = row - i;
  }
  w->blkidx = m;
}

void wrkmat_print(wrkmat *w, FILE *stream) {
  unsigned m = w->rows, n = w->cols;
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

uint8_t wrkmat_get(wrkmat *w, unsigned i, unsigned j) {
  return w->type[i] ? om_A(w->GF256, w->rowmap[i], j)
                    : binmat_get(w->GF2, i, j);
}

void wrkmat_set(wrkmat *w, unsigned i, unsigned j, uint8_t b) {
  if (w->type[i]) {
    om_A(w->GF256, w->rowmap[i], j) = b;
  } else if (b <= 1) {
    binmat_set(w->GF2, i, j, b);
  } else {
    printf("%s:%d -- unhandled set\n", __FILE__, __LINE__);
    abort();
  }
}

void wrkmat_axpy(wrkmat *w, unsigned i, unsigned j, uint8_t b) {
  if (w->type[i] == w->type[j]) {
    if (w->type[i]) {
      oaxpy(w->GF256, w->GF256, w->rowmap[i], w->rowmap[j], b);
      binmat_add(w->GF2, w->GF2, i, j);
    } else {
      binmat_add(w->GF2, w->GF2, i, j);
    }
  } else {
    // if target row is in gf256, axpy in place from gf2 row
    if (w->type[i]) {
      uint32_t *tmp = w->GF2->bits + w->GF2->stride * j;
      oaxpy_b32(w->GF256, tmp, w->rowmap[i], b);
    } else {
      if (w->blkidx >= w->GF256->rows) {
        printf("%s:%d -- unhandled axpy into gf2 row %d from gf256 row %d with "
               "beta %d\n",
               __FILE__, __LINE__, i, j, b);
        abort();
      }
      uint8_t *tmp = om_R(w->GF256, w->blkidx);
      binmat_fill(w->GF2, i, tmp);
      w->type[i] = 1; // row i is now a gf256 row
      w->rowmap[i] = w->blkidx;
      w->blkidx++;
      oaxpy(w->GF256, w->GF256, w->rowmap[i], w->rowmap[j], b);
    }
  }
}

void wrkmat_scal(wrkmat *w, unsigned i, uint8_t b) {
  if (w->type[i]) {
    oscal(w->GF256, w->rowmap[i], b);
  } else {
    printf("%s:%d -- unhandled scal row %d by beta %d ,  %d\n", __FILE__,
           __LINE__, i, b, OCT_INV[b]);
    abort();
  }
}
