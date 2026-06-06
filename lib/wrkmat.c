#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "oblas_lite.h"
#include "util.h"
#include "wrkmat.h"

wrkmat *wrkmat_new(int rows, int cols) {
  wrkmat *w = calloc(1, sizeof(wrkmat));
  w->rows = rows;
  w->cols = cols;

  w->gf2_stride = ((cols / 32) + ((cols % 32) ? 1 : 0));
  w->gf2_data = calloc(rows, w->gf2_stride * sizeof(uint32_t));
  w->row_ptrs = calloc(rows, sizeof(void *));
  w->type = calloc(rows, sizeof(uint8_t));

  for (int i = 0; i < rows; i++) {
    w->row_ptrs[i] = w->gf2_data + i * w->gf2_stride * sizeof(uint32_t);
    w->type[i] = 0;
  }

  size_t pool_stride = (cols + nanorq_oblas.align_size - 1) & ~(nanorq_oblas.align_size - 1);
  if (pool_stride == 0) pool_stride = nanorq_oblas.align_size;
  w->pool_data = obl_alloc(rows, pool_stride, nanorq_oblas.align_size);
  w->pool_next = 0;

  return w;
}

void wrkmat_free(wrkmat *w) {
  if (!w)
    return;
  if (w->row_ptrs)
    free(w->row_ptrs);
  if (w->type)
    free(w->type);
  if (w->gf2_data)
    free(w->gf2_data);
  if (w->pool_data)
    obl_free(w->pool_data);
  if (w->GF256) {
    if (w->GF256[0]) obl_free(w->GF256[0]);
    free(w->GF256);
  }
  free(w);
}

void wrkmat_assign_block(wrkmat *w, uint8_t **B, int i, int j, int m, int n) {
  w->GF256 = B;

  // overlay GF256 block over GF2
  for (int row = i; row < i + m; row++) {
    w->type[row] = 1;
    w->row_ptrs[row] = w->GF256[row - i];
  }
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
  return wrkmat_at(w, i, j);
}

void wrkmat_set(wrkmat *w, int i, int j, uint8_t b) {
  if (w->type[i]) {
    ((uint8_t *)w->row_ptrs[i])[j] = b;
  } else if (b <= 1) {
    uint32_t *row = (uint32_t *)w->row_ptrs[i];
    uint32_t mask = 1U << ((unsigned)j % 32);
    if (b) {
      row[(unsigned)j / 32] |= mask;
    } else {
      row[(unsigned)j / 32] &= ~mask;
    }
  } else {
    assert(0 && "unhandled set");
  }
}

static void wrkmat_promote(wrkmat *w, int i) {
  if (w->type[i] == 1)
    return;

  size_t stride = (w->cols + nanorq_oblas.align_size - 1) & ~(nanorq_oblas.align_size - 1);
  if (stride == 0) stride = nanorq_oblas.align_size;
  assert(w->pool_next < w->rows);
  uint8_t *promoted_row = w->pool_data + w->pool_next * stride;
  w->pool_next++;

  memset(promoted_row, 0, stride);
  nanorq_oblas.axpyb32(promoted_row, (uint32_t *)w->row_ptrs[i], 1, w->cols);

  w->row_ptrs[i] = promoted_row;
  w->type[i] = 1;
}

void wrkmat_axpy(wrkmat *w, int i, int j, int beta) {
  if (w->type[i] == w->type[j]) {
    if (w->type[i]) {
      nanorq_oblas.axpy((uint8_t *)w->row_ptrs[i], (uint8_t *)w->row_ptrs[j], beta, w->cols);
    } else {
      uint32_t * restrict ap = (uint32_t *)w->row_ptrs[i];
      uint32_t * restrict bp = (uint32_t *)w->row_ptrs[j];
      size_t stride = w->gf2_stride;
      for (size_t idx = 0; idx < stride; idx++) {
        ap[idx] ^= bp[idx];
      }
    }
  } else {
    if (w->type[i]) {
      nanorq_oblas.axpyb32((uint8_t *)w->row_ptrs[i], (uint32_t *)w->row_ptrs[j], beta, w->cols);
    } else {
      wrkmat_promote(w, i);
      nanorq_oblas.axpy((uint8_t *)w->row_ptrs[i], (uint8_t *)w->row_ptrs[j], beta, w->cols);
    }
  }
}

void wrkmat_scal(wrkmat *w, int i, int beta) {
  if (w->type[i]) {
    nanorq_oblas.scal((uint8_t *)w->row_ptrs[i], beta, w->cols);
  } else {
    assert(0 && "unhandled scal");
  }
}
