#ifndef WRKMAT_H
#define WRKMAT_H

#include <stdint.h>
#include <stdio.h>

typedef struct {
  void **row_ptrs;
  uint8_t *type;
  size_t rows;
  size_t cols;
  size_t gf2_stride;
  uint8_t *gf2_data;
  uint8_t **GF256;
  uint8_t *pool_data;
  size_t pool_next;
} wrkmat;

#define wrkmat_at(w, i, j)                                                     \
  ((w)->type[i] ? ((uint8_t *)(w)->row_ptrs[i])[j] : (uint8_t)((((uint32_t *)(w)->row_ptrs[i])[((unsigned)(j)) / 32] >> (((unsigned)(j)) % 32)) & 1))

wrkmat *wrkmat_new(int rows, int cols);
void wrkmat_free(wrkmat *w);
void wrkmat_assign_block(wrkmat *w, uint8_t **B, int i, int j, int m, int n);
void wrkmat_print(wrkmat *w, FILE *stream);

uint8_t wrkmat_get(wrkmat *w, int i, int j);
void wrkmat_set(wrkmat *w, int i, int j, uint8_t b);

void wrkmat_axpy(wrkmat *w, int i, int j, int beta);
void wrkmat_scal(wrkmat *w, int i, int beta);

#endif
