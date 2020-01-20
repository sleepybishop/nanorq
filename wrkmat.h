#ifndef WRKMAT_H
#define WRKMAT_H

#include <stdint.h>
#include <stdio.h>

#include "gf2.h"
#include "oblas.h"
#include "octmat.h"

typedef struct {
  gf2mat *GF2;
  octmat GF256;
  size_t rows;
  size_t cols;
  size_t blkidx;
  int *rowmap;
  int *type;
} wrkmat;

#define wrkmat_chk(w, i, j) (gf2_at(w->GF2, i, j))
#define wrkmat_at(w, i, j)                                                     \
  (w->type[i] ? om_A(w->GF256, w->rowmap[i], j) : gf2_at(w->GF2, i, j))
#define wrkmat_nnz(w, i, s, e) gf2mat_nnz(w->GF2, i, s, e)

wrkmat *wrkmat_new(int rows, int cols);
void wrkmat_free(wrkmat *w);
void wrkmat_assign_block(wrkmat *w, octmat *B, int i, int j, int m, int n);
void wrkmat_print(wrkmat *w, FILE *stream);

uint8_t wrkmat_get(wrkmat *w, int i, int j);
void wrkmat_set(wrkmat *w, int i, int j, uint8_t b);

void wrkmat_axpy(wrkmat *w, int i, int j, int beta);
void wrkmat_scal(wrkmat *w, int i, int beta);

#endif
