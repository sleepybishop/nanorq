#ifndef WRKMAT_H
#define WRKMAT_H

#include <stdint.h>
#include <stdio.h>

#include "oblas.h"
#include "octmat.h"

typedef struct {
  octmat A;
  size_t rows;
  size_t cols;
} wrkmat;

#define wrkmat_at(w, i, j) (om_A(w->A, i, j))
#define wrkmat_chk(w, i, j) (!!om_A(w->A, i, j))

wrkmat *wrkmat_new(size_t rows, size_t cols);
void wrkmat_free(wrkmat *w);
void wrkmat_assign_block(wrkmat *A, octmat *B, int i, int j, int m, int n);
void wrkmat_print(wrkmat *w, FILE *stream);

uint8_t wrkmat_get(wrkmat *w, int i, int j);
void wrkmat_set(wrkmat *w, int i, int j, uint8_t b);

void wrkmat_axpy(wrkmat *w, int i, int j, int beta);
void wrkmat_scal(wrkmat *w, int i, int beta);
int wrkmat_nnz(wrkmat *w, int i, int s, int e);

#endif
