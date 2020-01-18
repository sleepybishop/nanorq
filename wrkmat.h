#ifndef WRKMAT_H
#define WRKMAT_H

#include <stdint.h>
#include <stdio.h>

#include "oblas.h"
#include "octmat.h"

typedef struct {
  octmat A;
  int rows;
  int cols;
} wrkmat;

#define wrkmat_at(w, i, j) (om_A(w->A, i, j))

wrkmat *wrkmat_new(int rows, int cols);
void wrkmat_free(wrkmat *w);
void wrkmat_assign_block(wrkmat *A, octmat *B, int i, int j, int m, int n);
void wrkmat_print(wrkmat *w, FILE *stream);

uint8_t wrkmat_get(wrkmat *w, int i, int j);
void wrkmat_set(wrkmat *w, int i, int j, uint8_t b);

void wrkmat_axpy(wrkmat *w, int i, int j, int beta);
void wrkmat_scal(wrkmat *w, int i, int beta);

#endif
