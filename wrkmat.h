#ifndef WRKMAT_H
#define WRKMAT_H

#include <stdint.h>
#include <stdio.h>

#include "binmat.h"
#include "oblas.h"
#include "octmat.h"

typedef struct {
  binmat *GF2;
  octmat *GF256;
  unsigned rows;
  unsigned cols;
  unsigned blkidx;
  unsigned *rowmap;
  unsigned *type;
} wrkmat;

#define GF2ROW(m, i) (m->bits + i * m->stride)
#define GF2AT(m, i, j) ((GF2ROW(m, i)[(j) / 32] >> ((j) % 32)) & 1)

#define wrkmat_at(w, i, j)                                                     \
  (w->type[i] ? om_A(w->GF256, w->rowmap[i], j) : GF2AT(w->GF2, i, j))

wrkmat *wrkmat_new(unsigned rows, unsigned cols);
void wrkmat_free(wrkmat *w);
void wrkmat_assign_block(wrkmat *w, octmat *B, unsigned i, unsigned j,
                         unsigned m, unsigned n);
void wrkmat_print(wrkmat *w, FILE *stream);

uint8_t wrkmat_get(wrkmat *w, unsigned i, unsigned j);
void wrkmat_set(wrkmat *w, unsigned i, unsigned j, uint8_t b);

void wrkmat_axpy(wrkmat *w, unsigned i, unsigned j, uint8_t b);
void wrkmat_scal(wrkmat *w, unsigned i, uint8_t b);

#endif
