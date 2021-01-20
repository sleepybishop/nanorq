#ifndef SPMAT_H
#define SPMAT_H

#include <stdint.h>
#include <stdio.h>

#include "util.h"

typedef struct {
  unsigned rows;
  unsigned cols;
  uint_vec *idxs;
} spmat;

spmat *spmat_new(unsigned rows, unsigned cols);
void spmat_free(spmat *s);

void spmat_clear_row(spmat *s, unsigned i);
void spmat_push(spmat *s, unsigned i, unsigned j);
spmat *spmat_transpose(spmat *s);
unsigned spmat_nnz(spmat *s, unsigned row, unsigned start, unsigned end);

#endif
