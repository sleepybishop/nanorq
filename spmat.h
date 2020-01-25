#ifndef SPMAT_H
#define SPMAT_H

#include <stdint.h>
#include <stdio.h>

#include "util.h"

typedef struct {
  size_t rows;
  size_t cols;
  int_vec *idxs;
} spmat;

spmat *spmat_new(int rows, int cols);
void spmat_free(spmat *s);

void spmat_clear_row(spmat *s, int i);
void spmat_push(spmat *s, int i, int j);
spmat *spmat_transpose(spmat *s);
int spmat_nnz(spmat *m, int row, int start, int end);

#endif
