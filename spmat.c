#include <stdlib.h>
#include <string.h>

#include "spmat.h"

spmat *spmat_new(unsigned rows, unsigned cols) {
  spmat *s = calloc(1, sizeof(spmat));
  s->rows = rows;
  s->cols = cols;

  s->idxs = calloc(sizeof(uint_vec), rows);

  for (unsigned i = 0; i < rows; i++) {
    kv_init(s->idxs[i]);
    kv_resize(unsigned, s->idxs[i], 10);
  }

  return s;
}

void spmat_free(spmat *s) {
  if (!s)
    return;
  if (s->idxs) {
    for (unsigned i = 0; i < s->rows; i++)
      kv_destroy(s->idxs[i]);
    free(s->idxs);
  }
  free(s);
}

void spmat_clear_row(spmat *s, unsigned i) { kv_size(s->idxs[i]) = 0; }

void spmat_push(spmat *s, unsigned i, unsigned j) {
  kv_push(unsigned, s->idxs[i], j);
}

spmat *spmat_transpose(spmat *s) {
  spmat *t = spmat_new(s->cols, s->rows);
  for (unsigned i = 0; i < s->rows; i++) {
    for (unsigned it = 0; it < kv_size(s->idxs[i]); it++) {
      spmat_push(t, kv_A(s->idxs[i], it), i);
    }
  }
  return t;
}

unsigned spmat_nnz(spmat *s, unsigned row, unsigned start, unsigned end) {
  unsigned nz = 0;
  uint_vec rs = s->idxs[row];
  for (unsigned it = 0; it < kv_size(rs); it++) {
    unsigned col = kv_A(rs, it);
    if (col >= start && col < end) {
      nz++;
    }
  }
  return nz;
}
