#include <stdlib.h>
#include <string.h>

#include "spmat.h"

spmat *spmat_new(int rows, int cols) {
  spmat *s = calloc(1, sizeof(spmat));
  s->rows = rows;
  s->cols = cols;

  s->idxs = calloc(sizeof(int_vec), rows);

  for (int i = 0; i < rows; i++) {
    kv_init(s->idxs[i]);
    kv_resize(int, s->idxs[i], 10);
  }

  return s;
}

void spmat_free(spmat *s) {
  if (!s)
    return;
  if (s->idxs) {
    for (int i = 0; i < s->rows; i++)
      kv_destroy(s->idxs[i]);
    free(s->idxs);
  }
  free(s);
}

void spmat_clear_row(spmat *s, int i) { kv_size(s->idxs[i]) = 0; }

void spmat_push(spmat *s, int i, int j) { kv_push(int, s->idxs[i], j); }

spmat *spmat_transpose(spmat *s) {
  spmat *t = spmat_new(s->cols, s->rows);
  for (int i = 0; i < s->rows; i++) {
    for (int it = 0; it < kv_size(s->idxs[i]); it++) {
      spmat_push(t, kv_A(s->idxs[i], it), i);
    }
  }
  return t;
}

int spmat_nnz(spmat *s, int row, int start, int end) {
  int nz = 0;
  int_vec rs = s->idxs[row];
  for (int it = 0; it < kv_size(rs); it++) {
    int col = kv_A(rs, it);
    if (col >= start && col < end) {
      nz++;
    }
  }
  return nz;
}
