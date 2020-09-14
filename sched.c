#include "sched.h"

schedule *sched_new(int rows, int cols, int estimated_ops) {
  schedule *S = calloc(1, sizeof(schedule));
  S->rows = rows;
  S->cols = cols;
  S->c = calloc(cols, sizeof(int));
  S->ci = calloc(cols, sizeof(int));
  S->d = calloc(rows, sizeof(int));
  S->di = calloc(rows, sizeof(int));
  S->nz = calloc(rows, sizeof(int));
  // init permutation vectors
  for (int j = 0; j < cols; j++) {
    S->c[j] = j;
    S->ci[j] = j;
  }
  for (int i = 0; i < rows; i++) {
    S->d[i] = i;
    S->di[i] = i;
  }

  kv_init(S->ops);
  kv_resize(sched_op, S->ops, estimated_ops);

  return S;
}

void sched_free(schedule *S) {
  if (!S)
    return;
  if (S->c)
    free(S->c);
  if (S->d)
    free(S->d);
  if (S->ci)
    free(S->ci);
  if (S->di)
    free(S->di);
  if (S->nz)
    free(S->nz);
  if (kv_max(S->ops))
    kv_destroy(S->ops);
  free(S);
}

void sched_push(schedule *S, int i, int j, int beta) {
  sched_op op = {.i = i, .j = j, .beta = beta};
  kv_push(sched_op, S->ops, op);
}

void sched_rebuild_permutations(schedule *S) {
  for (int i = 0; i < S->rows; i++) {
    S->di[S->d[i]] = i;
  }
  for (int j = 0; j < S->cols; j++) {
    S->c[S->ci[j]] = j;
  }
}
