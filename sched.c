#include "sched.h"

schedule *sched_new(int rows, int cols, int estimated_ops) {
  schedule *S = calloc(1, sizeof(schedule));
  S->c = calloc(cols, sizeof(int));
  S->ci = calloc(cols, sizeof(int));
  S->d = calloc(rows, sizeof(int));
  S->di = calloc(rows, sizeof(int));

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
  kv_resize(struct sch_op, S->ops, estimated_ops);

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
  if (kv_max(S->ops))
    kv_destroy(S->ops);
  free(S);
}

void sched_push(schedule *S, int i, int j, int beta) {
  struct sch_op op = {.i = i, .j = j, .beta = beta};
  kv_push(struct sch_op, S->ops, op);
}
