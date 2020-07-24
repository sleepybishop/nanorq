#ifndef NANORQ_SCHED_H
#define NANORQ_SCHED_H

#include "util.h"

typedef struct {
  uint8_t beta;
  uint32_t i;
  uint32_t j;
} sched_op;

typedef kvec_t(sched_op) oplist;

typedef struct {
  int *c;  /* column permutation */
  int *d;  /* row permutation */
  int *ci; /* inverse map of c */
  int *di; /* inverse map of d */
  int *nz;
  oplist ops; /* list of decoding operations */
  int i;      /* dim of X submatrix */
  int u;      /* remaining cols */

  int Xs; /* start of operations on X */
  int Xe; /* end of operations on X */
} schedule;

schedule *sched_new(int rows, int cols, int estimated_ops);
void sched_free(schedule *S);
void sched_push(schedule *S, int i, int j, int beta);

#endif
