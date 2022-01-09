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
  unsigned rows;
  unsigned cols;
  int *c;  /* column permutation */
  int *d;  /* row permutation */
  int *ci; /* inverse map of c */
  int *di; /* inverse map of d */
  unsigned *nz;
  oplist ops; /* list of decoding operations */
  unsigned i; /* dim of X submatrix */
  unsigned u; /* remaining cols */

  unsigned marks[2]; /* checkpoints */
} schedule;

schedule *sched_new(unsigned rows, unsigned cols, unsigned estimated_ops);
void sched_free(schedule *S);
void sched_push(schedule *S, unsigned i, unsigned j, uint8_t beta);

#endif
