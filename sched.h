#ifndef NANORQ_SCHED_H
#define NANORQ_SCHED_H

#include "util.h"

struct sch_op {
  uint8_t beta;
  uint16_t i;
  uint16_t j;
};

typedef kvec_t(struct sch_op) oplist;

typedef struct sch {
  int *c;
  int *d;
  int *ci;
  int *di;
  oplist ops;
} schedule;

schedule *sched_new(int rows, int cols, int estimated_ops);
void sched_free(schedule *S);
void sched_push(schedule *S, int i, int j, int beta);

#endif
