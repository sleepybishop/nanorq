#include <stdio.h>
#include <stdlib.h>

#include <nanorq.h>
#include "kvec.h"

typedef struct {
    u8 u;
    u32 i;
    u32 j;
} sched_op;

typedef kvec_t(sched_op) oplist;

typedef struct {
    u32 cp[2];
    u32 cpidx;
    oplist ops;
} schedule;

void ops_push(void *arg, u32 i, u16 j, u8 u);
void ops_run(nanorq *rq, u8_vec *D, schedule *S);
void ops_mix(nanorq *rq, u8_vec *D, u32 esi, u8 *ptr);
