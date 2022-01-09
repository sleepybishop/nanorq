#ifndef NANORQ_H
#define NANORQ_H

#include "params.h"
#include "precode.h"
#include <stddef.h>

typedef struct {
    pc W;
    params P;
    u32 overhead;
} nanorq;

/* returns a new encoder configured with given parameters */
int nanorq_encoder_new(uint32_t K, uint32_t overhead, nanorq *rq);

/* return the amount of memory required for preconditioning */
size_t nanorq_calculate_prepare_memory(nanorq *rq);

/* prepare rq matrix for inversion */
int nanorq_prepare(nanorq *rq, uint8_t *prep_mem, size_t pm_len);

/* replace symbols in rq matrix for given row */
int nanorq_replace_symbol(nanorq *rq, u32 row, u32 esi);

/* return the amount of memory required for inversion */
size_t nanorq_calculate_work_memory(nanorq *rq);

/* precalculate rq matrix inversion */
int nanorq_precalculate(nanorq *rq, uint8_t *work_mem, size_t wm_len);

/* set callback for choosing a row during preconditioning */
void nanorq_set_choose_callback(nanorq *rq, void *arg, u32 (*on_choose)(void *, pc *, u32, u32, u32, u32));

/* set callback for when data matrix operations are computed */
void nanorq_set_op_callback(nanorq *rq, void *arg, void (*on_op)(void *, u32, u16, u8));

#endif
