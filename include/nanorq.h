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
void nanorq_prepare(nanorq *rq, uint8_t *prep_mem, size_t pm_len);

/* fetch the codebook for the given packet */
void nanorq_get_packet_mix(nanorq *rq, u32 esi, u32_vec *mix);

/* replace symbols in rq matrix for given row */
void nanorq_replace_symbol(nanorq *rq, u32 row, u32 esi);

/* return the amount of memory required for inversion */
size_t nanorq_calculate_work_memory(nanorq *rq);

/* precalculate rq matrix inversion */
int nanorq_precalculate(nanorq *rq, uint8_t *work_mem, size_t wm_len);

/* set callback for choosing a row during preconditioning */
void nanorq_set_choose_callback(nanorq *rq, void *arg, u32 (*on_choose)(void *, pc *, u32, u32, u32, u32));

/* set callback for when data matrix operations are computed */
void nanorq_set_op_callback(nanorq *rq, void *arg, void (*on_op)(void *, u32, u16, u8));

/* get the offset to the GENC rows in precode matrix */
static uint32_t nanorq_get_pc_genc_offset(nanorq *rq);

/* get the number of rows for the precode matrix */
static uint32_t nanorq_get_pc_rows(nanorq *rq);

/* get the number of columns for the precode matrix */
static uint32_t nanorq_get_pc_cols(nanorq *rq);

/* clone the precode matrix row pivot vector to supplied array */
static void nanorq_clone_pc_rows_pv(nanorq *rq, u32_vec *pv);

/* clone the precode matrix column pivot vector to supplied array */
static void nanorq_clone_pc_cols_pv(nanorq *rq, u32_vec *pv);

uint32_t nanorq_get_pc_genc_offset(nanorq *rq)
{
    return rq->P.S + rq->P.H;
}

inline uint32_t nanorq_get_pc_rows(nanorq *rq)
{
    return rq->W.rows;
}
inline uint32_t nanorq_get_pc_cols(nanorq *rq)
{
    return rq->W.cols;
}
inline void nanorq_clone_pc_rows_pv(nanorq *rq, u32_vec *pv)
{
    for (u32 i = 0; i < rq->W.rows; i++)
        uv_A(*pv, i) = uv_A(rq->W.di, i);
}
inline void nanorq_clone_pc_cols_pv(nanorq *rq, u32_vec *pv)
{
    for (u32 i = 0; i < rq->W.cols; i++)
        uv_A(*pv, i) = uv_A(rq->W.c, i);
}

#endif
