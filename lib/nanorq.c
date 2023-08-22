#include "nanorq.h"
#include "precode.h"
#include "tuple.h"

/*
 * T: size of each symbol in bytes (should be aligned to Al)
 * K: number of symbols per block
 */
int nanorq_encoder_new(u32 K, u32 overhead, nanorq *rq)
{
    assert(K > 1 && K <= K_max);
    rq->P = params_init(K);
    rq->overhead = overhead;
    return 0;
}

size_t nanorq_calculate_prepare_memory(nanorq *rq)
{
    params *P = &rq->P;
    u32 mem = 0, snz = 3 * DC(P->B, P->S) + 3;
    u32 tmp[GENC_MAX];

    /* c/ci, d/di & nz/cnz */
    mem += 3 * PAD(sizeof(u32) * P->L);
    mem += 3 * PAD(sizeof(u32) * (P->L + rq->overhead));
    /* NZT (one empty vec) */
    mem += 3 * PAD(sizeof(u32_vec));
    mem += 2 * PAD(sizeof(u32) * (P->L + rq->overhead));
    /* A */
    mem += (P->L + rq->overhead) * sizeof(u32_vec);
    u32 memb4a = mem;
    mem += P->S * PAD(snz * sizeof(u32));
    mem += (P->Kprime + rq->overhead) * PAD(GENC_MAX * sizeof(u32));
    /* AT */
    mem += (P->L) * sizeof(u32_vec);
    mem += (mem - memb4a) + P->L * PAD(sizeof(u32));

    return (size_t)mem;
}

static void assign_prepare_memory(nanorq *rq, u8 *mem, size_t len)
{
    params *P = &rq->P;
    pc *W = &rq->W;
    u32 snz = 3 * DC(P->B, P->S) + 3;
    u8 *ptr = mem;

    W->rows = P->L + rq->overhead;
    W->cols = P->L;
    W->prep_mem.base = mem;
    W->prep_mem.used = 0;
    W->prep_mem.max = len;

    ptr += u32_vec_init(&W->c, ptr, W->cols, W->cols, 0);
    ptr += u32_vec_init(&W->ci, ptr, W->cols, W->cols, 0);
    ptr += u32_vec_init(&W->d, ptr, W->rows, W->rows, 0);
    ptr += u32_vec_init(&W->di, ptr, W->rows, W->rows, 0);

    ptr += u32_vec_init(&W->cnz, ptr, W->cols, W->cols, 0);
    ptr += u32_vec_init(&W->nz, ptr, W->rows, W->rows, 0);
    uv_zero(W->cnz);
    uv_zero(W->nz);

    W->NZT = (u32_vec *)(ptr);
    ptr += 3 * PAD(sizeof(u32_vec));
    for (u32 i = 1; i < 3; i++)
        ptr += u32_vec_init(&W->NZT[i], ptr, 0, W->rows, 0);

    W->A = (u32_vec *)(ptr);
    ptr += W->rows * sizeof(u32_vec);
    for (u32 i = 0; i < W->rows; i++)
        W->A[i].n = W->A[i].m = 0;
    for (u32 i = 0; i < P->S; i++)
        ptr += u32_vec_init(&W->A[i], ptr, 0, snz, 0);
    for (u32 i = P->S + P->H, esi = 0; i < W->rows; i++, esi++)
        ptr += u32_vec_init(&W->A[i], ptr, 0, GENC_MAX, 0);

    W->AT = (u32_vec *)(ptr);
    ptr += W->cols * sizeof(u32_vec);
    for (u32 i = 0; i < W->cols; i++)
        W->AT[i].n = W->AT[i].m = 0;

    W->prep_mem.used = ptr - mem;

    assert(W->prep_mem.used <= W->prep_mem.max);

    W->cb.on_choose_arg = 0x0;
    W->cb.on_choose = precode_matrix_choose;
    W->cb.on_op_arg = 0x0;
    W->cb.on_op = precode_matrix_on_op;
}

void nanorq_prepare(nanorq *rq, uint8_t *prep_mem, size_t pm_len)
{
    assign_prepare_memory(rq, prep_mem, pm_len);
    precode_matrix_gen(&rq->P, &rq->W);
    precode_matrix_prepare(&rq->P, &rq->W);
}

void nanorq_get_packet_mix(nanorq *rq, u32 esi, u32_vec *mix)
{
    params *P = &rq->P;
    u32 X = esi + (P->Kprime - P->K);
    uv_clear(*mix);
    params_set_idxs(P, X, mix);
}

void nanorq_replace_symbol(nanorq *rq, u32 row, u32 esi)
{
    params *P = &rq->P;
    pc *W = &rq->W;
    row += P->H + P->S;
    nanorq_get_packet_mix(rq, esi, &W->A[row]);
}

size_t nanorq_calculate_work_memory(nanorq *rq)
{
    params *P = &rq->P;
    pc *W = &rq->W;
    u32 mem = 0;
    /* U */
    mem += PAD(W->rows * sizeof(u32) * DC(W->u, 8 * sizeof(u32)));
    /* field maps */
    mem += 2 * PAD(sizeof(u32) * W->rows);
    /* UL */
    mem += PAD(2 * P->H * W->u);
    /* HDPC */
    mem += PAD(P->H * (P->Kprime + P->S));
    return (size_t)mem;
}

static void assign_work_memory(nanorq *rq, u8 *mem, size_t len)
{
    params *P = &rq->P;
    pc *W = &rq->W;
    u8 *ptr = mem;
    u32 tmp = 0;

    W->work_mem.base = mem;
    W->work_mem.used = 0;
    W->work_mem.max = len;

    tmp = W->rows * DC(W->u, 8 * sizeof(u32));
    ptr += u32_vec_init(&W->U, ptr, tmp, tmp, tmp / W->rows);
    uv_zero(W->U);

    ptr += u32_vec_init(&W->F.rowmap, ptr, W->rows, W->rows, 0);
    ptr += u32_vec_init(&W->F.type, ptr, W->rows, W->rows, 0);
    uv_zero(W->F.rowmap);
    uv_zero(W->F.type);

    tmp = 2 * P->H * W->u;
    ptr += u8_vec_init(&W->UL, ptr, tmp, tmp, W->u);
    uv_zero(W->UL);

    tmp = P->H * (P->Kprime + P->S);
    ptr += u8_vec_init(&W->HDPC, ptr, tmp, tmp, P->Kprime + P->S);

    W->work_mem.used = ptr - mem;
    assert(W->work_mem.used <= W->work_mem.max);
}

int nanorq_precalculate(nanorq *rq, u8 *work_mem, size_t wm_len)
{
    assign_work_memory(rq, work_mem, wm_len);
    return precode_matrix_invert(&rq->P, &rq->W);
}

void nanorq_set_op_callback(nanorq *rq, void *arg, void (*on_op)(void *, u32, u16, u8))
{
    pc *W = &rq->W;
    W->cb.on_op_arg = arg;
    W->cb.on_op = on_op;
}

void nanorq_set_choose_callback(nanorq *rq, void *arg, u32 (*on_choose)(void *, pc *, u32, u32, u32, u32))
{
    pc *W = &rq->W;
    W->cb.on_choose_arg = arg;
    W->cb.on_choose = on_choose;
}
