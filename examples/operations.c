#include "operations.h"

static void ops_apply(u8_vec *D, schedule *S, u32 idx)
{
    sched_op op = kv_A(S->ops, idx);
    if (op.u)
        obl_axpy(uv_R(*D, op.i), uv_R(*D, op.j), op.u, D->s);
    else
        obl_scal(uv_R(*D, op.i), op.j, D->s);
}

static void ops_apply_schedule(u8_vec *D, schedule *S)
{
    for (u32 idx = 0; idx < S->cp[1]; idx++)
        ops_apply(D, S, idx);
    ops_apply(D, S, S->cp[0]);
    for (u32 idx = S->cp[0] - 1; idx < S->cp[0]; idx--)
        ops_apply(D, S, idx);
    for (u32 idx = S->cp[1]; idx < kv_size(S->ops); idx++)
        ops_apply(D, S, idx);
    for (u32 idx = 0; idx <= S->cp[0]; idx++)
        ops_apply(D, S, idx);
}

static void ops_permute(u8_vec *D, u32 P[], u32 n)
{
    for (size_t i = 0; i < n; i++) {
        u32 at = i, mark = UINT32_MAX;
        while (P[at] < UINT32_MAX) {
            obl_swap(uv_R(*D, i), uv_R(*D, P[at]), D->s);
            u32 tmp = P[at];
            P[at] = mark;
            at = tmp;
        }
    }
}

void ops_push(void *arg, u32 i, u16 j, u8 u)
{
    schedule *S = (schedule *)arg;
    sched_op op = {.i = i, .j = j, .u = u};
    if (i == 0 && j == 0 && u == 0) {
        S->cp[S->cpidx++] = kv_size(S->ops) - 1;
    } else {
        kv_push(sched_op, S->ops, op);
    }
}

void ops_run(nanorq *rq, u8_vec *D, schedule *S)
{
    ops_apply_schedule(D, S);
    u32 rows = nanorq_get_pc_rows(rq);
    u32 cols = nanorq_get_pc_cols(rq);
    u32 rm_base[rows];
    u32 cm_base[cols];
    u32_vec rpv, cpv;
    u32_vec_init(&rpv, (u8 *)rm_base, rows, rows, 0);
    u32_vec_init(&cpv, (u8 *)cm_base, cols, cols, 0);
    nanorq_clone_pc_rows_pv(rq, &rpv);
    nanorq_clone_pc_cols_pv(rq, &cpv);
    ops_permute(D, uv_R(rpv, 0), rows);
    ops_permute(D, uv_R(cpv, 0), cols);
}

void ops_mix(nanorq *rq, u8_vec *D, u32 esi, u8 *ptr)
{
    u32 base[GENC_MAX];
    u32_vec mix;
    u32_vec_init(&mix, (u8 *)base, 0, GENC_MAX, 0);

    nanorq_get_packet_mix(rq, esi, &mix);
    for (u32 i = 0; i < D->s; i++)
        ptr[i] = 0;
    for (u32 it = 0; it < uv_size(mix); it++) {
        u32 row = uv_A(mix, it);
        obl_axpy(ptr, uv_R(*D, row), 1, D->s);
    }
}
