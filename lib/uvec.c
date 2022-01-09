#ifndef UVEC_H
#define UVEC_H

#include "uvec.h"

u32 u8_vec_init(u8_vec *v, u8 *a, u32 n, u32 m, u32 s)
{
    v->a = a;
    v->s = s;
    v->n = n;
    v->m = m;
    return PAD(m * sizeof(u8));
}

u32 u32_vec_init(u32_vec *v, u8 *a, u32 n, u32 m, u32 s)
{
    v->a = (u32 *)a;
    v->s = s;
    v->n = n;
    v->m = m;
    return PAD(m * sizeof(u32));
}

u8 bm_get(u32_vec *m, u32 i, u32 j)
{
    u32 *a = m->a + i * m->s;
    u32 q = j / 32;
    u32 r = j % 32;
    return (a[q] >> r) & 1;
}

void bm_set(u32_vec *m, u32 i, u32 j)
{
    u32 *a = m->a + i * m->s;
    u32 q = j / 32;
    u32 r = j % 32;
    a[q] = (a[q] & ~(1U << r)) | (1U << r);
}

void bm_add(u32_vec *v, u32 i, u32 j)
{
    u32 *ap = v->a + i * v->s;
    u32 *bp = v->a + j * v->s;
    for (u32 idx = 0; idx < v->s; idx++)
        ap[idx] ^= bp[idx];
}

void bm_fill(u32_vec *m, u32 i, u8 *dst)
{
    u32 *a = m->a + i * m->s;
    for (u32 idx = 0; idx < m->s; idx++) {
        u32 tmp = a[idx];
        while (tmp > 0) {
            u32 tz = __builtin_ctz(tmp);
            u32 q = tz;
            tmp = tmp & (tmp - 1);
            dst[q + idx * 32] |= 1U << (tz % 32);
        }
    }
}

#endif
