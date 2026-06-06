#include <stdbool.h>
#ifndef UVEC_H
#define UVEC_H

#include "uvec.h"

struct oblas_impl {
  void (*axpy)(uint8_t *restrict a, uint8_t *restrict b, uint8_t u, unsigned k);
  void (*scal)(uint8_t *restrict a, uint8_t u, unsigned k);
  void (*axiy)(uint8_t *restrict a, uint8_t *restrict b, uint8_t u, unsigned k);
  void (*axpyb32)(uint8_t *restrict a, uint32_t *restrict b, uint8_t u,
                  unsigned k);
  size_t align_size;
};
extern struct oblas_impl nanorq_oblas;

bool u8_vec_init(u8_vec *v, arena *a, u32 n, u32 m, u32 s) {
  v->a = alloc(a, sizeof(u8), nanorq_oblas.align_size, m);
  if (!v->a && m > 0)
    return false;
  v->s = s;
  v->n = n;
  v->m = m;
  return true;
}

bool u32_vec_init(u32_vec *v, arena *a, u32 n, u32 m, u32 s) {
  v->a = alloc(a, sizeof(u32), nanorq_oblas.align_size, m);
  if (!v->a && m > 0)
    return false;
  v->s = s;
  v->n = n;
  v->m = m;
  return true;
}

u8 bm_get(u32_vec *v, u32 i, u32 j) {
  u32 *a = v->a + i * v->s;
  u32 q = j / 32;
  u32 r = j % 32;
  return (a[q] >> r) & 1;
}

void bm_set(u32_vec *v, u32 i, u32 j) {
  u32 *a = v->a + i * v->s;
  u32 q = j / 32;
  u32 r = j % 32;
  a[q] ^= (1U << r);
}

void bm_add(u32_vec *v, u32 i, u32 j) {
  u32 s = v->s;
  u32 *restrict ap = v->a + i * s;
  u32 *restrict bp = v->a + j * s;
  for (u32 idx = 0; idx < s; idx++)
    ap[idx] ^= bp[idx];
}

void bm_fill(u32_vec *v, u32 i, u8 *dst) {
  u32 s = v->s;
  u32 *restrict ap = v->a + i * s;
  for (u32 idx = 0; idx < s; idx++) {
    u32 tmp = ap[idx];
    while (tmp > 0) {
      u32 tz = __builtin_ctz(tmp);
      tmp = tmp & (tmp - 1);
      dst[tz + idx * 32] = 1;
    }
  }
}

u32 bm_gap(u32_vec *v, u32 i, u32 until) {
  u32 s = v->s;
  u32 *restrict ap = v->a + i * s;
  for (u32 idx = 0; idx < DC(until, 32); idx++) {
    u32 tmp = ~ap[idx];
    while (tmp > 0) {
      u32 tz = __builtin_ctz(tmp);
      tmp = tmp & (tmp - 1);
      u32 at = tz + idx * 32;
      if (at < until)
        return at;
    }
  }
  return UINT32_MAX;
}

#endif
