#ifndef NANORQ_UTIL_H
#define NANORQ_UTIL_H

#include <assert.h>
#include <stdint.h>
#include <string.h>

#include <octmat.h>

#include "kvec.h"

#define div_ceil(A, B) ((A) / (B) + ((A) % (B) ? 1 : 0))
#define div_floor(A, B) ((A) / (B))

#define TMPSWAP(type, a, b)                                                    \
  do {                                                                         \
    type __tmp = a;                                                            \
    a = b;                                                                     \
    b = __tmp;                                                                 \
  } while (0)

typedef struct {
  uint32_t esi;
  octmat row;
} repair_sym;

typedef kvec_t(repair_sym) repair_vec;
typedef kvec_t(unsigned) uint_vec;

#endif
