#ifndef NANORQ_UTIL_H
#define NANORQ_UTIL_H

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

struct repair_sym {
  uint32_t esi;
  octmat row;
};

typedef kvec_t(struct repair_sym) repair_vec;
typedef kvec_t(int) int_vec;

#endif
