#ifndef NANORQ_UTIL_H
#define NANORQ_UTIL_H

#include <stdint.h>
#include <string.h>

#include <octmat.h>
#include <kvec.h>

#define div_ceil(A, B) ((A) / (B) + ((A) % (B) ? 1 : 0))
#define div_floor(A, B) ((A) / (B))

struct pair {
  uint16_t first;
  uint16_t second;
};

typedef kvec_t(struct pair) pair_vec;
typedef kvec_t(uint16_t) uint16_vec;

struct repair_sym {
  uint32_t esi;
  octmat row;
};

typedef kvec_t(struct repair_sym) repair_vec;

#endif
