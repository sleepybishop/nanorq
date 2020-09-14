#ifndef NANORQ_TUPLE_H
#define NANORQ_TUPLE_H

#include "params.h"

typedef struct {
  uint32_t d;
  uint32_t a;
  uint32_t b;
  uint32_t d1;
  uint32_t a1;
  uint32_t b1;
} tuple;

tuple gen_tuple(uint32_t X, params *P);

#endif
