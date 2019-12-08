#ifndef NANORQ_TUPLE_H
#define NANORQ_TUPLE_H

#include "params.h"

typedef struct {
  uint16_t d;
  uint16_t a;
  uint16_t b;
  uint16_t d1;
  uint16_t a1;
  uint16_t b1;
} tuple;

tuple gen_tuple(uint32_t X, params *P);

#endif
