#ifndef NANORQ_PARAMS_H
#define NANORQ_PARAMS_H

#include "table2.h"
#include "util.h"
#include <stdbool.h>

typedef struct {
  uint16_t Kprime;
  uint16_t S;
  uint16_t H;
  uint16_t W;
  uint16_t L;
  uint16_t P;
  uint16_t P1;
  uint16_t U;
  uint16_t B;
  uint16_t J;
} params;

params params_init(uint16_t symbols);
uint_vec params_get_idxs(uint32_t X, params *P);

#endif
