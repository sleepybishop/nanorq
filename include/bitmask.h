#ifndef BITMASK_H
#define BITMASK_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "kvec.h"

typedef kvec_t(uint32_t) bitmask;

bitmask bitmask_new(size_t initial);
void bitmask_set(bitmask *bm, size_t id);
void bitmask_clear(bitmask *bm, size_t id);
bool bitmask_check(bitmask *bm, size_t id);
size_t bitmask_popcount(bitmask *bm);
size_t bitmask_gaps(bitmask *bm, size_t until);
void bitmask_reset(bitmask *bm);
void bitmask_free(bitmask *bm);
void bitmask_print(bitmask *bm);

#endif
