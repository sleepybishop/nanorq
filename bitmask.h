#ifndef BITMASK_H
#define BITMASK_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "kvec.h"

struct bitmask {
  kvec_t(uint32_t) mask;
};

struct bitmask *bitmask_new(size_t initial);
void bitmask_set(struct bitmask *bm, size_t id);
void bitmask_clear(struct bitmask *bm, size_t id);
bool bitmask_check(struct bitmask *bm, size_t id);
size_t bitmask_popcount(struct bitmask *bm);
size_t bitmask_gaps(struct bitmask *bm, size_t until);
void bitmask_free(struct bitmask *bm);
void bitmask_print(struct bitmask *bm);

#endif
