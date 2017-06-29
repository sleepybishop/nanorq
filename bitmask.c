#include "bitmask.h"
#include <stdio.h>

#define IDXBITS 32

struct bitmask *bitmask_new(size_t initial_size) {
  struct bitmask *bm = calloc(1, sizeof(struct bitmask));
  size_t max_idx = (initial_size / IDXBITS) + 1;
  while (max_idx >= kv_size(bm->mask))
    kv_push(uint32_t, bm->mask, 0);

  return bm;
}

void bitmask_set(struct bitmask *bm, size_t id) {
  size_t idx = id / IDXBITS;
  while (idx >= kv_size(bm->mask))
    kv_push(uint32_t, bm->mask, 0);

  uint32_t add_mask = 1 << (id % IDXBITS);
  kv_A(bm->mask, idx) |= add_mask;
}

void bitmask_clear(struct bitmask *bm, size_t id) {
  size_t idx = id / IDXBITS;
  while (idx >= kv_size(bm->mask))
    kv_push(uint32_t, bm->mask, 0);

  uint32_t clear_mask = 1 << (id % IDXBITS);

  kv_A(bm->mask, idx) &= ~clear_mask;
}

bool bitmask_check(struct bitmask *bm, size_t id) {
  size_t idx = id / IDXBITS;
  if (idx >= kv_size(bm->mask))
    return false;

  uint32_t check_mask = 1 << (id % IDXBITS);
  return (kv_A(bm->mask, idx) & check_mask) != 0;
}

size_t bitmask_popcount(struct bitmask *bm) {
  size_t popcount = 0, idx;
  for (idx = 0; idx < kv_size(bm->mask); idx++) {
    popcount += __builtin_popcountll(kv_A(bm->mask, idx));
  }
  return popcount;
}

size_t bitmask_gaps(struct bitmask *bm, size_t until) {
  size_t gaps = 0, idx;
  size_t until_idx = (until / IDXBITS);

  until_idx = (until_idx > kv_size(bm->mask)) ? kv_size(bm->mask) : until_idx;
  for (idx = 0; idx < until_idx; idx++) {
    uint32_t target = kv_A(bm->mask, idx);
    gaps += __builtin_popcountll(~target);
  }
  if (until % IDXBITS) {
    uint32_t until_mask = (1 << (until % IDXBITS)) - 1;
    uint32_t target = kv_A(bm->mask, idx) | ~until_mask;
    gaps += __builtin_popcountll(~target);
  }

  return gaps;
}

void bitmask_free(struct bitmask *bm) {
  if (bm) {
    kv_destroy(bm->mask);
    free(bm);
  }
}

void bitmask_print(struct bitmask *bm) {
  void *memory = (void *)bm->mask.a;
  size_t extent = (bm->mask.n) * sizeof(uint32_t);
  FILE *fp = stdout;
  char c = '|', e = '\n';
  char DIGITS_BIN[] = "01";

  char *offset = (char *)(memory);
  while (extent--) {
    unsigned bits = 0;
    while (bits < 8) {
      putc(DIGITS_BIN[(*offset >> bits++) & 1], fp);
    }
    if ((extent) && (c)) {
      putc(c, fp);
    }
    offset++;
  }
  if (e) {
    putc(e, fp);
  }
  return;
}
