#ifndef NANORQ_CHOOSER_H
#define NANORQ_CHOOSER_H

#include <stdbool.h>

#include "graph.h"
#include "util.h"

struct tracking_pair {
  bool is_hdpc;
  size_t row_degree;
};

struct chooser {
  kvec_t(struct tracking_pair) tracking;
  pair_vec r_rows;
  bool only_two_ones;
};

struct chooser chooser_init(uint16_t tp_size);
void chooser_clear(struct chooser *ch);
void chooser_add_tracking_pair(struct chooser *ch, bool is_hdpc,
                               size_t row_degree);

uint16_t chooser_non_zero(struct chooser *ch, octmat *A, struct graph *G,
                          uint16_t i, uint16_t sub_rows, uint16_t sub_cols);
uint16_t chooser_pick(struct chooser *ch, struct graph *G, uint16_t i,
                      uint16_t sub_rows, uint16_t non_zero);

#endif
