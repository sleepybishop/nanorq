#ifndef NANORQ_GRAPH_H
#define NANORQ_GRAPH_H

#include <stdbool.h>
#include <stdint.h>

#include "util.h"

struct graph {
  pair_vec edges;
  uint16_t max_edges;
};

struct graph *graph_new(uint16_t size);
void graph_link(struct graph *g, uint16_t node_a, uint16_t node_b);
bool graph_is_max(struct graph *g, uint16_t id);
uint16_t graph_find(struct graph *g, uint16_t id);
void graph_free(struct graph *g);

#endif
