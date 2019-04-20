#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "oblas.h"
#include "chooser.h"

struct chooser chooser_init(uint16_t tp_size) {
  struct chooser ch = {0};

  kv_init(ch.tracking);
  kv_resize(struct tracking_pair, ch.tracking, tp_size);

  kv_init(ch.r_rows);
  kv_size(ch.r_rows) = 0;

  ch.only_two_ones = false;
  return ch;
}

void chooser_clear(struct chooser *ch) {
  kv_destroy(ch->r_rows);
  kv_destroy(ch->tracking);
}

void chooser_add_tracking_pair(struct chooser *ch, bool is_hdpc,
                               size_t row_degree) {
  struct tracking_pair tp = {is_hdpc, row_degree};
  kv_push(struct tracking_pair, ch->tracking, tp);
}

uint16_t chooser_non_zero(struct chooser *ch, octmat *A, struct graph *G,
                          uint16_t i, uint16_t sub_rows, uint16_t sub_cols) {
  uint16_t non_zero = sub_cols + 1;

  for (int row = 0; row < sub_rows; row++) {
    bool next_row = false;
    int nnz = 0;
    int ones = 0;
    int ones_idx[] = {0, 0};

    onnz(om_P(*A), row + i, i, sub_cols + i, A->cols, &nnz, &ones, ones_idx);

    if (nnz > non_zero) {
      next_row = true;
      break;
    }

    if (next_row || nnz == 0)
      continue;

    // this row has same nnz as last row
    if (non_zero == nnz) {
      if (!ch->only_two_ones || ones == 2) {
        struct pair rp = {row, ones_idx[0]};
        kv_push(struct pair, ch->r_rows, rp);
      }
    } else {
      // this row has lower nnz than last row
      non_zero = nnz;
      if (ch->only_two_ones && non_zero == 1) {
        ch->only_two_ones = false;
      }
      kv_size(ch->r_rows) = 0;
      struct pair rp = {row, ones_idx[0]};
      kv_push(struct pair, ch->r_rows, rp);
    }

    if (ones == 2) {
      if (non_zero == 2) {
        if (kv_A(ch->tracking, row + i).is_hdpc == 0) {
          graph_link(G, ones_idx[0], ones_idx[1]);
        }
        if (!ch->only_two_ones) {
          ch->only_two_ones = true;
          kv_size(ch->r_rows) = 0;
          struct pair rp = {row, ones_idx[0]};
          kv_push(struct pair, ch->r_rows, rp);
        }
      }
    }
  }
  return non_zero;
}

uint16_t chooser_pick(struct chooser *ch, struct graph *G, uint16_t i,
                      uint16_t sub_rows, uint16_t non_zero) {
  uint16_t chosen = sub_rows;

  if (non_zero != 2) {
    uint16_t min_row = sub_rows;
    uint16_t min_row_hdpc = min_row;
    size_t min_degree = SIZE_MAX;
    size_t min_degree_hdpc = min_degree;
    for (size_t rp_idx = 0; rp_idx < kv_size(ch->r_rows); rp_idx++) {
      uint16_t row = kv_A(ch->r_rows, rp_idx).first;
      if (kv_A(ch->tracking, row + i).is_hdpc) {
        if (kv_A(ch->tracking, row + i).row_degree < min_degree_hdpc) {
          min_degree_hdpc = kv_A(ch->tracking, row + i).row_degree;
          min_row_hdpc = row;
        }
      } else {
        if (kv_A(ch->tracking, row + i).row_degree < min_degree) {
          min_degree = kv_A(ch->tracking, row + i).row_degree;
          min_row = row;
        }
      }
    }
    if (min_row != sub_rows) {
      chosen = min_row;
    } else {
      chosen = min_row_hdpc;
    }
  } else {
    if (ch->only_two_ones) {
      for (size_t rp_idx = 0; rp_idx < kv_size(ch->r_rows); rp_idx++) {
        if (graph_is_max(G, kv_A(ch->r_rows, rp_idx).second)) {
          chosen = kv_A(ch->r_rows, rp_idx).first;
          break;
        }
      }
    }
    if (chosen == sub_rows) {
      chosen = kv_A(ch->r_rows, 0).first;
    }
  }
  return chosen;
}
