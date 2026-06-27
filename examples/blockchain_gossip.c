#include "nanorq.h"
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BLOCK_SIZE_SYMBOLS 100
#define SYMBOL_SIZE 256
#define NUM_NODES 5

typedef struct {
  int id;
  nanorq *decoder;
  struct ioctx *dec_io;
  uint8_t *decoded_payload;
  bool has_decoded;
  uint32_t rx_count;
} p2p_node;

p2p_node *node_new(int id, uint64_t common, uint32_t specific) {
  p2p_node *n = calloc(1, sizeof(p2p_node));
  n->id = id;
  n->decoder = nanorq_decoder_new(common, specific);
  assert(n->decoder != NULL);
  nanorq_set_max_esi(n->decoder, 1000000);
  n->decoded_payload = calloc(1, BLOCK_SIZE_SYMBOLS * SYMBOL_SIZE);
  n->dec_io =
      ioctx_from_mem(n->decoded_payload, BLOCK_SIZE_SYMBOLS * SYMBOL_SIZE);
  return n;
}

void node_receive(p2p_node *n, uint8_t *coefs, uint8_t *payload, uint32_t tag) {
  if (n->has_decoded)
    return; /* ignore if synced */

  n->rx_count++;

  /* try adding the symbol */
  int res = nanorq_decoder_add_recoded_symbol(n->decoder, payload, tag, coefs,
                                              n->dec_io);
  if (res == NANORQ_SYM_ADDED || res == NANORQ_SYM_DUP) {
    size_t missing = nanorq_num_missing(n->decoder, 0);
    size_t repair = nanorq_num_repair(n->decoder, 0);
    if (repair >= missing) {
      if (nanorq_repair_block(n->decoder, n->dec_io, 0)) {
        n->has_decoded = true;
        printf("[Node %d] BLOCK SYNCED! (Received %d packets)\n", n->id,
               n->rx_count);
      }
    }
  }
}

int main() {
  uint32_t payload_len = BLOCK_SIZE_SYMBOLS * SYMBOL_SIZE;
  uint8_t *block_data = malloc(payload_len);
  for (uint32_t i = 0; i < payload_len; i++)
    block_data[i] = rand() % 256;

  struct ioctx *src_io = ioctx_from_mem(block_data, payload_len);

  /* miner who mined the block */
  nanorq *miner =
      nanorq_encoder_new_ex(payload_len, SYMBOL_SIZE, BLOCK_SIZE_SYMBOLS, 1, 1);
  printf("[Miner] Mined new block! Size: %d KB. Encoding with TSNC...\n\n",
         payload_len / 1024);

  uint64_t common = nanorq_oti_common(miner);
  uint32_t specific = nanorq_oti_scheme_specific(miner);

  /* p2p network topology */
  p2p_node *nodes[NUM_NODES];
  for (int i = 0; i < NUM_NODES; i++)
    nodes[i] = node_new(i, common, specific);

  uint8_t *coef_buf = malloc(BLOCK_SIZE_SYMBOLS);
  uint8_t *payload_buf = malloc(SYMBOL_SIZE);

  int rounds = 0;
  bool all_synced = false;

  /* simulate network ticks */
  uint32_t tag_idx = 0;
  while (!all_synced && rounds < 500) {
    rounds++;
    uint32_t current_tag = nanorq_tag(0, tag_idx++);

    /* miner streams to node 0 */
    nanorq_generate_recoded_symbol(miner, src_io, 0, rounds, coef_buf,
                                   payload_buf);
    node_receive(nodes[0], coef_buf, payload_buf, current_tag);

    /* gossip recoded packets to next node */
    for (int i = 0; i < NUM_NODES - 1; i++) {
      p2p_node *sender = nodes[i];
      p2p_node *receiver = nodes[i + 1];

      /* generate recoded packet if sender has buffered packets */
      if (sender->rx_count > 0 && !receiver->has_decoded) {
        /* simulate 10% packet loss */
        if ((rand() % 100) > 10) {
          uint32_t recoded_tag = nanorq_tag(0, tag_idx++);
          if (sender->has_decoded) {
            nanorq_generate_recoded_symbol(sender->decoder, sender->dec_io, 0,
                                           tag_idx, coef_buf, payload_buf);
          } else {
            nanorq_generate_recoded_symbol(sender->decoder, NULL, 0, 0,
                                           coef_buf, payload_buf);
          }
          node_receive(receiver, coef_buf, payload_buf, recoded_tag);
        }
      }
    }

    /* check completion */
    all_synced = true;
    for (int i = 0; i < NUM_NODES; i++) {
      if (!nodes[i]->has_decoded)
        all_synced = false;
    }
  }

  printf("\n[Network] Global consensus reached in %d ticks!\n", rounds);

  /* verification */
  int failed = 0;
  for (int i = 0; i < NUM_NODES; i++) {
    if (memcmp(block_data, nodes[i]->decoded_payload, payload_len) != 0) {
      failed++;
      printf("[Verification] Node %d FAILED! Differences:\n", i);
      for (uint32_t sym = 0; sym < BLOCK_SIZE_SYMBOLS; sym++) {
        if (memcmp(block_data + sym * SYMBOL_SIZE,
                   nodes[i]->decoded_payload + sym * SYMBOL_SIZE,
                   SYMBOL_SIZE) != 0) {
          printf("  Symbol %d mismatch:\n", sym);
          printf("    Src: ");
          for (int k = 0; k < 16; k++)
            printf("%02x ", block_data[sym * SYMBOL_SIZE + k]);
          printf("\n    Dec: ");
          for (int k = 0; k < 16; k++)
            printf("%02x ", nodes[i]->decoded_payload[sym * SYMBOL_SIZE + k]);
          printf("\n");
          break; /* only print first mismatch per node */
        }
      }
    }

    nanorq_free(nodes[i]->decoder);
    nodes[i]->dec_io->destroy(nodes[i]->dec_io);
    free(nodes[i]->decoded_payload);
    free(nodes[i]);
  }

  if (failed == 0) {
    printf("[Verification] SUCCESS: All nodes correctly decoded!\n");
  } else {
    printf("[Verification] FAILED on %d nodes!\n", failed);
  }

  free(coef_buf);
  free(payload_buf);
  free(block_data);
  src_io->destroy(src_io);
  nanorq_free(miner);

  return failed ? 1 : 0;
}
