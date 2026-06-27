#include "nanorq.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SYMBOL_SIZE 128
#define NUM_SYMBOLS 100

int main() {
  /* generate source data */
  uint32_t payload_len = NUM_SYMBOLS * SYMBOL_SIZE;
  uint8_t *source_data = malloc(payload_len);
  for (uint32_t i = 0; i < payload_len; i++) {
    source_data[i] = rand() % 256;
  }

  struct ioctx *src_io = ioctx_from_mem(source_data, payload_len);

  printf("[Encoder] Initializing TSNC Encoder with K=%d, T=%d\n", NUM_SYMBOLS,
         SYMBOL_SIZE);
  nanorq *enc =
      nanorq_encoder_new_ex(payload_len, SYMBOL_SIZE, NUM_SYMBOLS, 1, 1);

  /* Get oti common / scheme specific parameters for decoder */
  nanorq *tmp_rq =
      nanorq_encoder_new_ex(payload_len, SYMBOL_SIZE, NUM_SYMBOLS, 1, 1);
  assert(tmp_rq != NULL);
  uint64_t common = nanorq_oti_common(tmp_rq);
  uint32_t specific = nanorq_oti_scheme_specific(tmp_rq);
  nanorq_free(tmp_rq);

  printf("[Decoder] Initializing nanorq Decoder with TSNC capabilities\n");
  nanorq *dec = nanorq_decoder_new(common, specific);
  assert(dec != NULL);
  nanorq_set_max_esi(dec, 1000000);

  uint8_t *decoded_payload = calloc(1, payload_len);
  struct ioctx *dec_io = ioctx_from_mem(decoded_payload, payload_len);

  uint32_t packets_transmitted = 0;
  uint32_t packets_useful = 0;
  uint32_t packets_dropped_loss = 0;
  double loss_rate = 0.20; /* 20% packet loss */

  uint8_t *coef_buf = malloc(NUM_SYMBOLS);
  uint8_t *payload_buf = malloc(SYMBOL_SIZE);

  uint32_t esi = 0;
  nanorq_generate_symbols(enc, 0, src_io);

  bool has_decoded = false;
  while (!has_decoded) {
    memset(payload_buf, 0, SYMBOL_SIZE);
    nanorq_generate_recoded_symbol(enc, src_io, 0, esi, coef_buf, payload_buf);
    uint32_t tag = nanorq_tag(0, esi);
    esi++;
    packets_transmitted++;

    /* simulate network loss */
    double r = (double)rand() / RAND_MAX;
    if (r < loss_rate) {
      packets_dropped_loss++;
      continue;
    }

    /* try adding recoded/explicit TSNC symbol to nanorq decoder */
    int res = nanorq_decoder_add_recoded_symbol(dec, payload_buf, tag, coef_buf,
                                                dec_io);
    if (res == NANORQ_SYM_ADDED || res == NANORQ_SYM_DUP) {
      packets_useful++;
      size_t missing = nanorq_num_missing(dec, 0);
      size_t repair = nanorq_num_repair(dec, 0);
      if (repair % 10 == 0) {
        printf("  -> Decoder received %zu repair symbols, missing %zu...\n",
               repair, missing);
      }
      if (repair >= missing) {
        if (nanorq_repair_block(dec, dec_io, 0)) {
          has_decoded = true;
        }
      }
    }
  }

  printf("\n[Decoder] Decoding successful via nanorq hybrid solver!\n");
  printf("  Packets Transmitted: %u\n", packets_transmitted);
  printf("  Packets Lost (Network): %u\n", packets_dropped_loss);
  printf("  Packets Received (Useful): %u\n", packets_useful);

  if (memcmp(source_data, decoded_payload, payload_len) == 0) {
    printf("\n[Verification] SUCCESS: Decoded payloads exactly match the "
           "original source data!\n");
  } else {
    printf("\n[Verification] FAILED: Payload mismatch!\n");
    uint32_t matched_symbols = 0;
    uint32_t permuted_matches = 0;
    for (uint32_t i = 0; i < NUM_SYMBOLS; i++) {
      if (memcmp(source_data + i * SYMBOL_SIZE,
                 decoded_payload + i * SYMBOL_SIZE, SYMBOL_SIZE) == 0) {
        matched_symbols++;
      } else {
        /* check if this decoded symbol matches any source symbol */
        int found = -1;
        for (uint32_t j = 0; j < NUM_SYMBOLS; j++) {
          if (memcmp(source_data + j * SYMBOL_SIZE,
                     decoded_payload + i * SYMBOL_SIZE, SYMBOL_SIZE) == 0) {
            found = j;
            break;
          }
        }
        if (found != -1) {
          permuted_matches++;
          printf("  Decoded Symbol %d matches Source Symbol %d!\n", i, found);
        }
      }
    }
    printf("  Source vs Decoded payloads comparison (first 5 symbols, 16 bytes "
           "each):\n");
    for (uint32_t i = 0; i < 5 && i < NUM_SYMBOLS; i++) {
      printf("    Symbol %d Src: ", i);
      for (int k = 0; k < 16; k++)
        printf("%02x ", source_data[i * SYMBOL_SIZE + k]);
      printf("\n");
      printf("    Symbol %d Dec: ", i);
      for (int k = 0; k < 16; k++)
        printf("%02x ", decoded_payload[i * SYMBOL_SIZE + k]);
      printf("\n");
    }
    printf("  Total exactly matched symbols: %u / %d\n", matched_symbols,
           NUM_SYMBOLS);
    printf("  Total permuted matched symbols: %u / %d\n", permuted_matches,
           NUM_SYMBOLS);
  }

  free(coef_buf);
  free(payload_buf);
  free(source_data);
  free(decoded_payload);
  dec_io->destroy(dec_io);
  src_io->destroy(src_io);
  nanorq_free(enc);
  nanorq_free(dec);

  return 0;
}
