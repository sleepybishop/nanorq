#include "nanorq.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SYMBOL_SIZE 128
#define NUM_SYMBOLS 50

int main() {
  uint32_t payload_len = NUM_SYMBOLS * SYMBOL_SIZE;
  uint8_t *source_data = malloc(payload_len);
  for (uint32_t i = 0; i < payload_len; i++)
    source_data[i] = rand() % 256;

  struct ioctx *src_io = ioctx_from_mem(source_data, payload_len);

  printf("[Source Node] Initializing Encoder (K=%d)\n", NUM_SYMBOLS);
  nanorq *enc =
      nanorq_encoder_new_ex(payload_len, SYMBOL_SIZE, NUM_SYMBOLS, 1, 1);
  assert(enc != NULL);
  uint64_t common = nanorq_oti_common(enc);
  uint32_t specific = nanorq_oti_scheme_specific(enc);

  printf("[Relay 1] Initializing Recoder\n");
  nanorq *relay1 = nanorq_decoder_new(common, specific);
  nanorq_set_max_esi(relay1, 1000000);

  printf("[Relay 2] Initializing Recoder\n");
  nanorq *relay2 = nanorq_decoder_new(common, specific);
  nanorq_set_max_esi(relay2, 1000000);

  printf("[Dest Node] Initializing nanorq Decoder with TSNC capabilities\n");
  nanorq *dec = nanorq_decoder_new(common, specific);
  assert(dec != NULL);
  nanorq_set_max_esi(dec, 1000000);

  uint8_t *decoded_payload = calloc(1, payload_len);
  struct ioctx *dec_io = ioctx_from_mem(decoded_payload, payload_len);

  uint8_t *coef_buf = malloc(NUM_SYMBOLS);
  uint8_t *payload_buf = malloc(SYMBOL_SIZE);

  /* phase 1: source to relays */
  printf("\n--- Phase 1: Source transmitting to Relays ---\n");
  uint32_t esi = 0;
  uint32_t r1_esi = 0, r2_esi = 0;
  for (int i = 0; i < 60; i++) {
    nanorq_generate_recoded_symbol(enc, src_io, 0, esi++, coef_buf,
                                   payload_buf);
    /* relay 1 gets first half */
    if (i < 30) {
      uint32_t tag = nanorq_tag(0, r1_esi++);
      nanorq_decoder_add_recoded_symbol(relay1, payload_buf, tag, coef_buf,
                                        NULL);
    } else {
      /* relay 2 gets second half */
      uint32_t tag = nanorq_tag(0, r2_esi++);
      nanorq_decoder_add_recoded_symbol(relay2, payload_buf, tag, coef_buf,
                                        NULL);
    }
  }
  printf("  Relay 1 buffered packets\n");
  printf("  Relay 2 buffered packets\n");

  /* phase 2: relays re-encode and transmit to destination */
  printf("\n--- Phase 2: Relays Multipath Re-encoding to Destination ---\n");
  uint32_t relay1_tx = 0, relay2_tx = 0;

  bool has_decoded = false;
  uint32_t dest_tag_idx = 0;
  while (!has_decoded) {
    /* relay 1 transmits with 30% loss */
    nanorq_generate_recoded_symbol(relay1, NULL, 0, relay1_tx, coef_buf,
                                   payload_buf);
    relay1_tx++;
    if ((double)rand() / RAND_MAX > 0.30) {
      uint32_t tag = nanorq_tag(0, dest_tag_idx++);
      nanorq_decoder_add_recoded_symbol(dec, payload_buf, tag, coef_buf,
                                        dec_io);
    }

    size_t missing = nanorq_num_missing(dec, 0);
    size_t repair = nanorq_num_repair(dec, 0);
    if (repair >= missing && missing > 0) {
      if (nanorq_repair_block(dec, dec_io, 0)) {
        has_decoded = true;
        break;
      }
    }

    /* relay 2 transmits with 50% loss */
    nanorq_generate_recoded_symbol(relay2, NULL, 0, relay2_tx, coef_buf,
                                   payload_buf);
    relay2_tx++;
    if ((double)rand() / RAND_MAX > 0.50) {
      uint32_t tag = nanorq_tag(0, dest_tag_idx++);
      nanorq_decoder_add_recoded_symbol(dec, payload_buf, tag, coef_buf,
                                        dec_io);
    }

    missing = nanorq_num_missing(dec, 0);
    repair = nanorq_num_repair(dec, 0);
    if (repair >= missing && missing > 0) {
      if (nanorq_repair_block(dec, dec_io, 0)) {
        has_decoded = true;
        break;
      }
    }
  }

  printf("\n[Dest Node] Decoding successful!\n");
  printf("  Relay 1 Transmitted: %u\n", relay1_tx);
  printf("  Relay 2 Transmitted: %u\n", relay2_tx);

  if (memcmp(source_data, decoded_payload, payload_len) == 0) {
    printf("\n[Verification] SUCCESS: Destination decoded payload exactly "
           "matches Source!\n");
  } else {
    printf("\n[Verification] FAILED: Payload mismatch!\n");
  }

  free(coef_buf);
  free(payload_buf);
  free(source_data);
  free(decoded_payload);
  dec_io->destroy(dec_io);
  src_io->destroy(src_io);
  nanorq_free(enc);
  nanorq_free(dec);
  nanorq_free(relay1);
  nanorq_free(relay2);

  return 0;
}
