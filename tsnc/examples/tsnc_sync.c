#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../tsnc.h"

#define SYMBOL_SIZE 128
#define NUM_SYMBOLS 100

int main() {
    /* generate source data */
    uint32_t payload_len = NUM_SYMBOLS * SYMBOL_SIZE;
    uint8_t *source_data = malloc(payload_len);
    for (uint32_t i = 0; i < payload_len; i++) {
        source_data[i] = rand() % 256;
    }

    printf("[Encoder] Initializing TSNC Encoder with K=%d, T=%d\n", NUM_SYMBOLS, SYMBOL_SIZE);
    tsnc_encoder *enc = tsnc_encoder_new(NUM_SYMBOLS, SYMBOL_SIZE, source_data);

    printf("[Decoder] Initializing TSNC Decoder\n");
    tsnc_decoder *dec = tsnc_decoder_new(NUM_SYMBOLS, SYMBOL_SIZE);

    uint32_t packets_transmitted = 0;
    uint32_t packets_useful = 0;
    uint32_t packets_dropped_loss = 0;
    double loss_rate = 0.20; /* 20% packet loss */

    uint8_t *coef_buf = malloc(NUM_SYMBOLS);
    uint8_t *payload_buf = malloc(SYMBOL_SIZE);

    uint32_t esi = 0;
    while (!tsnc_decoder_is_decoded(dec)) {
        tsnc_encode_symbol_sparse(enc, esi++, coef_buf, payload_buf);
        packets_transmitted++;

        /* simulate network loss */
        double r = (double)rand() / RAND_MAX;
        if (r < loss_rate) {
            packets_dropped_loss++;
            continue;
        }

        /* try decoding */
        if (tsnc_decoder_add_symbol(dec, coef_buf, payload_buf)) {
            packets_useful++;
            if (dec->rank % 10 == 0) {
                printf("  -> Decoder received %d/%d independent symbols...\n", dec->rank, NUM_SYMBOLS);
            }
        }
    }

    printf("\n[Decoder] Decoding successful!\n");
    printf("  Packets Transmitted: %u\n", packets_transmitted);
    printf("  Packets Lost (Network): %u\n", packets_dropped_loss);
    printf("  Packets Linearly Dependent (Wasted): %u\n", packets_transmitted - packets_dropped_loss - packets_useful);
    printf("  Packets Independent (Useful): %u\n", packets_useful);

    uint8_t *decoded_data = tsnc_decoder_get_data(dec);
    assert(decoded_data != NULL);

    if (memcmp(source_data, decoded_data, payload_len) == 0) {
        printf("\n[Verification] SUCCESS: Decoded payloads exactly match the original source data!\n");
    } else {
        printf("\n[Verification] FAILED: Payload mismatch!\n");
    }

    free(coef_buf);
    free(payload_buf);
    free(source_data);
    tsnc_encoder_free(enc);
    tsnc_decoder_free(dec);

    return 0;
}
