#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../tsnc.h"

#define SYMBOL_SIZE 128
#define NUM_SYMBOLS 50

int main() {
    uint32_t payload_len = NUM_SYMBOLS * SYMBOL_SIZE;
    uint8_t *source_data = malloc(payload_len);
    for (uint32_t i = 0; i < payload_len; i++) source_data[i] = rand() % 256;

    printf("[Source Node] Initializing Encoder (K=%d)\n", NUM_SYMBOLS);
    tsnc_encoder *enc = tsnc_encoder_new(NUM_SYMBOLS, SYMBOL_SIZE, source_data);

    printf("[Relay 1] Initializing Recoder (buffer capacity=30)\n");
    tsnc_recoder *relay1 = tsnc_recoder_new(NUM_SYMBOLS, SYMBOL_SIZE, 30);
    
    printf("[Relay 2] Initializing Recoder (buffer capacity=30)\n");
    tsnc_recoder *relay2 = tsnc_recoder_new(NUM_SYMBOLS, SYMBOL_SIZE, 30);

    printf("[Dest Node] Initializing Decoder\n");
    tsnc_decoder *dec = tsnc_decoder_new(NUM_SYMBOLS, SYMBOL_SIZE);

    uint8_t *coef_buf = malloc(NUM_SYMBOLS);
    uint8_t *payload_buf = malloc(SYMBOL_SIZE);

    /* phase 1: source to relays */
    printf("\n--- Phase 1: Source transmitting to Relays ---\n");
    uint32_t esi = 0;
    for (int i=0; i<60; i++) {
        tsnc_encode_symbol_sparse(enc, esi++, coef_buf, payload_buf);
        /* relay 1 gets first half */
        if (i < 30) {
            tsnc_recoder_save(relay1, coef_buf, payload_buf);
        } else {
            /* relay 2 gets second half */
            tsnc_recoder_save(relay2, coef_buf, payload_buf);
        }
    }
    printf("  Relay 1 buffered %d packets\n", relay1->current_packets);
    printf("  Relay 2 buffered %d packets\n", relay2->current_packets);

    /* phase 2: relays re-encode and transmit to destination */
    printf("\n--- Phase 2: Relays Multipath Re-encoding to Destination ---\n");
    uint32_t relay1_tx = 0, relay2_tx = 0;
    
    while (!tsnc_decoder_is_decoded(dec)) {
        /* relay 1 transmits with 30% loss */
        tsnc_recoder_reencode(relay1, coef_buf, payload_buf);
        relay1_tx++;
        if ((double)rand() / RAND_MAX > 0.30) {
            tsnc_decoder_add_symbol(dec, coef_buf, payload_buf);
        }

        if (tsnc_decoder_is_decoded(dec)) break;

        /* relay 2 transmits with 50% loss */
        tsnc_recoder_reencode(relay2, coef_buf, payload_buf);
        relay2_tx++;
        if ((double)rand() / RAND_MAX > 0.50) {
            tsnc_decoder_add_symbol(dec, coef_buf, payload_buf);
        }
    }

    printf("\n[Dest Node] Decoding successful!\n");
    printf("  Relay 1 Transmitted: %u\n", relay1_tx);
    printf("  Relay 2 Transmitted: %u\n", relay2_tx);

    uint8_t *decoded_data = tsnc_decoder_get_data(dec);
    if (memcmp(source_data, decoded_data, payload_len) == 0) {
        printf("\n[Verification] SUCCESS: Destination decoded payload exactly matches Source!\n");
    } else {
        printf("\n[Verification] FAILED: Payload mismatch!\n");
    }

    free(coef_buf); free(payload_buf); free(source_data);
    tsnc_encoder_free(enc); tsnc_decoder_free(dec);
    tsnc_recoder_free(relay1); tsnc_recoder_free(relay2);

    return 0;
}
