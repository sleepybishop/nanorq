#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include "../tsnc.h"

#define BLOCK_SIZE_SYMBOLS 100
#define SYMBOL_SIZE 256
#define NUM_NODES 5

typedef struct {
    int id;
    tsnc_decoder *decoder;
    tsnc_recoder *recoder;
    bool has_decoded;
    uint32_t rx_count;
} p2p_node;

p2p_node* node_new(int id) {
    p2p_node *n = calloc(1, sizeof(p2p_node));
    n->id = id;
    n->decoder = tsnc_decoder_new(BLOCK_SIZE_SYMBOLS, SYMBOL_SIZE);
    /* buffer 50 packets for recoding */
    n->recoder = tsnc_recoder_new(BLOCK_SIZE_SYMBOLS, SYMBOL_SIZE, 50);
    return n;
}

void node_receive(p2p_node *n, uint8_t *coefs, uint8_t *payload) {
    if (n->has_decoded) return; /* ignore if synced */

    n->rx_count++;
    /* save to recoder buffer */
    tsnc_recoder_save(n->recoder, coefs, payload);
    
    /* decode the block */
    if (tsnc_decoder_add_symbol(n->decoder, coefs, payload)) {
        if (tsnc_decoder_is_decoded(n->decoder)) {
            n->has_decoded = true;
            printf("[Node %d] BLOCK SYNCED! (Received %d packets)\n", n->id, n->rx_count);
        }
    }
}

int main() {
    uint32_t payload_len = BLOCK_SIZE_SYMBOLS * SYMBOL_SIZE;
    uint8_t *block_data = malloc(payload_len);
    for (uint32_t i = 0; i < payload_len; i++) block_data[i] = rand() % 256;

    /* miner who mined the block */
    tsnc_encoder *miner = tsnc_encoder_new(BLOCK_SIZE_SYMBOLS, SYMBOL_SIZE, block_data);
    printf("[Miner] Mined new block! Size: %d KB. Encoding with TSNC...\n\n", payload_len / 1024);

    /* p2p network topology */
    p2p_node *nodes[NUM_NODES];
    for (int i = 0; i < NUM_NODES; i++) nodes[i] = node_new(i);

    uint8_t *coef_buf = malloc(BLOCK_SIZE_SYMBOLS);
    uint8_t *payload_buf = malloc(SYMBOL_SIZE);

    int rounds = 0;
    bool all_synced = false;

    /* simulate network ticks */
    while (!all_synced && rounds < 500) {
        rounds++;
        
        /* miner streams to node 0 */
        tsnc_encode_symbol_sparse(miner, rounds, coef_buf, payload_buf);
        node_receive(nodes[0], coef_buf, payload_buf);

        /* gossip recoded packets to next node */
        for (int i = 0; i < NUM_NODES - 1; i++) {
            p2p_node *sender = nodes[i];
            p2p_node *receiver = nodes[i+1];

            /* generate recoded packet if sender has buffered packets */
            if (sender->recoder->current_packets > 0 && !receiver->has_decoded) {
                /* simulate 10% packet loss */
                if ((rand() % 100) > 10) {
                    tsnc_recoder_reencode(sender->recoder, coef_buf, payload_buf);
                    node_receive(receiver, coef_buf, payload_buf);
                }
            }
        }

        /* check completion */
        all_synced = true;
        for (int i = 0; i < NUM_NODES; i++) {
            if (!nodes[i]->has_decoded) all_synced = false;
        }
    }

    printf("\n[Network] Global consensus reached in %d ticks!\n", rounds);

    /* verification */
    for (int i = 0; i < NUM_NODES; i++) {
        uint8_t *decoded = tsnc_decoder_get_data(nodes[i]->decoder);
        assert(memcmp(block_data, decoded, payload_len) == 0);
        
        tsnc_decoder_free(nodes[i]->decoder);
        tsnc_recoder_free(nodes[i]->recoder);
        free(nodes[i]);
    }

    free(coef_buf); free(payload_buf); free(block_data);
    tsnc_encoder_free(miner);

    return 0;
}
