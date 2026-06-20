#include "tsnc.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../include/util.h" /* for nanorq_oblas and gf2_8_inv */
#include "../include/rand.h" /* for pseudo-random generator */


static bool tsnc_initialized = false;

void tsnc_init(void) {
    if (!tsnc_initialized) {
        oblas_get_impl(&nanorq_oblas);
        srand(time(NULL));
        tsnc_initialized = true;
    }
}

tsnc_encoder *tsnc_encoder_new(uint32_t K, uint32_t T, const uint8_t *data) {
    tsnc_init();
    tsnc_encoder *enc = calloc(1, sizeof(tsnc_encoder));
    enc->K = K;
    enc->T = T;
    enc->source_data = malloc(K * T);
    memcpy(enc->source_data, data, K * T);
    return enc;
}

void tsnc_encoder_free(tsnc_encoder *enc) {
    free(enc->source_data);
    free(enc);
}

tsnc_recoder *tsnc_recoder_new(uint32_t K, uint32_t T, uint32_t max_packets) {
    tsnc_init();
    tsnc_recoder *rec = calloc(1, sizeof(tsnc_recoder));
    rec->K = K;
    rec->T = T;
    rec->max_packets = max_packets;
    rec->current_packets = 0;
    rec->coef_buffer = calloc(max_packets, K);
    rec->payload_buffer = calloc(max_packets, T);
    return rec;
}

void tsnc_recoder_save(tsnc_recoder *rec, const uint8_t *coefs, const uint8_t *payload) {
    if (rec->current_packets >= rec->max_packets) {
        /* ring buffer overwrite */
        uint32_t idx = rand() % rec->max_packets;
        memcpy(rec->coef_buffer + idx * rec->K, coefs, rec->K);
        memcpy(rec->payload_buffer + idx * rec->T, payload, rec->T);
    } else {
        memcpy(rec->coef_buffer + rec->current_packets * rec->K, coefs, rec->K);
        memcpy(rec->payload_buffer + rec->current_packets * rec->T, payload, rec->T);
        rec->current_packets++;
    }
}

void tsnc_recoder_reencode(tsnc_recoder *rec, uint8_t *out_coefs, uint8_t *out_payload) {
    memset(out_coefs, 0, rec->K);
    memset(out_payload, 0, rec->T);
    
    if (rec->current_packets == 0) return;

    for (uint32_t i = 0; i < rec->current_packets; i++) {
        /* Biased GF(256) Coefficient Distribution:
         * We strictly operate in a single GF(256) field end-to-end.
         * To save CPU cycles on multiplications, we intentionally skew 
         * the selection of the multiplicative identity (1) to 80%.
         * This is NOT an inner GF(2) code, but rather a sparse/biased GF(256) matrix. */
        uint8_t c;
        int r = rand() % 100;
        if (r < 80) {
            c = 1; /* 80% GF(256) identity element (fast XOR) */
        } else if (r < 90) {
            c = 0; /* 10% sparsity */
        } else {
            c = (rand() % 254) + 2; /* 10% high-entropy GF(256) coefficients */
        }

        if (c > 0) {
            nanorq_oblas.axpy(out_coefs, rec->coef_buffer + i * rec->K, c, rec->K);
            nanorq_oblas.axpy(out_payload, rec->payload_buffer + i * rec->T, c, rec->T);
        }
    }
}

void tsnc_recoder_free(tsnc_recoder *rec) {
    free(rec->coef_buffer);
    free(rec->payload_buffer);
    free(rec);
}

#include "../include/tuple.h" /* for degree_dist */

static uint32_t tsnc_deg(uint32_t v, uint32_t K) {
    for (uint32_t d = 0; d < degree_dist_size; d++) {
        if (v < degree_dist[d]) {
            return (d < K) ? d : K;
        }
    }
    return K; /* fallback */
}

void tsnc_encode_symbol_sparse(tsnc_encoder *enc, uint32_t esi, uint8_t *out_coefs, uint8_t *out_payload) {
    memset(out_payload, 0, enc->T);
    memset(out_coefs, 0, enc->K);
    
    /* get degree d */
    uint32_t v = rnd_get(esi, 0, 1048576); 
    uint32_t d = tsnc_deg(v, enc->K);
    if (d == 0) d = 1;
    
    /* select d unique indices to mix */
    uint32_t a = 1 + rnd_get(esi, 1, enc->K - 1);
    uint32_t b = rnd_get(esi, 2, enc->K);
    
    for (uint32_t j = 0; j < d; j++) {
        uint32_t idx = b % enc->K;
        
        /* Biased GF(256) Distribution: Skew towards identity element (1) */
        uint8_t coef;
        if (rnd_get(esi, 3 + j, 100) < 90) {
            coef = 1; 
        } else {
            coef = rnd_get(esi, 3 + j, 254) + 2; /* [2, 255] */
        }
        
        out_coefs[idx] = coef;
        b = (b + a) % enc->K;
    }
    
    /* mix the data with gf(256) oblas */
    for (uint32_t i = 0; i < enc->K; i++) {
        if (out_coefs[i] > 0) {
            nanorq_oblas.axpy(out_payload, enc->source_data + i * enc->T, out_coefs[i], enc->T);
        }
    }
}


tsnc_decoder *tsnc_decoder_new(uint32_t K, uint32_t T) {
    tsnc_init();
    tsnc_decoder *dec = calloc(1, sizeof(tsnc_decoder));
    dec->K = K;
    dec->T = T;
    dec->rank = 0;
    dec->matrix = calloc(K, K);
    dec->payloads = calloc(K, T);
    dec->pivot_flags = calloc(K, sizeof(bool));
    return dec;
}

bool tsnc_decoder_add_symbol(tsnc_decoder *dec, const uint8_t *coefs, const uint8_t *payload) {
    if (dec->rank == dec->K) return true;

    uint8_t *cur_coef = malloc(dec->K);
    memcpy(cur_coef, coefs, dec->K);
    uint8_t *cur_payload = malloc(dec->T);
    memcpy(cur_payload, payload, dec->T);

    /* fully reduce against existing pivots */
    for (uint32_t i = 0; i < dec->K; i++) {
        if (cur_coef[i] != 0 && dec->pivot_flags[i]) {
            uint8_t mult = cur_coef[i];
            nanorq_oblas.axpy(cur_coef, dec->matrix + i * dec->K, mult, dec->K);
            nanorq_oblas.axpy(cur_payload, dec->payloads + i * dec->T, mult, dec->T);
        }
    }

    /* find the new pivot element */
    int new_pivot = -1;
    for (uint32_t i = 0; i < dec->K; i++) {
        if (cur_coef[i] != 0) {
            new_pivot = i;
            break;
        }
    }

    if (new_pivot == -1) {
        free(cur_coef);
        free(cur_payload);
        return false; /* linearly dependent */
    }

    uint32_t i = new_pivot;

    /* scale the new pivot row */
    uint8_t inv = GF2_8_INV[cur_coef[i]];
    nanorq_oblas.scal(cur_coef, inv, dec->K);
    nanorq_oblas.scal(cur_payload, inv, dec->T);

    /* eliminate new pivot from established rows */
    for (uint32_t j = 0; j < dec->K; j++) {
        if (j == i) continue;
        if (dec->pivot_flags[j] && dec->matrix[j * dec->K + i] != 0) {
            uint8_t mult = dec->matrix[j * dec->K + i];
            nanorq_oblas.axpy(dec->matrix + j * dec->K, cur_coef, mult, dec->K);
            nanorq_oblas.axpy(dec->payloads + j * dec->T, cur_payload, mult, dec->T);
        }
    }

    /* lock row in the matrix */
    memcpy(dec->matrix + i * dec->K, cur_coef, dec->K);
    memcpy(dec->payloads + i * dec->T, cur_payload, dec->T);
    dec->pivot_flags[i] = true;
    dec->rank++;
    
    free(cur_coef);
    free(cur_payload);
    return true;
}

bool tsnc_decoder_is_decoded(tsnc_decoder *dec) {
    return dec->rank == dec->K;
}

uint8_t *tsnc_decoder_get_data(tsnc_decoder *dec) {
    if (!tsnc_decoder_is_decoded(dec)) return NULL;
    return dec->payloads; /* diagonalized payloads match source */
}

void tsnc_decoder_free(tsnc_decoder *dec) {
    free(dec->matrix);
    free(dec->payloads);
    free(dec->pivot_flags);
    free(dec);
}

