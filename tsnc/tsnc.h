#ifndef TSNC_H
#define TSNC_H

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t K; /* source symbols */
    uint32_t T; /* symbol size */
    uint8_t *source_data; 
} tsnc_encoder;

typedef struct {
    uint32_t K;
    uint32_t T;
    uint32_t rank;
    uint8_t *matrix; /* coefficient matrix */
    uint8_t *payloads; /* payload matrix */
    bool *pivot_flags; /* pivot flags */
} tsnc_decoder;

typedef struct {
    uint32_t K;
    uint32_t T;
    uint32_t max_packets;
    uint32_t current_packets;
    uint8_t *coef_buffer;
    uint8_t *payload_buffer;
} tsnc_recoder;

/* initialize oblas */
void tsnc_init(void);

tsnc_encoder *tsnc_encoder_new(uint32_t K, uint32_t T, const uint8_t *data);
void tsnc_encode_symbol_sparse(tsnc_encoder *enc, uint32_t esi, uint8_t *out_coefs, uint8_t *out_payload);
void tsnc_encoder_free(tsnc_encoder *enc);

tsnc_recoder *tsnc_recoder_new(uint32_t K, uint32_t T, uint32_t max_packets);
void tsnc_recoder_save(tsnc_recoder *rec, const uint8_t *coefs, const uint8_t *payload);
void tsnc_recoder_reencode(tsnc_recoder *rec, uint8_t *out_coefs, uint8_t *out_payload);
void tsnc_recoder_free(tsnc_recoder *rec);

tsnc_decoder *tsnc_decoder_new(uint32_t K, uint32_t T);
bool tsnc_decoder_add_symbol(tsnc_decoder *dec, const uint8_t *coefs, const uint8_t *payload);
bool tsnc_decoder_is_decoded(tsnc_decoder *dec);
uint8_t *tsnc_decoder_get_data(tsnc_decoder *dec);
void tsnc_decoder_free(tsnc_decoder *dec);


#ifdef __cplusplus
}
#endif

#endif // TSNC_H
