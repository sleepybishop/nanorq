#ifndef NANORQ_H
#define NANORQ_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "io.h"

static const uint64_t NANORQ_MAX_TRANSFER = 946270874880ULL; // ~881 GB

typedef struct nanorq nanorq;

// returns a new encoder configured with given parameters
nanorq *nanorq_encoder_new(uint64_t len, uint16_t T, uint8_t Al);
nanorq *nanorq_encoder_new_ex(uint64_t len, uint16_t T, uint16_t K, uint16_t Z,
                              uint8_t Al);

// returns success of generating symbols for a given sbn
bool nanorq_generate_symbols(nanorq *rq, uint8_t sbn, struct ioctx *io);

// frees up any resources used by a decoder/encoder
void nanorq_free(nanorq *rq);

// returns basic parameters to initialize a decoder
uint64_t nanorq_oti_common(nanorq *rq);

// returns extended parameters to initialize a decoder
uint32_t nanorq_oti_scheme_specific(nanorq *rq);

// returns the total size of the payload to be transferred
uint64_t nanorq_transfer_length(nanorq *rq);

// returns the size symbols in bytes (T)
uint16_t nanorq_symbol_size(nanorq *rq);

// returns number of blocks (SBN's)
uint8_t nanorq_blocks(nanorq *rq);

// returns the number of symbol rows per sbn block
uint16_t nanorq_block_symbols(nanorq *rq, uint8_t sbn);

// returns a compound symbol identifier comprised of sbn and esi
uint32_t nanorq_fid(uint8_t sbn, uint32_t esi);

// return the max number of repair symbols allowed
uint32_t nanorq_encoder_max_repair(nanorq *rq, uint8_t sbn);

// return the number of bytes written for a given sbn and esi encode request
uint64_t nanorq_encode(nanorq *rq, void *data, uint32_t esi, uint8_t sbn,
                       struct ioctx *io);

// cleanup encoder resouces of a given block
void nanorq_encode_cleanup(nanorq *rq, uint8_t sbn);

// returns a new decoder initialized with given parameters
nanorq *nanorq_decoder_new(uint64_t common, uint32_t specific);

// returns the success of adding a symbol to the decoder
bool nanorq_decoder_add_symbol(nanorq *rq, void *data, uint32_t fid);

// returns number of symbol gaps in decoder for given block
uint32_t nanorq_num_missing(nanorq *rq, uint8_t sbn);

// returns number of repair symbols in decoder for given block
uint32_t nanorq_num_repair(nanorq *rq, uint8_t sbn);

// returns the number of bytes written from decoding a given sbn
uint64_t nanorq_decode_block(nanorq *rq, struct ioctx *io, uint8_t sbn);

// cleanup decoder resouces of a given block
void nanorq_decode_cleanup(nanorq *rq, uint8_t sbn);

#endif
