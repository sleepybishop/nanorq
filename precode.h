#ifndef NANORQ_PRECODE_H
#define NANORQ_PRECODE_H

#include "bitmask.h"
#include "params.h"

void precode_matrix_gen(params *P, octmat *A, uint16_t overhead);

octmat precode_matrix_intermediate1(params *P, octmat *A, octmat *D);
bool precode_matrix_intermediate2(octmat *M, octmat *A, octmat *D, params *P,
                                  repair_vec *repair_bin, struct bitmask *mask,
                                  uint16_t num_symbols, uint16_t overhead);

octmat precode_matrix_encode(params *P, octmat *C, uint32_t isi);

bool precode_matrix_decode(params *P, octmat *X, repair_vec *repair_bin,
                           struct bitmask *mask);

#endif
