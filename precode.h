#ifndef NANORQ_PRECODE_H
#define NANORQ_PRECODE_H

#include "bitmask.h"
#include "params.h"
#include "rand.h"
#include "sched.h"
#include "spmat.h"
#include "wrkmat.h"

wrkmat *precode_matrix_gen(params *P, int overhead);

bool precode_matrix_intermediate1(params *P, wrkmat *A, octmat *D);
bool precode_matrix_intermediate2(params *P, wrkmat *A, octmat *D, octmat *M,
                                  repair_vec *repair_bin, struct bitmask *mask,
                                  int num_symbols, int overhead);

void precode_matrix_fill_slot(params *P, octmat *D, uint32_t isi, uint8_t *ptr,
                              size_t len);
bool precode_matrix_decode(params *P, octmat *D, octmat *M,
                           repair_vec *repair_bin, struct bitmask *mask);

#endif
