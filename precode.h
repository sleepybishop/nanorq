#ifndef NANORQ_PRECODE_H
#define NANORQ_PRECODE_H

#include "params.h"
#include "rand.h"
#include "sched.h"
#include "spmat.h"
#include "wrkmat.h"

spmat *precode_matrix_gen(params *P, int overhead);
bool precode_matrix_intermediate(params *P, spmat *A, octmat *D);
schedule *precode_matrix_invert(params *P, spmat *A);

#endif
