#ifndef NANORQ_UTIL_H
#define NANORQ_UTIL_H

#include "obl/oblas_lite.h"
#include "uvec.h"
#include <assert.h>

#define OCT_EXP GF2_8_EXP
#define OCT_LOG GF2_8_LOG
#define OCT_INV GF2_8_INV

#define TMPSWAP(t, a, b)                                                                                                           \
    do {                                                                                                                           \
        t __tmp = a;                                                                                                               \
        a = b;                                                                                                                     \
        b = __tmp;                                                                                                                 \
    } while (0)

#endif
