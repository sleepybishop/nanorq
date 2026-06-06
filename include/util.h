#ifndef NANORQ_UTIL_H
#define NANORQ_UTIL_H

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "oblas_lite.h"
#include "kvec.h"

extern struct oblas_impl nanorq_oblas;

static inline uint8_t **mat_new(size_t rows, size_t cols) {
    size_t stride = (cols + nanorq_oblas.align_size - 1) & ~(nanorq_oblas.align_size - 1);
    if (stride == 0) stride = nanorq_oblas.align_size;
    uint8_t **m = calloc(rows, sizeof(uint8_t*));
    if (rows == 0) return m;
    uint8_t *data = obl_alloc(rows, stride, nanorq_oblas.align_size);
    for (size_t i = 0; i < rows; i++) {
        m[i] = data + i * stride;
    }
    return m;
}

static inline void mat_free(uint8_t **m, size_t rows) {
    if (!m) return;
    if (rows > 0 && m[0]) {
        uint8_t *base = m[0];
        for (size_t i = 1; i < rows; i++) {
            if (m[i] < base) base = m[i];
        }
        obl_free(base);
    }
    free(m);
}

#define div_ceil(A, B) ((A) / (B) + ((A) % (B) ? 1 : 0))
#define div_floor(A, B) ((A) / (B))

#define TMPSWAP(type, a, b)                                                    \
  do {                                                                         \
    type __tmp = a;                                                            \
    a = b;                                                                     \
    b = __tmp;                                                                 \
  } while (0)

typedef struct {
  uint32_t esi;
  uint8_t *row;
} repair_sym;

typedef kvec_t(repair_sym) repair_vec;
typedef kvec_t(unsigned) uint_vec;

#endif
