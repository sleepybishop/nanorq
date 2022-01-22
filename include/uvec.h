#include <stdint.h>

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t i64;

#define DECL_VEC(type)                                                                                                             \
    typedef struct {                                                                                                               \
        u32 m;                                                                                                                     \
        u32 n;                                                                                                                     \
        u32 s;                                                                                                                     \
        type *a;                                                                                                                   \
    } type##_vec

DECL_VEC(u32);
DECL_VEC(u8);

#define DC(A, B) ((A) / (B) + ((A) % (B) ? 1 : 0))
#ifndef OBLAS_ALIGN
#define OBLAS_ALIGN sizeof(void *)
#endif
#define PAD(x) (DC(x, OBLAS_ALIGN) * OBLAS_ALIGN)

#define uv_push(v, x)                                                                                                              \
    do {                                                                                                                           \
        assert(((v).n) < (v).m);                                                                                                   \
        (v).a[(v).n++] = (x);                                                                                                      \
    } while (0)

#if !defined(NDEBUG)
#define BC(v, i) (i)
#else
#define BC(v, i) (((i) < (v).m) ? (i) : (*(int *)0 = 0))
#endif

#define uv_A(v, i) ((v).a[BC((v), (i))])
#define uv_R(v, i) ((v).a + BC((v), (i) * (v).s))
#define uv_E(v, i, j) ((v).a[BC((v), (i) * (v).s + (j))])

#define uv_pop(v) ((v).a[--(v).n])
#define uv_clear(v) ((v).n = 0)
#define uv_size(v) ((v).n)

u32 u8_vec_init(u8_vec *v, u8 *a, u32 n, u32 m, u32 s);
u32 u32_vec_init(u32_vec *v, u8 *a, u32 n, u32 m, u32 s);

#if !defined(NDEBUG)
#define bm_at(v, i, j) bm_get(v, i, j)
#else
#define bm_at(v, i, j) (((v)->a[(i) * (v)->s + ((j) / 32)] >> ((j) % 32)) & 1)
#endif

u8 bm_get(u32_vec *v, u32 i, u32 j);
void bm_set(u32_vec *v, u32 i, u32 j);
void bm_add(u32_vec *a, u32 i, u32 j);
void bm_fill(u32_vec *m, u32 i, u8 *dst);
