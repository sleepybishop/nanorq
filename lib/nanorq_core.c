#include "nanorq_core.h"
#include "precode.h"
#include "tuple.h"
#include "util.h"

#define ALIGN_VAL(x, align) (((x) + (align) - 1) & ~((align) - 1))
#define PAD_MEM(x) ALIGN_VAL(x, get_align_size())

struct oblas_impl nanorq_oblas = {.align_size = 32};

/*
 * T: size of each symbol in bytes (should be aligned to Al)
 * K: number of symbols per block
 */
bool nanorq_core_encoder_new(u32 K, u32 overhead, nanorq_core *rq) {
  if (!rq)
    return false;
  if (K < 1 || K > K_max)
    return false;
  rq->P = params_init(K);
  rq->overhead = overhead;
  oblas_get_impl(&nanorq_oblas);
  return true;
}

uint32_t nanorq_core_recommended_stride(uint32_t T) {
  uint32_t align = (uint32_t)get_align_size();
  uint32_t padded = ALIGN_VAL(T, align);
  /* avoid cache-line aliasing: if padded is an exact multiple of 64, add one
   * alignment unit */
  if (padded % 64 == 0)
    padded += align;
  return padded;
}
void nanorq_core_place_symbol(nanorq_core *rq, uint8_t *D, uint32_t stride,
                              uint32_t esi, const uint8_t *src, uint32_t T) {
  uint32_t SH = nanorq_core_get_pc_genc_offset(rq);
  uint8_t *dst = D + (SH + esi) * stride;
  memcpy(dst, src, T);
}
size_t nanorq_core_calculate_prepare_memory(nanorq_core *rq) {
  params *P = &rq->P;
  u32 mem = 0, snz = 3 * DC(P->B, P->S) + 3;

  /* c/ci, d/di & nz/cnz */
  mem += 3 * PAD_MEM(sizeof(u32) * P->L);
  mem += 3 * PAD_MEM(sizeof(u32) * (P->L + rq->overhead));
  /* NZT (one empty vec) */
  mem += 3 * PAD_MEM(sizeof(u32_vec));
  mem += 2 * PAD_MEM(sizeof(u32) * (P->L + rq->overhead));
  /* A */
  mem += (P->L + rq->overhead) * sizeof(u32_vec);
  u32 memb4a = mem;
  mem += P->S * PAD_MEM(snz * sizeof(u32));
  mem += (P->Kprime + rq->overhead) * PAD_MEM(GENC_MAX * sizeof(u32));
  /* AT */
  mem += (P->L) * sizeof(u32_vec);
  mem += (mem - memb4a) + P->L * PAD_MEM(sizeof(u32));

  return (size_t)mem;
}

static bool assign_prepare_memory(nanorq_core *rq, u8 *mem, size_t len) {
  if (!rq || !mem || len == 0)
    return false;
  params *P = &rq->P;
  pc *W = &rq->W;
  u32 snz = 3 * DC(P->B, P->S) + 3;

  W->rows = P->L + rq->overhead;
  W->cols = P->L;

  arena *a = &W->prep_mem;
  a->beg = mem;
  a->end = mem + len;
  if (!u32_vec_init(&W->c, a, W->cols, W->cols, 0))
    return false;
  if (!u32_vec_init(&W->ci, a, W->cols, W->cols, 0))
    return false;
  if (!u32_vec_init(&W->d, a, W->rows, W->rows, 0))
    return false;
  if (!u32_vec_init(&W->di, a, W->rows, W->rows, 0))
    return false;
  if (!u32_vec_init(&W->cnz, a, W->cols, W->cols, 0))
    return false;
  if (!u32_vec_init(&W->nz, a, W->rows, W->rows, 0))
    return false;

  W->NZT = alloc_array(a, u32_vec, 3);
  if (!W->NZT)
    return false;
  for (u32 i = 1; i < 3; i++)
    if (!u32_vec_init(&W->NZT[i], a, 0, W->rows, 0))
      return false;

  W->A = alloc_array(a, u32_vec, W->rows);
  if (!W->A)
    return false;
  for (u32 i = 0; i < P->S; i++)
    if (!u32_vec_init(&W->A[i], a, 0, snz, 0))
      return false;
  for (u32 i = P->S + P->H, esi = 0; i < W->rows; i++, esi++)
    if (!u32_vec_init(&W->A[i], a, 0, GENC_MAX, 0))
      return false;

  W->AT = alloc_array(a, u32_vec, W->cols);
  if (!W->AT)
    return false;

  W->cb.on_choose_arg = 0x0;
  W->cb.on_choose = precode_matrix_choose;
  W->cb.on_op_arg = 0x0;
  W->cb.on_op = precode_matrix_on_op;

  W->AT_mem_beg = NULL;
  return true;
}

bool nanorq_core_prepare(nanorq_core *rq, uint8_t *prep_mem, size_t pm_len) {
  if (!rq || !prep_mem || pm_len == 0)
    return false;
  if (!assign_prepare_memory(rq, prep_mem, pm_len))
    return false;
  precode_matrix_gen(&rq->P, &rq->W);
  return precode_matrix_prepare(&rq->P, &rq->W);
}

uint32_t nanorq_core_get_packet_mix(nanorq_core *rq, u32 esi,
                                    uint32_t *mix_idxs, uint32_t mix_max) {
  if (!rq || !mix_idxs || mix_max < GENC_MAX)
    return 0;
  params *P = &rq->P;
  u32 X = esi;
  if (esi >= P->K)
    X += (P->Kprime - P->K);
  u32_vec mix = {.m = mix_max, .n = 0, .s = 0, .a = mix_idxs};
  params_set_idxs(P, X, &mix);
  return mix.n;
}

void nanorq_core_replace_symbol(nanorq_core *rq, u32 row, u32 esi) {
  if (!rq)
    return;
  params *P = &rq->P;
  pc *W = &rq->W;
  if (!W->A)
    return;
  row += P->H + P->S;
  if (row >= W->rows)
    return;
  u32 X = esi;
  if (esi >= P->K)
    X += (P->Kprime - P->K);
  uv_clear(W->A[row]);
  params_set_idxs(P, X, &W->A[row]);
}

bool nanorq_core_patch_matrix(nanorq_core *rq) {
  if (!rq)
    return false;
  pc *W = &rq->W;
  if (!W->cnz.a || !W->nz.a || !W->NZT)
    return false;

  for (u32 i = 0; i < W->cols; i++) {
    uv_A(W->cnz, i) = 0;
  }
  for (u32 i = 0; i < W->rows; i++) {
    uv_A(W->nz, i) = 0;
  }
  uv_clear(W->NZT[1]);
  uv_clear(W->NZT[2]);

  return precode_matrix_prepare(&rq->P, W);
}

size_t nanorq_core_calculate_work_memory(nanorq_core *rq) {
  params *P = &rq->P;
  u32 mem = 0;
  u32 max_u = P->L;
  u32 rows = P->L + rq->overhead;
  /* U */
  u32 u_stride = DC(max_u, 32);
  mem += PAD_MEM(rows * sizeof(u32) * u_stride);
  /* field maps */
  mem += 2 * PAD_MEM(sizeof(u32) * rows);
  /* UL */
  mem += PAD_MEM(2 * P->H * (u_stride * 32));
  /* HDPC */
  mem += PAD_MEM(P->H * ((P->Kprime + P->S + 31) & ~31));
  /* add w->a to bump allocator */
  mem += PAD_MEM(P->Kprime + P->S);
  return (size_t)mem;
}

static bool assign_work_memory(nanorq_core *rq, u8 *mem, size_t len) {
  if (!rq || !mem || len == 0)
    return false;
  params *P = &rq->P;
  pc *W = &rq->W;
  u32 tmp = 0;

  arena *a = &W->work_mem;
  a->beg = mem;
  a->end = mem + len;

  u32 u_stride = DC(W->u, 32);
  tmp = W->rows * u_stride;
  if (!u32_vec_init(&W->U, a, tmp, tmp, u_stride))
    return false;
  if (!u32_vec_init(&W->F.rowmap, a, W->rows, W->rows, 0))
    return false;
  if (!u32_vec_init(&W->F.type, a, W->rows, W->rows, 0))
    return false;

  u32 u_aligned = u_stride * 32;
  tmp = 2 * P->H * u_aligned;
  if (!u8_vec_init(&W->UL, a, tmp, tmp, u_aligned))
    return false;

  u32 k_aligned = (P->Kprime + P->S + 31) & ~31;
  tmp = P->H * k_aligned;
  if (!u8_vec_init(&W->HDPC, a, tmp, tmp, k_aligned))
    return false;
  return true;
}

bool nanorq_core_precalculate(nanorq_core *rq, u8 *work_mem, size_t wm_len) {
  if (!rq || !work_mem || wm_len == 0)
    return false;
  if (!assign_work_memory(rq, work_mem, wm_len))
    return false;
  return precode_matrix_invert(&rq->P, &rq->W);
}

void nanorq_core_set_op_callback(nanorq_core *rq, void *arg,
                                 void (*on_op)(void *, u32, u32, u8)) {
  pc *W = &rq->W;
  W->cb.on_op_arg = arg;
  W->cb.on_op = on_op;
}

void nanorq_core_set_choose_callback(nanorq_core *rq, void *arg,
                                     u32 (*on_choose)(void *, pc *, u32, u32,
                                                      u32, u32)) {
  pc *W = &rq->W;
  W->cb.on_choose_arg = arg;
  W->cb.on_choose = on_choose;
}

void nanorq_core_init_matrix(nanorq_core *rq, uint8_t *D, uint32_t stride) {
  uint32_t rows = nanorq_core_get_pc_rows(rq);
  memset(D, 0, (size_t)rows * stride);
}
