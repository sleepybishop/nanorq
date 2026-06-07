#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nanorq.h"
#include "nanorq_core.h"
#include "nanorq_ops.h"
#include "util.h"

static inline size_t div_ceil(size_t a, size_t b) { return (a + b - 1) / b; }

static inline size_t div_floor(size_t a, size_t b) { return a / b; }

struct oti_common {
  size_t F;  /* input size in bytes */
  size_t T;  /* the symbol size in octets, which MUST be a multiple of Al */
  size_t Al; /* byte alignment, 0 < Al <= 8, 4 is recommended */
};

struct oti_scheme {
  size_t Z;  /* number of source blocks */
  size_t N;  /* number of sub-blocks in each source block */
  size_t Kt; /* the total number of symbols required to represent input */
};

struct partition {
  size_t IL; /* size of long blocks */
  size_t IS; /* size of short blocks*/
  size_t JL; /* number of long blocks */
  size_t JS; /* number of short blocks */
};

struct source_block {
  size_t sbloc;
  size_t part_tot;
  struct partition part;
};

typedef struct {
  uint32_t esi;
  uint8_t *row;
} repair_sym;

typedef struct {
  repair_sym *a;
  size_t n;
  size_t m;
} repair_vec;

typedef struct {
  uint32_t *words;
  size_t size;
} compat_bitmask;

struct block_encoder {
  uint16_t K;
  bool loaded;
  bool inverted;

  // Core low level solver:
  nanorq_core core;

  // Dynamic scratch memory for core solver:
  uint8_t *prep_mem;
  size_t prep_len;
  uint8_t *work_mem;
  size_t work_len;

  uint8_t *D;      // Matrix data
  uint32_t stride; // Recommended row stride

  schedule S; // Operations schedule

  repair_vec repair_bin;
  compat_bitmask repair_mask;
};

struct nanorq {
  struct oti_common common;
  struct oti_scheme scheme;
  struct partition src_part; /* (KL, KS, ZL, ZS) = Partition[Kt, Z] */
  struct partition sub_part; /* (TL, TS, NL, NS) = Partition[T/Al, N] */
  params P;
  uint32_t max_esi;
  struct block_encoder *encoders[Z_max];

  // Cached precalculated components
  nanorq_core *precalc_core;
  uint8_t *precalc_prep_mem;
  uint8_t *precalc_work_mem;
  schedule *precalc_S;
};

static void repair_vec_push(repair_vec *v, repair_sym val) {
  if (v->n >= v->m) {
    v->m = v->m == 0 ? 8 : v->m * 2;
    v->a = (repair_sym *)realloc(v->a, v->m * sizeof(repair_sym));
  }
  v->a[v->n++] = val;
}

static void repair_vec_free(repair_vec *v) {
  if (v->a) {
    for (size_t i = 0; i < v->n; i++) {
      obl_free(v->a[i].row);
    }
    free(v->a);
  }
  v->a = NULL;
  v->n = v->m = 0;
}

static compat_bitmask compat_bitmask_new(size_t size) {
  compat_bitmask b;
  size_t num_words = (size + 31) / 32;
  b.words = (uint32_t *)calloc(num_words, sizeof(uint32_t));
  b.size = size;
  return b;
}

static void compat_bitmask_free(compat_bitmask *b) {
  free(b->words);
  b->words = NULL;
  b->size = 0;
}

static void compat_bitmask_set(compat_bitmask *b, size_t bit) {
  if (bit < b->size) {
    b->words[bit / 32] |= (1U << (bit % 32));
  }
}

static bool compat_bitmask_check(compat_bitmask *b, size_t bit) {
  if (bit < b->size) {
    return (b->words[bit / 32] & (1U << (bit % 32))) != 0;
  }
  return false;
}

static void compat_bitmask_reset(compat_bitmask *b) {
  size_t num_words = (b->size + 31) / 32;
  memset(b->words, 0, num_words * sizeof(uint32_t));
}

static size_t compat_bitmask_gaps(compat_bitmask *b, size_t limit) {
  size_t gaps = 0;
  for (size_t i = 0; i < limit; i++) {
    if (!compat_bitmask_check(b, i)) {
      if (gaps == limit)
        return limit; // prevent overflow warning
      gaps++;
    }
  }
  return gaps;
}

static struct oti_scheme gen_scheme_specific(struct oti_common *common, int K,
                                             int Z) {
  size_t Kn = K;
  struct oti_scheme ret = {0};
  ret.Kt = div_ceil(common->F, common->T);

  if (K == 0) {
    Kn = ret.Kt;
    if (Z == 0) {
      Z = 16; // if num_sbn's not specified default to at least this
      while (div_ceil(ret.Kt, Z) > K_max)
        Z++;
      while (div_ceil(ret.Kt, Z) < 10 && Z > 1)
        Z--;
    }
  }
  if (Z > 0 && K == 0) {
    Kn = div_ceil(ret.Kt, Z);
  }
  ret.Z = div_ceil(ret.Kt, Kn);
  ret.N = 1; // disable interleaving

  return ret;
}

static struct partition fill_partition(size_t I, uint16_t J) {
  struct partition p = {0, 0, 0, 0};
  if (J == 0)
    return p;
  p.IL = (size_t)(div_ceil(I, J));
  p.IS = (size_t)(div_floor(I, J));
  p.JL = (size_t)(I - p.IS * J);
  p.JS = J - p.JL;

  if (p.JL == 0)
    p.IL = 0;
  return p;
}

static struct source_block get_source_block(nanorq *rq, uint8_t sbn,
                                            uint16_t symbol_size) {
  struct source_block blk;
  blk.part = rq->sub_part;
  blk.sbloc = 0;
  blk.part_tot = rq->sub_part.IL * rq->sub_part.JL;

  if (sbn < rq->src_part.JL) {
    blk.sbloc = sbn * rq->src_part.IL * symbol_size;
  } else if (sbn - rq->src_part.JL < rq->src_part.JS) {
    blk.sbloc = (rq->src_part.IL * rq->src_part.JL) * symbol_size +
                (sbn - rq->src_part.JL) * rq->src_part.IS * symbol_size;
  }

  return blk;
}

static size_t get_symbol_offset(struct source_block *blk, size_t pos,
                                uint16_t K, uint32_t esi) {
  size_t i;
  if (pos < blk->part_tot) {
    size_t sub_blk_id = pos / blk->part.IL;
    i = blk->sbloc + sub_blk_id * K * blk->part.IL + esi * blk->part.IL +
        pos % blk->part.IL;
  } else {
    size_t pos_part2 = pos - blk->part_tot;
    size_t sub_blk_id = pos_part2 / blk->part.IS;
    i = blk->sbloc + (blk->part_tot * K) + sub_blk_id * K * blk->part.IS +
        esi * blk->part.IS + pos_part2 % blk->part.IS;
  }
  return i;
}

static size_t transfer_esi(nanorq *rq, uint8_t sbn, uint32_t esi, uint16_t K,
                           uint8_t *ptr, size_t len, struct ioctx *io,
                           int out) {
  size_t transfer = 0;
  int col = 0, symbol_size = rq->common.T / rq->common.Al;
  struct source_block blk = get_source_block(rq, sbn, symbol_size);
  for (int i = 0; i < symbol_size;) {
    size_t offset = get_symbol_offset(&blk, i, K, esi) * rq->common.Al;
    size_t sublen = (i < blk.part_tot) ? blk.part.IL : blk.part.IS;
    size_t stride = sublen * rq->common.Al;
    i += sublen;

    if (offset >= rq->common.F)
      continue;
    if (io->seek(io, offset)) {
      if ((offset + stride) >= rq->common.F)
        stride = (rq->common.F - offset);
      if (out)
        transfer += io->write(io, ptr + col, stride);
      else
        transfer += io->read(io, ptr + col, stride);
      col += stride;
    }
  }
  return transfer;
}

static struct block_encoder *get_block_encoder(nanorq *rq, uint8_t sbn) {
  if (rq->encoders[sbn])
    return rq->encoders[sbn];

  struct block_encoder *enc = (struct block_encoder *)calloc(1, sizeof(struct block_encoder));
  enc->K = nanorq_block_symbols(rq, sbn);

  enc->repair_mask = compat_bitmask_new(enc->K);

  rq->encoders[sbn] = enc;
  return enc;
}

nanorq *nanorq_encoder_new_ex(size_t len, uint16_t T, uint16_t K, uint16_t Z,
                              uint8_t Al) {
  uint8_t alignments[] = {1, 2, 4, 8};

  for (int a = sizeof(alignments) - 1; a >= 0; a--) {
    if (Al >= alignments[a]) {
      Al = alignments[a];
      break;
    }
  }
  if (Al == 0)
    Al = 1;

  if (T < Al) {
    T = Al;
  } else {
    T -= T % Al;
  }

  while (div_ceil(len, T) > Z_max * K_max) {
    if ((uint32_t)T * Al > 65535)
      return NULL;
    T *= Al;
  }
  if (len == 0 || len > NANORQ_MAX_TRANSFER) {
    return NULL;
  }

  oblas_get_impl(&nanorq_oblas);

  nanorq *rq = (nanorq *)calloc(1, sizeof(nanorq));
  rq->common.F = len;
  rq->common.T = T;
  rq->common.Al = Al;

  rq->scheme = gen_scheme_specific(&rq->common, K, Z);

  if (rq->scheme.Z == 0 || rq->scheme.N == 0 || rq->scheme.Z > Z_max ||
      div_ceil(rq->scheme.Kt, rq->scheme.Z) > K_max) {
    free(rq);
    return NULL;
  }

  rq->src_part = fill_partition(rq->scheme.Kt, rq->scheme.Z);
  rq->sub_part = fill_partition(rq->common.T / rq->common.Al, rq->scheme.N);
  rq->P = params_init(nanorq_block_symbols(rq, 0));

  rq->max_esi = (1 << 24) - 1;

  return rq;
}

nanorq *nanorq_encoder_new(size_t len, uint16_t T, uint8_t Al) {
  return nanorq_encoder_new_ex(len, T, 0, 0, Al);
}

void nanorq_free(nanorq *rq) {
  if (!rq)
    return;
  size_t num_sbn = nanorq_blocks(rq);
  for (size_t sbn = 0; sbn < num_sbn; sbn++)
    nanorq_encoder_cleanup(rq, sbn);

  if (rq->precalc_core) {
    if (rq->precalc_prep_mem)
      free(rq->precalc_prep_mem);
    if (rq->precalc_work_mem)
      free(rq->precalc_work_mem);
    if (rq->precalc_S) {
      if (rq->precalc_S->ops.a)
        free(rq->precalc_S->ops.a);
      free(rq->precalc_S);
    }
    free(rq->precalc_core);
  }

  free(rq);
}

uint64_t nanorq_oti_common(nanorq *rq) {
  uint64_t ret = 0;
  ret |= ((uint64_t)rq->common.F) << 24;
  ret |= (rq->common.T - 1) & 0xffff;
  return ret;
}

uint32_t nanorq_oti_scheme_specific(nanorq *rq) {
  uint32_t ret = 0;
  ret |= (rq->scheme.Z - 1) << 24;
  ret |= (rq->scheme.N - 1) << 8;
  ret |= rq->common.Al;
  return ret;
}

size_t nanorq_transfer_length(nanorq *rq) { return rq->common.F; }

size_t nanorq_symbol_size(nanorq *rq) { return rq->common.T; }

nanorq *nanorq_decoder_new(uint64_t common, uint32_t scheme) {
  uint64_t F = common >> 24;
  uint16_t T = (common & 0xffff) + 1;

  if (F == 0 || F > NANORQ_MAX_TRANSFER)
    return NULL;

  oblas_get_impl(&nanorq_oblas);

  nanorq *rq = (nanorq *)calloc(1, sizeof(nanorq));
  rq->common.F = F;
  rq->common.T = T;

  rq->scheme.Z = ((scheme >> 24) & 0x00ff) + 1;
  rq->scheme.N = ((scheme >> 8) & 0xffff) + 1;
  rq->common.Al = scheme & 0xff;
  rq->scheme.Kt = div_ceil(rq->common.F, rq->common.T);

  if (rq->scheme.Z == 0)
    rq->scheme.Z = Z_max;

  if (rq->scheme.N == 0) {
    rq->scheme.N = 1;
  }

  if (rq->common.Al == 0 || rq->common.T < rq->common.Al ||
      rq->common.T % rq->common.Al != 0 ||
      div_ceil(div_ceil(rq->common.F, rq->common.T), rq->scheme.Z) > K_max) {
    free(rq);
    return NULL;
  }

  rq->src_part = fill_partition(rq->scheme.Kt, rq->scheme.Z);
  rq->sub_part = fill_partition(rq->common.T / rq->common.Al, rq->scheme.N);
  rq->P = params_init(nanorq_block_symbols(rq, 0));

  rq->max_esi = (1 << 24) - 1;
  return rq;
}

size_t nanorq_block_symbols(nanorq *rq, uint8_t sbn) {
  if (sbn < rq->src_part.JL)
    return rq->src_part.IL;
  if (sbn - rq->src_part.JL < rq->src_part.JS)
    return rq->src_part.IS;
  return 0;
}

size_t nanorq_max_blocks(nanorq *rq) { return Z_max; }

size_t nanorq_blocks(nanorq *rq) {
  return (size_t)(rq->src_part.JL + rq->src_part.JS);
}

bool nanorq_precalculate(nanorq *rq) {
  if (rq->precalc_core)
    return true; // Already precalculated

  uint16_t K = nanorq_block_symbols(rq, 0);

  rq->precalc_core = (nanorq_core *)calloc(1, sizeof(nanorq_core));
  if (!nanorq_core_encoder_new(K, 0, rq->precalc_core)) {
    free(rq->precalc_core);
    rq->precalc_core = NULL;
    return false;
  }

  size_t prep_len = nanorq_core_calculate_prepare_memory(rq->precalc_core);
  rq->precalc_prep_mem = (uint8_t *)malloc(prep_len);
  if (!nanorq_core_prepare(rq->precalc_core, rq->precalc_prep_mem, prep_len)) {
    return false;
  }

  size_t work_len = nanorq_core_calculate_work_memory(rq->precalc_core);
  rq->precalc_work_mem = (uint8_t *)malloc(work_len);

  rq->precalc_S = (schedule *)calloc(1, sizeof(schedule));
  size_t sched_bytes = ops_estimate_schedule_bytes(K);
  schedule_init(rq->precalc_S, (uint8_t *)malloc(sched_bytes), sched_bytes);

  nanorq_core_set_op_callback(rq->precalc_core, rq->precalc_S, ops_push);
  if (!nanorq_core_precalculate(rq->precalc_core, rq->precalc_work_mem,
                                work_len)) {
    return false;
  }

  return true;
}

bool nanorq_generate_symbols(nanorq *rq, uint8_t sbn, struct ioctx *io) {
  struct block_encoder *enc = get_block_encoder(rq, sbn);
  if (enc == NULL)
    return false;

  if (enc->inverted)
    return true;

  if (!enc->loaded) {
    if (!enc->D) {
      if (!nanorq_core_encoder_new(enc->K, 0, &enc->core)) {
        return false;
      }
      u32 rows = nanorq_core_get_pc_rows(&enc->core);
      enc->stride = nanorq_core_recommended_stride(rq->common.T);
      enc->D = (uint8_t *)obl_alloc(rows, enc->stride, nanorq_oblas.align_size);
      nanorq_core_init_matrix(&enc->core, enc->D, enc->stride);
    }

    for (int esi = 0; esi < enc->K; esi++) {
      uint8_t *tmp_buf = (uint8_t *)malloc(rq->common.T);
      if (!tmp_buf)
        return false;
      size_t got =
          transfer_esi(rq, sbn, esi, enc->K, tmp_buf, rq->common.T, io, 0);
      if (got < rq->common.T) {
        memset(tmp_buf + got, 0, rq->common.T - got);
      }
      nanorq_core_place_symbol(&enc->core, enc->D, enc->stride, esi, tmp_buf,
                               rq->common.T);
      free(tmp_buf);
    }
    enc->loaded = true;
  }

  if (enc->K == nanorq_block_symbols(rq, 0) && rq->precalc_core) {
    ops_run(rq->precalc_core, enc->D, enc->stride, rq->precalc_S);
    enc->inverted = true;
    return true;
  }

  enc->prep_len = nanorq_core_calculate_prepare_memory(&enc->core);
  enc->prep_mem = (uint8_t *)malloc(enc->prep_len);
  if (!nanorq_core_prepare(&enc->core, enc->prep_mem, enc->prep_len)) {
    return false;
  }

  enc->work_len = nanorq_core_calculate_work_memory(&enc->core);
  enc->work_mem = (uint8_t *)malloc(enc->work_len);

  size_t sched_bytes = ops_estimate_schedule_bytes(enc->K);
  schedule_init(&enc->S, (uint8_t *)malloc(sched_bytes), sched_bytes);

  nanorq_core_set_op_callback(&enc->core, &enc->S, ops_push);
  if (!nanorq_core_precalculate(&enc->core, enc->work_mem, enc->work_len)) {
    return false;
  }

  ops_run(&enc->core, enc->D, enc->stride, &enc->S);
  enc->inverted = true;

  return true;
}

size_t nanorq_encode(nanorq *rq, void *data, uint32_t esi, uint8_t sbn,
                     struct ioctx *io) {
  struct block_encoder *enc = get_block_encoder(rq, sbn);
  if (enc == NULL)
    return 0;

  if (esi < enc->K) {
    if (enc->inverted) {
      uint8_t *tmp = (uint8_t *)malloc(enc->stride);
      if (!tmp)
        return 0;
      ops_mix(&enc->core, enc->D, enc->stride, esi, tmp);
      memcpy(data, tmp, rq->common.T);
      free(tmp);
      return rq->common.T;
    } else {
      transfer_esi(rq, sbn, esi, enc->K, (uint8_t *)data, rq->common.T, io, 0);
      return rq->common.T;
    }
  } else {
    if (esi > ((1 << 24) - 1))
      return 0;
    if (!enc->inverted) {
      if (!nanorq_generate_symbols(rq, sbn, io)) {
        return 0;
      }
    }
    uint8_t *tmp = (uint8_t *)malloc(enc->stride);
    if (!tmp)
      return 0;
    ops_mix(&enc->core, enc->D, enc->stride, esi, tmp);
    memcpy(data, tmp, rq->common.T);
    free(tmp);
    return rq->common.T;
  }
}

void nanorq_encoder_cleanup(nanorq *rq, uint8_t sbn) {
  if (!rq->encoders[sbn])
    return;
  struct block_encoder *enc = rq->encoders[sbn];
  if (enc->D) {
    obl_free(enc->D);
  }
  if (enc->prep_mem) {
    free(enc->prep_mem);
  }
  if (enc->work_mem) {
    free(enc->work_mem);
  }
  if (enc->S.ops.a) {
    free(enc->S.ops.a);
  }
  repair_vec_free(&enc->repair_bin);
  compat_bitmask_free(&enc->repair_mask);
  free(enc);
  rq->encoders[sbn] = NULL;
}

void nanorq_encoder_reset(nanorq *rq, uint8_t sbn) {
  if (!rq->encoders[sbn])
    return;
  struct block_encoder *enc = rq->encoders[sbn];
  enc->loaded = false;
  enc->inverted = false;
  if (enc->D) {
    nanorq_core_init_matrix(&enc->core, enc->D, enc->stride);
  }
  if (enc->prep_mem) {
    free(enc->prep_mem);
    enc->prep_mem = NULL;
  }
  if (enc->work_mem) {
    free(enc->work_mem);
    enc->work_mem = NULL;
  }
  if (enc->S.ops.a) {
    free(enc->S.ops.a);
    enc->S.ops.a = NULL;
  }
  repair_vec_free(&enc->repair_bin);
  compat_bitmask_reset(&enc->repair_mask);
}

bool nanorq_set_max_esi(nanorq *rq, uint32_t max_esi) {
  if (!rq || max_esi >= (1 << 24) || max_esi < rq->P.Kprime)
    return false;
  rq->max_esi = max_esi;
  return true;
}

int nanorq_decoder_add_symbol(nanorq *rq, void *data, uint32_t tag,
                              struct ioctx *io) {
  uint8_t sbn = (tag >> 24) & 0xff;
  uint32_t esi = (tag & 0x00ffffff);

  struct block_encoder *dec = get_block_encoder(rq, sbn);
  if (dec == NULL || esi >= (1 << 24) || esi > rq->max_esi)
    return NANORQ_SYM_ERR;

  if (compat_bitmask_gaps(&dec->repair_mask, dec->K) == 0) {
    return NANORQ_SYM_IGN;
  }

  if (esi < dec->K) {
    if (compat_bitmask_check(&dec->repair_mask, esi))
      return NANORQ_SYM_DUP;
  } else {
    for (size_t i = 0; i < dec->repair_bin.n; i++) {
      if (dec->repair_bin.a[i].esi == esi) {
        return NANORQ_SYM_DUP;
      }
    }
  }

  if (!dec->D) {
    if (!nanorq_core_encoder_new(dec->K, 0, &dec->core)) {
      return NANORQ_SYM_ERR;
    }
    u32 rows = nanorq_core_get_pc_rows(&dec->core);
    dec->stride = nanorq_core_recommended_stride(rq->common.T);
    dec->D = (uint8_t *)obl_alloc(rows, dec->stride, nanorq_oblas.align_size);
    nanorq_core_init_matrix(&dec->core, dec->D, dec->stride);
  }

  if (esi < dec->K) {
    nanorq_core_place_symbol(&dec->core, dec->D, dec->stride, esi, (const uint8_t *)data,
                             rq->common.T);
    transfer_esi(rq, sbn, esi, dec->K, (uint8_t *)data, rq->common.T, io, 1);
    compat_bitmask_set(&dec->repair_mask, esi);
  } else {
    repair_sym rs;
    rs.esi = esi;
    rs.row = (uint8_t *)obl_alloc(1, dec->stride, nanorq_oblas.align_size);
    memcpy(rs.row, data, rq->common.T);
    repair_vec_push(&dec->repair_bin, rs);
  }

  return NANORQ_SYM_ADDED;
}

size_t nanorq_num_missing(nanorq *rq, uint8_t sbn) {
  struct block_encoder *dec = get_block_encoder(rq, sbn);
  if (dec == NULL)
    return 0;
  return compat_bitmask_gaps(&dec->repair_mask, dec->K);
}

size_t nanorq_num_repair(nanorq *rq, uint8_t sbn) {
  struct block_encoder *dec = get_block_encoder(rq, sbn);
  if (dec == NULL)
    return 0;
  return dec->repair_bin.n;
}

uint32_t nanorq_tag(uint8_t sbn, uint32_t esi) {
  uint32_t ret = (uint32_t)(sbn) << 24;
  ret |= esi & 0x00ffffff;
  return ret;
}

bool nanorq_repair_block(nanorq *rq, struct ioctx *io, uint8_t sbn) {
  struct block_encoder *dec = get_block_encoder(rq, sbn);
  if (dec == NULL)
    return false;

  size_t num_gaps = compat_bitmask_gaps(&dec->repair_mask, dec->K);
  if (num_gaps == 0) {
    return true;
  }

  size_t num_repair = dec->repair_bin.n;
  if (num_repair < num_gaps) {
    return false;
  }

  size_t overhead = num_repair - num_gaps;

  size_t old_rows = nanorq_core_get_pc_rows(&dec->core);
  size_t new_rows = old_rows + overhead;
  if (new_rows < old_rows) {
    return false; // overflow
  }
  uint8_t *new_D = (uint8_t *)obl_alloc(new_rows, dec->stride, nanorq_oblas.align_size);
  if (!new_D) {
    return false;
  }
  if (dec->D) {
    memcpy(new_D, dec->D, old_rows * dec->stride);
    obl_free(dec->D);
  }
  memset(new_D + old_rows * dec->stride, 0, overhead * dec->stride);
  dec->D = new_D;

  if (!nanorq_core_encoder_new(dec->K, overhead, &dec->core)) {
    return false;
  }

  dec->prep_len = nanorq_core_calculate_prepare_memory(&dec->core);
  dec->prep_mem = (uint8_t *)malloc(dec->prep_len);
  if (!nanorq_core_prepare(&dec->core, dec->prep_mem, dec->prep_len)) {
    free(dec->prep_mem);
    dec->prep_mem = NULL;
    return false;
  }

  size_t rep_idx = 0;
  for (size_t gap = 0; gap < dec->K; gap++) {
    if (compat_bitmask_check(&dec->repair_mask, gap)) {
      continue;
    }
    repair_sym rs = dec->repair_bin.a[rep_idx++];
    nanorq_core_replace_symbol(&dec->core, gap, rs.esi);
    nanorq_core_place_symbol(&dec->core, dec->D, dec->stride, gap, rs.row,
                             rq->common.T);
  }

  for (size_t extra = 0; extra < overhead; extra++) {
    repair_sym rs = dec->repair_bin.a[rep_idx++];
    nanorq_core_replace_symbol(&dec->core, dec->core.P.Kprime + extra, rs.esi);
    nanorq_core_place_symbol(&dec->core, dec->D, dec->stride,
                             dec->core.P.Kprime + extra, rs.row, rq->common.T);
  }

  if (!nanorq_core_patch_matrix(&dec->core)) {
    return false;
  }

  dec->work_len = nanorq_core_calculate_work_memory(&dec->core);
  dec->work_mem = (uint8_t *)malloc(dec->work_len);

  size_t sched_bytes = ops_estimate_schedule_bytes(dec->K);
  schedule_init(&dec->S, (uint8_t *)malloc(sched_bytes), sched_bytes);

  nanorq_core_set_op_callback(&dec->core, &dec->S, ops_push);
  if (!nanorq_core_precalculate(&dec->core, dec->work_mem, dec->work_len)) {
    return false;
  }

  ops_run(&dec->core, dec->D, dec->stride, &dec->S);

  uint8_t *recovered = (uint8_t *)malloc(dec->stride);
  if (!recovered)
    return false;

  for (size_t gap = 0; gap < dec->K; gap++) {
    if (compat_bitmask_check(&dec->repair_mask, gap)) {
      continue;
    }
    ops_mix(&dec->core, dec->D, dec->stride, gap, recovered);
    transfer_esi(rq, sbn, gap, dec->K, recovered, rq->common.T, io, 1);
    compat_bitmask_set(&dec->repair_mask, gap);
  }

  free(recovered);
  return true;
}
