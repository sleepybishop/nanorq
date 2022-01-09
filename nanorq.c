#include <stdio.h>

#include "bitmask.h"
#include "nanorq.h"
#include "precode.h"
#include "tuple.h"

#define APPLYROW(D, w, a, k)                                                   \
  do {                                                                         \
    uint8_t *tmp = om_R(*D, w);                                                \
    for (int i = 0; i < k; i++)                                                \
      a[i] ^= tmp[i];                                                          \
  } while (0);

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

struct block_encoder {
  uint16_t K;
  bool loaded;
  bool inverted;
  octmat D;
  repair_vec repair_bin;
  bitmask repair_mask;
};

struct nanorq {
  struct oti_common common;
  struct oti_scheme scheme;
  struct partition src_part; /* (KL, KS, ZL, ZS) = Partition[Kt, Z] */
  struct partition sub_part; /* (TL, TS, NL, NS) = Partition[T/Al, N] */
  params P;
  uint32_t max_esi;
  schedule *S;
  struct block_encoder *encoders[Z_max];
};

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
  struct source_block ret;
  ret.part = rq->sub_part;
  ret.sbloc = 0;
  ret.part_tot = rq->sub_part.IL * rq->sub_part.JL;

  if (sbn < rq->src_part.JL) {
    ret.sbloc = sbn * rq->src_part.IL * symbol_size;
  } else if (sbn - rq->src_part.JL < rq->src_part.JS) {
    ret.sbloc = (rq->src_part.IL * rq->src_part.JL) * symbol_size +
                (sbn - rq->src_part.JL) * rq->src_part.IS * symbol_size;
  }

  return ret;
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

static struct block_encoder *get_block_encoder(nanorq *rq, uint8_t sbn) {
  if (rq->encoders[sbn])
    return rq->encoders[sbn];

  struct block_encoder *enc = calloc(1, sizeof(struct block_encoder));
  enc->K = nanorq_block_symbols(rq, sbn);

  int spare = 0;
  if (rq->max_esi) {
    enc->repair_mask = bitmask_new(rq->max_esi);
    spare = rq->max_esi - enc->K;
  }
  om_resize(&enc->D, rq->P.L + spare, rq->common.T);

  rq->encoders[sbn] = enc;
  return enc;
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

static bool load_symbol_matrix(nanorq *rq, uint8_t sbn, struct ioctx *io) {
  struct block_encoder *enc = get_block_encoder(rq, sbn);
  if (enc == NULL)
    return false;
  for (int esi = 0, row = rq->P.S + rq->P.H; esi < enc->K; esi++, row++)
    transfer_esi(rq, sbn, esi, enc->K, om_R(enc->D, row), enc->D.cols, io, 0);
  return true;
}

void decode_row(params *P, octmat *D, uint32_t row, uint8_t *ptr, size_t len) {
  memset(ptr, 0, len);
  tuple t = gen_tuple(row, P);

  // ptr might be unaligned so xor directly instead of using oblas
  APPLYROW(D, t.b, ptr, len);
  for (unsigned j = 1; j < t.d; j++) {
    t.b = (t.b + t.a) % P->W;
    APPLYROW(D, t.b, ptr, len);
  }
  while (t.b1 >= P->P)
    t.b1 = (t.b1 + t.a1) % P->P1;

  APPLYROW(D, P->W + t.b1, ptr, len);
  for (unsigned j = 1; j < t.d1; j++) {
    t.b1 = (t.b1 + t.a1) % P->P1;
    while (t.b1 >= P->P)
      t.b1 = (t.b1 + t.a1) % P->P1;
    APPLYROW(D, P->W + t.b1, ptr, len);
  }
}

bool nanorq_generate_symbols(nanorq *rq, uint8_t sbn, struct ioctx *io) {
  struct block_encoder *enc = get_block_encoder(rq, sbn);
  if (enc == NULL)
    return false;

  if (enc->inverted)
    return true;
  if (!enc->loaded)
    enc->loaded = load_symbol_matrix(rq, sbn, io);
  if (!enc->loaded)
    return false;

  schedule *S = NULL;
  if (rq->S) {
    S = rq->S;
  } else {
    spmat *A = precode_matrix_gen(&rq->P, 0);
    S = precode_matrix_invert(&rq->P, A);
  }
  if (S == NULL)
    return false;
  precode_matrix_intermediate(&rq->P, &enc->D, S);
  if (rq->S == NULL)
    sched_free(S);
  enc->inverted = true;
  return true;
}

/*
 * len: total transfer size in bytes
 * T: size of each symbol in bytes (should be aligned to Al)
 * K: number of symbols per block (set Z to zero);
 * Z: number of source blocks     (set K to zero);
 * Al: symbol alignment size
 */
nanorq *nanorq_encoder_new_ex(size_t len, uint16_t T, uint16_t K, uint16_t Z,
                              uint8_t Al) {
  nanorq *rq = NULL;

  uint8_t alignments[] = {1, 2, 4, 8};

  if (len > NANORQ_MAX_TRANSFER) {
    return NULL;
  }

  // match supported aligment
  for (int a = sizeof(alignments) - 1; a >= 0; a--) {
    if (Al >= alignments[a]) {
      Al = alignments[a];
      break;
    }
  }
  if (Al == 0)
    Al = 1;

  // enforce T is multiple of aligment
  if (T < Al) {
    T = Al;
  } else {
    T -= T % Al;
  }

  // grow symbol size to minimum needed by transfer size
  // Z_max : max number of source blocks
  // K_max : max number of symbols per source block
  while (div_ceil(len, T) > Z_max * K_max)
    T *= Al;

  rq = calloc(1, sizeof(nanorq));
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

  return rq;
}

nanorq *nanorq_encoder_new(size_t len, uint16_t T, uint8_t Al) {
  return nanorq_encoder_new_ex(len, T, 0, 0, Al);
}

void nanorq_free(nanorq *rq) {
  int num_sbn = nanorq_blocks(rq);
  if (rq) {
    if (rq->S)
      sched_free(rq->S);
    for (int sbn = 0; sbn < num_sbn; sbn++)
      nanorq_encoder_cleanup(rq, sbn);
    free(rq);
  }
}

uint64_t nanorq_oti_common(nanorq *rq) {
  uint64_t ret = 0;
  /* T is decremented by one to avoid overflow */
  ret |= ((uint64_t)rq->common.F) << 24; /* transfer length */
  ret |= (rq->common.T - 1) & 0xffff;    /* symbol size */
  return ret;
}

uint32_t nanorq_oti_scheme_specific(nanorq *rq) {
  uint32_t ret = 0;
  /* Z and N are decremented by one to avoid overflow */
  ret |= (rq->scheme.Z - 1) << 24; /* number of source blocks */
  ret |= (rq->scheme.N - 1) << 8;  /* number of sub-blocks */
  ret |= rq->common.Al;            /* symbol alignment */
  return ret;
}

uint32_t nanorq_tag(uint8_t sbn, uint32_t esi) {
  uint32_t ret = (uint32_t)(sbn) << 24;
  ret |= esi & 0x00ffffff;
  return ret;
}

size_t nanorq_transfer_length(nanorq *rq) { return rq->common.F; }

size_t nanorq_symbol_size(nanorq *rq) { return rq->common.T; }

nanorq *nanorq_decoder_new(uint64_t common, uint32_t scheme) {
  uint64_t F = common >> 24;
  /* increment T by one since it was decremented to avoid overflow */
  uint16_t T = (common & 0xffff) + 1;

  nanorq *rq = NULL;

  if (F > NANORQ_MAX_TRANSFER)
    return NULL;

  rq = calloc(1, sizeof(nanorq));

  rq->common.F = F;
  rq->common.T = T;

  /* increment Z and N by one since they were decremented to avoid overflow */
  rq->scheme.Z = ((scheme >> 24) & 0x00ff) + 1;
  rq->scheme.N = ((scheme >> 8) & 0xffff) + 1;
  rq->common.Al = scheme & 0xff;
  rq->scheme.Kt = div_ceil(rq->common.F, rq->common.T);

  if (rq->scheme.Z == 0)
    rq->scheme.Z = Z_max;

  if (rq->scheme.N == 0) {
    rq->scheme.N = 1;
  }

  if (rq->common.T < rq->common.Al || rq->common.T % rq->common.Al != 0 ||
      div_ceil(div_ceil(rq->common.F, rq->common.T), rq->scheme.Z) > K_max) {
    free(rq);
    return NULL;
  }

  rq->src_part = fill_partition(rq->scheme.Kt, rq->scheme.Z);
  rq->sub_part = fill_partition(rq->common.T / rq->common.Al, rq->scheme.N);
  rq->P = params_init(nanorq_block_symbols(rq, 0));

  rq->max_esi = 2 * rq->P.Kprime;
  return rq;
}

/* JL: num(0), JS: num(1), IL: size(0), IS: size(1) */
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
  params P = params_init(nanorq_block_symbols(rq, 0));
  spmat *A = precode_matrix_gen(&P, 0);
  schedule *S = precode_matrix_invert(&P, A);
  if (S == NULL)
    return false;
  rq->S = S;
  return true;
}

size_t nanorq_encode(nanorq *rq, void *data, uint32_t esi, uint8_t sbn,
                     struct ioctx *io) {
  size_t written = 0;
  struct block_encoder *enc = get_block_encoder(rq, sbn);
  if (enc == NULL)
    return 0;

  if (esi < enc->K) {
    if (enc->inverted) {
      decode_row(&rq->P, &enc->D, esi, data, enc->D.cols);
      written += enc->D.cols;
    } else {
      if (!enc->loaded)
        enc->loaded = load_symbol_matrix(rq, sbn, io);
      if (enc->loaded) {
        memcpy(data, om_R(enc->D, rq->P.S + rq->P.H + esi), enc->D.cols);
        written += enc->D.cols;
      }
    }
  } else {
    if (esi > ((1 << 24) - 1))
      return 0;
    // esi is for repair symbol
    if (!enc->inverted)
      enc->inverted = nanorq_generate_symbols(rq, sbn, io);
    if (enc->inverted) {
      uint32_t isi = esi + (rq->P.Kprime - enc->K);
      decode_row(&rq->P, &enc->D, isi, data, enc->D.cols);
      written += enc->D.cols;
    }
  }
  return written;
}

void nanorq_encoder_cleanup(nanorq *rq, uint8_t sbn) {
  if (!rq->encoders[sbn])
    return;
  struct block_encoder *enc = rq->encoders[sbn];
  om_destroy(&enc->D);
  if (kv_size(enc->repair_bin) > 0) {
    for (int rs = 0; rs < kv_size(enc->repair_bin); rs++)
      om_destroy(&(kv_A(enc->repair_bin, rs).row));
    kv_destroy(enc->repair_bin);
  }
  if (kv_size(enc->repair_mask) > 0)
    bitmask_free(&enc->repair_mask);
  free(enc);
  rq->encoders[sbn] = NULL;
}

void nanorq_encoder_reset(nanorq *rq, uint8_t sbn) {
  if (!rq->encoders[sbn])
    return;
  struct block_encoder *enc = rq->encoders[sbn];
  enc->loaded = false;
  enc->inverted = false;
  if (om_P(enc->D))
    memset(om_P(enc->D), 0, enc->D.rows * enc->D.cols_al);
  if (kv_size(enc->repair_bin) > 0) {
    for (int rs = 0; rs < kv_size(enc->repair_bin); rs++)
      om_destroy(&(kv_A(enc->repair_bin, rs).row));
    kv_destroy(enc->repair_bin);
    kv_init(enc->repair_bin);
  }
  if (kv_size(enc->repair_mask) > 0)
    bitmask_reset(&enc->repair_mask);
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

  if (dec == NULL || esi > rq->max_esi)
    return NANORQ_SYM_ERR;

  if (bitmask_gaps(&dec->repair_mask, dec->K) == 0) {
    return NANORQ_SYM_IGN; // no repair needed.
  }

  if (bitmask_check(&dec->repair_mask, esi))
    return NANORQ_SYM_DUP; // already got this esi

  if (esi < dec->K) {
    // write original symbol to decode mat and output stream
    memcpy(om_R(dec->D, rq->P.S + rq->P.H + esi), data, dec->D.cols);
    transfer_esi(rq, sbn, esi, dec->K, data, dec->D.cols, io, 1);
  } else {
    // save repair symbol for precode patching
    repair_sym rs = {esi, OM_INITIAL};
    om_resize(&rs.row, 1, dec->D.cols);
    memcpy(om_R(rs.row, 0), data, dec->D.cols);
    kv_push(repair_sym, dec->repair_bin, rs);
  }
  bitmask_set(&dec->repair_mask, esi);

  return NANORQ_SYM_ADDED;
}

size_t nanorq_num_missing(nanorq *rq, uint8_t sbn) {
  struct block_encoder *dec = get_block_encoder(rq, sbn);
  if (dec == NULL)
    return 0;

  return bitmask_gaps(&dec->repair_mask, dec->K);
}

size_t nanorq_num_repair(nanorq *rq, uint8_t sbn) {
  struct block_encoder *dec = get_block_encoder(rq, sbn);
  if (dec == NULL)
    return 0;

  return kv_size(dec->repair_bin);
}

static void patch_precode_matrix(params *P, spmat *A, uint16_t K, int num_gaps,
                                 bitmask *mask, repair_vec *repair_bin) {
  size_t padding = P->Kprime - K;
  int rep_idx = 0;
  for (int gap = 0; gap < P->L && num_gaps > 0; gap++) {
    if (bitmask_check(mask, gap))
      continue;
    int row = gap + P->H + P->S;
    uint32_t esi = kv_A(*repair_bin, rep_idx++).esi + padding;

    spmat_clear_row(A, row);
    params_set_idxs(esi, P, &A->idxs[row]);
    num_gaps--;
  }

  for (int row = P->L; row < A->rows; row++) {
    uint32_t esi = kv_A(*repair_bin, rep_idx++).esi + padding;
    spmat_clear_row(A, row);
    params_set_idxs(esi, P, &A->idxs[row]);
  }
}

static void fill_symbol_matrix_gaps(params *P, octmat *D, uint16_t K,
                                    bitmask *repair_mask,
                                    repair_vec *repair_bin) {
  int rep_idx = 0, skip = P->S + P->H, num_repair = kv_size(*repair_bin);
  for (int gap = 0; gap < K && rep_idx < num_repair; gap++) {
    if (bitmask_check(repair_mask, gap))
      continue;
    int row = skip + gap;
    repair_sym rs = kv_A(*repair_bin, rep_idx++);
    memcpy(om_R(*D, row), om_P(rs.row), D->cols);
  }

  for (int row = P->L; rep_idx < num_repair; row++) {
    repair_sym rs = kv_A(*repair_bin, rep_idx++);
    memcpy(om_R(*D, row), om_P(rs.row), D->cols);
  }
}

static void decode_repair_rows(params *P, octmat *D, octmat *M, uint16_t K,
                               int num_gaps, bitmask *repair_mask) {
  om_resize(M, num_gaps, D->cols);
  for (int gap = 0, row = 0; gap < K && num_gaps > 0; gap++) {
    if (bitmask_check(repair_mask, gap))
      continue;
    decode_row(P, D, gap, om_R(*M, row), M->cols);
    row++;
    num_gaps--;
  }
}

static void write_repair_rows(nanorq *rq, uint8_t sbn, uint16_t K,
                              struct ioctx *io, octmat *M,
                              bitmask *repair_mask) {
  for (int row = 0, miss_row = 0; row < K && miss_row < M->rows; row++) {
    if (bitmask_check(repair_mask, row))
      continue;
    transfer_esi(rq, sbn, row, K, om_R(*M, miss_row), M->cols, io, 1);
    bitmask_set(repair_mask, row);
    miss_row++;
  }
}

bool nanorq_repair_block(nanorq *rq, struct ioctx *io, uint8_t sbn) {
  int overhead, num_repair, num_gaps;
  struct block_encoder *dec = get_block_encoder(rq, sbn);
  if (dec == NULL)
    return 0;

  params *P = &rq->P;
  octmat *D = &dec->D;
  octmat M = OM_INITIAL;
  bitmask *repair_mask = &dec->repair_mask;
  repair_vec *repair_bin = &dec->repair_bin;

  num_repair = kv_size(*repair_bin);
  num_gaps = bitmask_gaps(repair_mask, dec->K);
  if (num_gaps == 0)
    return true;
  if (num_repair < num_gaps)
    return false;
  overhead = num_repair - num_gaps;

  if (D->rows < P->L + overhead) {
    return false;
  }

  fill_symbol_matrix_gaps(P, D, dec->K, repair_mask, repair_bin);
  spmat *A = precode_matrix_gen(P, overhead);
  patch_precode_matrix(P, A, dec->K, num_gaps, repair_mask, repair_bin);

  schedule *S = precode_matrix_invert(P, A);
  if (S == NULL) {
    om_destroy(&M);
    return false;
  }
  precode_matrix_intermediate(P, D, S);
  sched_free(S);
  decode_repair_rows(P, D, &M, dec->K, num_gaps, repair_mask);
  write_repair_rows(rq, sbn, dec->K, io, &M, repair_mask);
  om_destroy(&M);

  return (nanorq_num_missing(rq, sbn) == 0);
}
