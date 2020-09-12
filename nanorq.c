#include <stdio.h>

#include "nanorq.h"
#include "precode.h"

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
  size_t al;
};

struct block_encoder {
  size_t sbn;
  size_t num_symbols;
  size_t symbol_size;
  bool loaded;
  bool inverted;
  params P;
  octmat D;
  repair_vec repair_bin;
  bitmask repair_mask;
};

struct nanorq {
  struct oti_common common;
  struct oti_scheme scheme;
  struct partition src_part; /* (KL, KS, ZL, ZS) = Partition[Kt, Z] */
  struct partition sub_part; /* (TL, TS, NL, NS) = Partition[T/Al, N] */
  struct block_encoder *encoders[Z_max];
};

static struct oti_scheme gen_scheme_specific(struct oti_common *common, int K,
                                             int Z) {
  size_t Kn = K;
  struct oti_scheme ret = {0};
  ret.Kt = div_ceil(common->F, common->T);

  if (K == 0) {
    Kn = ret.Kt;
    /*
     * user more sbns by unless otherwise specified
     * performance is better for small K
     *
     */
    if (Z == 0) {
      Z = 16;
      while (div_ceil(ret.Kt, Z) > K_max)
        Z++;
    }
  }
  if (Z > 0 && K == 0) {
    Kn = div_ceil(ret.Kt, Z);
  }
  ret.Z = div_ceil(ret.Kt, Kn);
  // disable interleaving
  ret.N = 1;

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
  ret.al = rq->common.Al;
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
                                uint16_t K, uint32_t symbol_id) {
  size_t i;

  if (pos < blk->part_tot) {
    size_t sub_blk_id = pos / blk->part.IL;
    i = blk->sbloc + sub_blk_id * K * blk->part.IL + symbol_id * blk->part.IL +
        pos % blk->part.IL;
  } else {
    size_t pos_part2 = pos - blk->part_tot;
    size_t sub_blk_id = pos_part2 / blk->part.IS;
    i = blk->sbloc + (blk->part_tot * K) + sub_blk_id * K * blk->part.IS +
        symbol_id * blk->part.IS + pos_part2 % blk->part.IS;
  }

  return i * blk->al;
}

static struct block_encoder *nanorq_block_encoder(nanorq *rq, uint8_t sbn,
                                                  int rp) {
  size_t num_symbols = nanorq_block_symbols(rq, sbn);
  size_t symbol_size = rq->common.T / rq->common.Al;

  if (rq->encoders[sbn])
    return rq->encoders[sbn];

  if (num_symbols == 0 || symbol_size == 0)
    return NULL;

  struct block_encoder *enc = calloc(1, sizeof(struct block_encoder));
  enc->sbn = sbn;
  enc->num_symbols = num_symbols;
  enc->symbol_size = symbol_size;
  enc->P = params_init(num_symbols);

  int matrows = enc->P.S + enc->P.H + enc->P.Kprime;
  if (rp) {
    enc->repair_mask = bitmask_new(num_symbols);
    matrows += (matrows / 5); // estimate 20 pct overhead
  }
  om_resize(&enc->D, matrows, symbol_size * rq->common.Al);

  rq->encoders[sbn] = enc;
  return enc;
}

static bool load_sbn(nanorq *rq, uint8_t sbn, struct ioctx *io) {
  struct block_encoder *enc = nanorq_block_encoder(rq, sbn, 0);
  if (enc == NULL)
    return false;

  params *P = &enc->P;
  octmat *D = &enc->D;
  struct source_block blk = get_source_block(rq, sbn, enc->symbol_size);
  int row = P->S + P->H;
  for (; row < P->S + P->H + enc->num_symbols; row++) {
    uint32_t esi = row - (P->S + P->H);
    int col = 0;
    for (int i = 0; i < enc->symbol_size;) {
      size_t offset = get_symbol_offset(&blk, i, enc->num_symbols, esi);
      size_t sublen = (i < blk.part_tot) ? blk.part.IL : blk.part.IS;
      size_t stride = sublen * rq->common.Al;
      size_t got = 0;
      i += sublen;
      if (io->seek(io, offset))
        got = io->read(io, om_R(*D, row) + col, stride);
      col += stride;
      for (int byte = got; byte < stride; byte++)
        om_A(*D, row, col++) = 0; // padding
    }
  }
  return true;
}

bool nanorq_generate_symbols(nanorq *rq, uint8_t sbn, struct ioctx *io) {
  struct block_encoder *enc = nanorq_block_encoder(rq, sbn, 0);
  if (enc == NULL)
    return false;

  if (enc->inverted)
    return true;
  if (!enc->loaded)
    enc->loaded = load_sbn(rq, sbn, io);
  if (!enc->loaded)
    return false;

  params *P = &enc->P;
  octmat *D = &enc->D;
  spmat *A = precode_matrix_gen(P, 0);
  if (!precode_matrix_intermediate1(P, A, D))
    return false;
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
nanorq *nanorq_encoder_new_ex(uint64_t len, uint16_t T, uint16_t K, uint16_t Z,
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

#ifdef NANORQ_DEBUG
  fprintf(stderr, "T: %06d AL: %d \n", T, Al);
  fprintf(stderr, "Z: %06d N : %06d Kt: %06d\n", (int)rq->scheme.Z,
          (int)rq->scheme.N, (int)rq->scheme.Kt);

  fprintf(stderr, "P1 %dx%d P2 %dx%d\n", rq->src_part.JL, rq->src_part.IL,
          rq->src_part.JS, rq->src_part.IS);
#endif

  return rq;
}

nanorq *nanorq_encoder_new(uint64_t len, uint16_t T, uint8_t Al) {
  return nanorq_encoder_new_ex(len, T, 0, 0, Al);
}

void nanorq_free(nanorq *rq) {
  int num_sbn = nanorq_blocks(rq);
  if (rq) {
    for (int sbn = 0; sbn < num_sbn; sbn++) {
      nanorq_encoder_cleanup(rq, sbn);
    }
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

uint64_t nanorq_transfer_length(nanorq *rq) { return rq->common.F; }

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

#ifdef NANORQ_DEBUG
  fprintf(stderr, "T: %06d AL: %d \n", rq->common.T, rq->common.Al);
  fprintf(stderr, "Z: %06d N : %06d Kt: %06d\n", (int)rq->scheme.Z,
          (int)rq->scheme.N, (int)rq->scheme.Kt);

  fprintf(stderr, "P1 %dx%d P2 %dx%d\n", rq->src_part.JL, rq->src_part.IL,
          rq->src_part.JS, rq->src_part.IS);
#endif

  return rq;
}

/*
+   num(0) JL
+   num(1) JS
+   size(0) IL
+   size(1) IS
*/

size_t nanorq_block_symbols(nanorq *rq, uint8_t sbn) {
  if (sbn < rq->src_part.JL)
    return rq->src_part.IL;
  if (sbn - rq->src_part.JL < rq->src_part.JS)
    return rq->src_part.IS;
  return 0;
}

size_t nanorq_encoder_max_repair(nanorq *rq, uint8_t sbn) {
  return (uint32_t)((1 << 20) - nanorq_block_symbols(rq, sbn));
}

size_t nanorq_blocks(nanorq *rq) {
  return (int)(rq->src_part.JL + rq->src_part.JS);
}

uint64_t nanorq_encode(nanorq *rq, void *data, uint32_t esi, uint8_t sbn,
                       struct ioctx *io) {
  uint64_t written = 0;
  struct block_encoder *enc = nanorq_block_encoder(rq, sbn, 0);
  if (enc == NULL)
    return 0;

  uint8_t tmp[enc->D.cols_al];
  params *P = &enc->P;
  memset(tmp, 0, enc->D.cols_al);
  if (esi < enc->num_symbols) {
    if (enc->inverted) {
      precode_matrix_fill_slot(P, &enc->D, esi, tmp, enc->D.cols);
      memcpy(data, tmp, enc->D.cols);
      written += enc->D.cols;
    } else {
      if (!enc->loaded)
        enc->loaded = load_sbn(rq, sbn, io);
      if (enc->loaded) {
        memcpy(data, om_R(enc->D, enc->P.S + enc->P.H + esi), enc->D.cols);
        written += enc->D.cols;
      }
    }
  } else {
    // esi is for repair symbol
    if (!enc->inverted)
      enc->inverted = nanorq_generate_symbols(rq, sbn, io);
    if (enc->inverted) {
      uint32_t isi = esi + (P->Kprime - enc->num_symbols);
      precode_matrix_fill_slot(P, &enc->D, isi, tmp, enc->D.cols);
      memcpy(data, tmp, enc->D.cols);
      written += enc->D.cols;
    }
  }
  return written;
}

void nanorq_encoder_cleanup(nanorq *rq, uint8_t sbn) {
  if (rq->encoders[sbn]) {
    struct block_encoder *enc = rq->encoders[sbn];
    om_destroy(&enc->D);
    if (kv_size(enc->repair_bin) > 0) {
      for (int rs = 0; rs < kv_size(enc->repair_bin); rs++)
        om_destroy(&(kv_A(enc->repair_bin, rs).row));
      kv_destroy(enc->repair_bin);
    }
    if (kv_size(enc->repair_mask) > 0) {
      bitmask_free(&enc->repair_mask);
    }
    free(enc);
    rq->encoders[sbn] = NULL;
  }
}

uint64_t nanorq_decode_write_esi(nanorq *rq, struct ioctx *io, uint8_t sbn,
                                 uint32_t esi, uint8_t *ptr, size_t dlen) {
  struct block_encoder *dec = nanorq_block_encoder(rq, sbn, 1);
  if (dec == NULL)
    return 0;

  struct source_block blk = get_source_block(rq, sbn, dec->symbol_size);

  uint64_t written = 0;
  int col = 0;
  for (int i = 0; i < dec->symbol_size;) {
    size_t offset = get_symbol_offset(&blk, i, dec->num_symbols, esi);
    size_t sublen = (i < blk.part_tot) ? blk.part.IL : blk.part.IS;
    size_t stride = sublen * rq->common.Al;
    i += sublen;

    if (offset >= rq->common.F)
      continue;
    if (io->seek(io, offset)) {
      size_t len = stride;
      if ((offset + stride) >= rq->common.F) {
        len = (rq->common.F - offset);
      }
      written += io->write(io, ptr + col, len);
      col += stride;
    }
  }
  return written;
}

bool nanorq_decoder_add_symbol(nanorq *rq, void *data, uint32_t tag,
                               struct ioctx *io) {
  uint8_t sbn = (tag >> 24) & 0xff;
  uint32_t esi = (tag & 0x00ffffff);

  struct block_encoder *dec = nanorq_block_encoder(rq, sbn, 1);

  if (dec == NULL)
    return false;

  size_t cols = dec->D.cols;

  if (esi >= (1 << 20))
    return false;

  if (bitmask_gaps(&dec->repair_mask, dec->num_symbols) == 0) {
    return true; // no gaps! no repair needed.
  }

  if (bitmask_check(&dec->repair_mask, esi))
    return true; // already got this esi

  if (esi < dec->num_symbols) {
    // write original symbol to decode mat and output stream
    memcpy(om_R(dec->D, dec->P.S + dec->P.H + esi), data, cols);
    nanorq_decode_write_esi(rq, io, sbn, esi, data, cols);
  } else {
    // save repair symbol for precode patching
    repair_sym rs = {esi, OM_INITIAL};
    om_resize(&rs.row, 1, cols);
    memcpy(om_R(rs.row, 0), data, cols);
    kv_push(repair_sym, dec->repair_bin, rs);
  }
  bitmask_set(&dec->repair_mask, esi);

  return true;
}

size_t nanorq_num_missing(nanorq *rq, uint8_t sbn) {
  size_t num_symbols = nanorq_block_symbols(rq, sbn);
  struct block_encoder *dec = nanorq_block_encoder(rq, sbn, 1);
  if (dec == NULL)
    return 0;

  return bitmask_gaps(&dec->repair_mask, num_symbols);
}

size_t nanorq_num_repair(nanorq *rq, uint8_t sbn) {
  struct block_encoder *dec = nanorq_block_encoder(rq, sbn, 1);
  if (dec == NULL)
    return 0;

  return kv_size(dec->repair_bin);
}

bool nanorq_repair_block(nanorq *rq, struct ioctx *io, uint8_t sbn) {
  struct block_encoder *dec = nanorq_block_encoder(rq, sbn, 1);
  if (dec == NULL)
    return 0;

  params *P = &dec->P;
  octmat M = OM_INITIAL;
  bool success = precode_matrix_decode(P, &dec->D, &M, &dec->repair_bin,
                                       &dec->repair_mask);
  if (!success) {
    om_destroy(&M);
    return false;
  }

  int miss_row = 0;
  for (int row = 0; row < dec->num_symbols && miss_row < M.rows; row++) {
    if (bitmask_check(&dec->repair_mask, row))
      continue;
    nanorq_decode_write_esi(rq, io, sbn, row, om_R(M, miss_row), M.cols);
    bitmask_set(&dec->repair_mask, row);
    miss_row++;
  }
  om_destroy(&M);

  return (nanorq_num_missing(rq, sbn) == 0);
}
