#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include <nanorq.h>

#include "kvec.h"

struct sym {
  uint32_t tag;
  uint8_t *data;
};

typedef kvec_t(struct sym) symvec;

uint64_t usecs() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (tv.tv_sec * (uint64_t)1000000 + tv.tv_usec);
}

void random_bytes(uint8_t *buf, uint64_t len) {
  for (int i = 0; i < len; i++) {
    buf[i] = rand();
  }
}

void dump_esi(nanorq *rq, struct ioctx *myio, int sbn, uint32_t esi,
              symvec *packets) {
  int packet_size = nanorq_symbol_size(rq);
  uint8_t *data = malloc(packet_size);
  uint64_t written = nanorq_encode(rq, (void *)data, esi, sbn, myio);

  if (written != packet_size) {
    free(data);
    fprintf(stderr, "failed to encode packet data for sbn %d esi %d.", sbn,
            esi);
    abort();
  } else {
    uint32_t tag = nanorq_tag(sbn, esi);
    struct sym s = {.tag = tag, .data = data};
    kv_push(struct sym, *packets, s);
  }
}

void dump_block(nanorq *rq, struct ioctx *myio, int sbn, symvec *packets,
                float overhead_pct) {
  float expected_loss = 6.0;
  int num_esi = nanorq_block_symbols(rq, sbn);
  int overhead = (int)(num_esi * overhead_pct) / 100;
  int num_dropped = 0, num_rep = 0;
  for (uint32_t esi = 0; esi < num_esi; esi++) {
    float dropped = ((float)(rand()) / (float)RAND_MAX) * (float)100.0;
    float drop_prob = expected_loss;
    if (dropped < drop_prob) {
      num_dropped++;
    } else {
      dump_esi(rq, myio, sbn, esi, packets);
    }
  }
  for (uint32_t esi = num_esi; esi < num_esi + num_dropped + overhead; esi++) {
    dump_esi(rq, myio, sbn, esi, packets);
    num_rep++;
  }

  nanorq_encoder_cleanup(rq, sbn);
}

void usage(char *prog) {
  fprintf(stderr,
          "usage:\n%s <packet_size> <num_packets> <overhead_pct> "
          "[<precalculate: (0,10)]\n",
          prog);
  exit(1);
}

double encode(uint64_t len, size_t packet_size, size_t num_packets,
              float overhead_pct, struct ioctx *myio, symvec *packets,
              uint64_t *oti_common, uint32_t *oti_scheme, bool precalc) {
  nanorq *rq = nanorq_encoder_new_ex(len, packet_size, num_packets, 0, 8);

  if (rq == NULL) {
    fprintf(stderr, "Could not initialize encoder.\n");
    return -1;
  }

  *oti_common = nanorq_oti_common(rq);
  *oti_scheme = nanorq_oti_scheme_specific(rq);

  if (precalc)
    nanorq_precalculate(rq);
  int num_sbn = nanorq_blocks(rq);
  double elapsed = 0.0;
  uint64_t t0 = usecs();
  for (int b = 0; b < num_sbn; b++) {
    nanorq_generate_symbols(rq, b, myio);
  }
  elapsed = (usecs() - t0) / 1000000.0;

  for (int sbn = 0; sbn < num_sbn; sbn++) {
    dump_block(rq, myio, sbn, packets, overhead_pct);
  }
  nanorq_free(rq);
  return elapsed;
}

double decode(uint64_t oti_common, uint32_t oti_scheme, struct ioctx *myio,
              symvec *packets) {

  nanorq *rq = nanorq_decoder_new(oti_common, oti_scheme);
  if (rq == NULL) {
    fprintf(stderr, "Could not initialize decoder.\n");
    return -1;
  }

  int num_sbn = nanorq_blocks(rq);

  for (int i = 0; i < kv_size(*packets); i++) {
    struct sym s = kv_A(*packets, i);
    if (!nanorq_decoder_add_symbol(rq, (void *)s.data, s.tag, myio)) {
      fprintf(stderr, "adding symbol %d to sbn %d failed.\n",
              s.tag & 0x00ffffff, (s.tag >> 24) & 0xff);
      abort();
    }
  }
  double elapsed = 0.0;
  for (int sbn = 0; sbn < num_sbn; sbn++) {
    uint64_t t0 = usecs();
    if (!nanorq_repair_block(rq, myio, sbn)) {
      fprintf(stderr, "decode of sbn %d failed.\n", sbn);
      abort();
    }
    elapsed += (usecs() - t0) / 1000000.0;
    nanorq_encoder_cleanup(rq, sbn);
  }
  nanorq_free(rq);
  return elapsed;
}

int run(size_t num_packets, size_t packet_size, float overhead_pct,
        bool precalc) {
  double elapsed;
  uint64_t objsize = 160 * 1024 * 1024;
  int num_sbn = objsize / (num_packets * packet_size);
  if (num_sbn > 256)
    num_sbn = 256;
  uint64_t oti_common = 0;
  uint32_t oti_scheme = 0;
  struct ioctx *myio;

  uint64_t sz = num_packets * packet_size * num_sbn;
  uint8_t *in = calloc(1, sz);
  uint8_t *out = calloc(1, sz);
  random_bytes(in, sz);

  myio = ioctx_from_mem(in, sz);
  if (!myio) {
    fprintf(stderr, "couldnt access mem at %p\n", in);
    free(in);
    free(out);
    return -1;
  }

  symvec packets;
  kv_init(packets);

  // encode
  elapsed = encode(sz, packet_size, num_packets, overhead_pct, myio, &packets,
                   &oti_common, &oti_scheme, precalc);

  if (elapsed < 0)
    abort();
  fprintf(stdout,
          "ENCODE | Symbol size: %u, symbol count = %u, encoded %.2f MB in "
          "%5.3fsecs, "
          "throughput: "
          "%6.1fMbit/s (%d sbns)\n",
          (unsigned)packet_size, (unsigned)num_packets,
          1.0 * sz / (1024 * 1024), elapsed,
          (8.0 * sz / (1024 * 1024 * elapsed)), num_sbn);

  myio->destroy(myio);

  // decode
  myio = ioctx_from_mem(out, sz);
  if (!myio) {
    fprintf(stderr, "couldnt access mem at %p\n", out);
    free(in);
    free(out);
    kv_destroy(packets);
    return -1;
  }

  elapsed = decode(oti_common, oti_scheme, myio, &packets);

  if (elapsed < 0)
    abort();
  fprintf(
      stdout,
      "DECODE | Symbol size: %u, symbol count = %u, decoded %.2f MB in "
      "%5.3fsecs using %3.1f%% overhead, throughput: %6.1fMbit/s (%d sbns)\n",
      (unsigned)packet_size, (unsigned)num_packets, 1.0 * sz / (1024 * 1024),
      elapsed, overhead_pct, (8.0 * sz / (1024 * 1024 * elapsed)), num_sbn);

  myio->destroy(myio);
  // verify
  for (int i = 0; i < sz; i++) {
    assert(in[i] == out[i]);
  }

  // cleanup
  if (kv_size(packets) > 0) {
    for (int i = 0; i < kv_size(packets); i++) {
      free(kv_A(packets, i).data);
    }
    kv_destroy(packets);
  }

  free(in);
  free(out);
  return 0;
}

int main(int argc, char *argv[]) {
  if (argc < 4)
    usage(argv[0]);

  srand((unsigned int)time(0));

  // determine chunks, symbol size
  size_t packet_size = strtol(argv[1], NULL, 10); // T
  size_t num_packets = strtol(argv[2], NULL, 10); // K
  float overhead_pct = strtof(argv[3], NULL);     // overhead pct
  bool precalc = ((argc > 4) && (strtod(argv[4], NULL) > 0)) ? true : false;

  run(num_packets, packet_size, overhead_pct, precalc);

  return 0;
}
