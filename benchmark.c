#include <assert.h>
#include <endian.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include <nanorq.h>

#define NUM_SBN 10

#include "kvec.h"

struct sym {
  uint32_t fid;
  uint8_t *data;
};

typedef kvec_t(struct sym) symvec;

uint64_t ms() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (tv.tv_sec * (uint64_t)1000000 + tv.tv_usec) / 1000;
}

void random_bytes(uint8_t *buf, size_t len) {
  for (int i = 0; i < len; i++) {
    buf[i] = rand();
  }
}

void dump_esi(nanorq *rq, struct ioctx *myio, uint8_t sbn, uint32_t esi,
              symvec *packets) {
  uint16_t packet_size = nanorq_symbol_size(rq);
  uint8_t *data = malloc(packet_size);
  uint64_t written = nanorq_encode(rq, (void *)data, esi, sbn, myio);

  if (written != packet_size) {
    free(data);
    fprintf(stderr, "failed to encode packet data for sbn %d esi %d.", sbn,
            esi);
    abort();
  } else {
    uint32_t fid = nanorq_fid(sbn, esi);
    struct sym s = {.fid = fid, .data = data};
    kv_push(struct sym, *packets, s);
  }
}

void dump_block(nanorq *rq, struct ioctx *myio, uint8_t sbn, symvec *packets,
                float overhead_pct) {
  float expected_loss = 6.0;
  uint32_t num_esi = nanorq_block_symbols(rq, sbn);
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

  nanorq_encode_cleanup(rq, sbn);
}

void usage(char *prog) {
  fprintf(stderr, "usage:\n%s <packet_size> <num_packets> <overhead_pct>\n", prog);
  exit(1);
}

int encode(size_t len, uint16_t packet_size, uint16_t num_packets,
           float overhead_pct, struct ioctx *myio, symvec *packets,
           uint64_t *oti_common, uint32_t *oti_scheme) {
  nanorq *rq = nanorq_encoder_new_ex(len, packet_size, num_packets, 0, 8);

  if (rq == NULL) {
    fprintf(stderr, "Coud not initialize encoder.\n");
    return -1;
  }

  *oti_common = nanorq_oti_common(rq);
  *oti_scheme = nanorq_oti_scheme_specific(rq);

  uint8_t num_sbn = nanorq_blocks(rq);
  for (uint8_t b = 0; b < num_sbn; b++) {
    nanorq_generate_symbols(rq, b, myio);
  }

  for (uint8_t sbn = 0; sbn < num_sbn; sbn++) {
    dump_block(rq, myio, sbn, packets, overhead_pct);
  }
  nanorq_free(rq);
  return 0;
}

int decode(uint64_t oti_common, uint32_t oti_scheme, struct ioctx *myio,
           symvec *packets) {

  nanorq *rq = nanorq_decoder_new(oti_common, oti_scheme);
  if (rq == NULL) {
    fprintf(stderr, "Coud not initialize decoder.\n");
    return -1;
  }

  uint8_t num_sbn = nanorq_blocks(rq);
  uint64_t written = 0;

  for (int i = 0; i < kv_size(*packets); i++) {
    struct sym s = kv_A(*packets, i);
    if (!nanorq_decoder_add_symbol(rq, (void *)s.data, s.fid)) {
      fprintf(stderr, "adding symbol %d failed.\n", s.fid);
      abort();
    }
  }
  for (int sbn = 0; sbn < num_sbn; sbn++) {
    written = nanorq_decode_block(rq, myio, sbn);
    if (written == 0) {
      fprintf(stderr, "decode of sbn %d failed.\n", sbn);
    }
    nanorq_decode_cleanup(rq, sbn);
  }
  nanorq_free(rq);
  return 0;
}

int run(uint16_t num_packets, uint16_t packet_size, float overhead_pct) {
  uint64_t t0, elapsed;
  size_t objsize = num_packets * packet_size * NUM_SBN;
  uint64_t oti_common;
  uint32_t oti_scheme;
  struct ioctx *myio;

  size_t sz = num_packets * packet_size * NUM_SBN;
  uint8_t *in = malloc(sz);
  uint8_t *out = malloc(sz);
  random_bytes(in, sz);

  myio = ioctx_from_mem(in, objsize);
  if (!myio) {
    fprintf(stderr, "couldnt access mem at %p\n", in);
    return -1;
  }

  symvec packets;
  kv_init(packets);

  // encode
  t0 = ms();
  encode(objsize, packet_size, num_packets, overhead_pct, myio, &packets,
         &oti_common, &oti_scheme);
  elapsed = ms() - t0;

  fprintf(stdout,
          "ENCODE | Symbol size: %d, symbol count = %d, encoded %.2f MB in "
          "%5.3fsecs, "
          "throughput: "
          "%6.1fMbit/s \n",
          packet_size, num_packets, 1.0 * objsize / (1024 * 1024),
          elapsed / 1000.0, (8.0 * objsize / (1.024 * 1024 * elapsed)));

  myio->destroy(myio);

  // decode
  myio = ioctx_from_mem(out, objsize);
  if (!myio) {
    fprintf(stderr, "couldnt access mem at %p\n", out);
    return -1;
  }

  t0 = ms();
  decode(oti_common, oti_scheme, myio, &packets);
  elapsed = ms() - t0;

  fprintf(stdout,
          "DECODE | Symbol size: %d, symbol count = %d, decoded %.2f MB in "
          "%5.3fsecs using %3.1f%% overhead, throughput: %6.1fMbit/s \n",
          packet_size, num_packets, 1.0 * objsize / (1024 * 1024),
          elapsed / 1000.0, overhead_pct,
          (8.0 * objsize / (1.024 * 1024 * elapsed)));

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
  uint16_t packet_size = strtol(argv[1], NULL, 10); // T
  uint16_t num_packets = strtol(argv[2], NULL, 10); // K
  float overhead_pct = strtof(argv[3], NULL);       // overhead pct

  run(num_packets, packet_size, overhead_pct);

  return 0;
}
