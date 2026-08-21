#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <nanorq.h>
#include <rqtables.h>

void dump_esi(nanorq *rq, struct ioctx *myio, FILE *oh, uint8_t sbn,
              uint32_t esi) {
  uint32_t tag = nanorq_tag(sbn, esi);
  size_t packet_size = nanorq_symbol_size(rq);
  uint8_t data[packet_size];
  memset(data, 0, packet_size);
  uint64_t written = nanorq_encode(rq, (void *)data, esi, sbn, myio);

  if (written != packet_size) {
    fprintf(stdout, "failed to encode packet data for sbn %d esi %d.", sbn,
            esi);
    abort();
  } else {
    if (fwrite(&tag, 1, sizeof(tag), oh) != sizeof(tag) ||
        fwrite(data, 1, packet_size, oh) != packet_size) {
      fprintf(stderr, "failed to write encoded packet\n");
      exit(EXIT_FAILURE);
    }
  }
}

void dump_block(nanorq *rq, struct ioctx *myio, FILE *oh, uint8_t sbn,
                float expected_loss, int overhead) {

  uint32_t num_esi = nanorq_block_symbols(rq, sbn);
  int num_dropped = 0, num_rep = 0;
  for (uint32_t esi = 0; esi < num_esi; esi++) {
    float dropped = ((float)(rand()) / (float)RAND_MAX) * (float)100.0;
    float drop_prob = expected_loss;
    if (dropped < drop_prob) {
      num_dropped++;
    } else {
      dump_esi(rq, myio, oh, sbn, esi);
    }
  }
  for (uint32_t esi = num_esi; esi < num_esi + num_dropped + overhead; esi++) {
    dump_esi(rq, myio, oh, sbn, esi);
    num_rep++;
  }
  nanorq_encoder_cleanup(rq, sbn);
  fprintf(stdout, "block %d is %d packets, dropped %d, created %d repair\n",
          sbn, num_esi, num_dropped, num_rep);
}

void usage(char *prog) {
  fprintf(stdout, "usage:\n%s <filename> <packet_size>\n", prog);
  exit(1);
}

int main(int argc, char *argv[]) {
  if (argc < 3)
    usage(argv[0]);

  char *infile = argv[1];
  float expected_loss = 6.0;
  int overhead = 5;
  if (argc >= 4) {
    char *loss_end = NULL;
    errno = 0;
    expected_loss = strtof(argv[3], &loss_end);
    if (errno != 0 || !loss_end || *loss_end != '\0')
      expected_loss = NAN;
  }
  if (argc >= 5) {
    char *overhead_end = NULL;
    errno = 0;
    long overhead_arg = strtol(argv[4], &overhead_end, 10);
    if (errno != 0 || !overhead_end || *overhead_end != '\0' ||
        overhead_arg < 0 || overhead_arg > INT_MAX)
      overhead = -1;
    else
      overhead = (int)overhead_arg;
  }
  if (!isfinite(expected_loss) || expected_loss < 0.0f ||
      expected_loss > 100.0f || overhead < 0 ||
      overhead > ((1 << 24) - 1 - 2 * K_max)) {
    fprintf(stderr,
            "loss must be 0..100 and overhead must fit the 24-bit ESI range\n");
    return EXIT_FAILURE;
  }
  struct ioctx *myio = ioctx_from_file(infile, 1);
  if (!myio) {
    fprintf(stdout, "couldnt access file %s\n", infile);
    return -1;
  }

  size_t filesize = myio->size(myio);

  // determine chunks, symbol size, memory usage from size
  errno = 0;
  char *packet_end = NULL;
  unsigned long packet_arg = strtoul(argv[2], &packet_end, 10);
  if (errno != 0 || !packet_end || *packet_end != '\0' || packet_arg == 0 ||
      packet_arg > UINT16_MAX) {
    fprintf(stderr, "packet_size must be an integer from 1 to %u\n",
            UINT16_MAX);
    myio->destroy(myio);
    return EXIT_FAILURE;
  }
  size_t packet_size = packet_arg; // T
  uint8_t align = 8;

  srand((unsigned int)time(0));

  nanorq *rq = nanorq_encoder_new(filesize, packet_size, align);

  if (rq == NULL) {
    fprintf(stdout, "Could not initialize encoder.\n");
    myio->destroy(myio);
    return -1;
  }

  int num_sbn = nanorq_blocks(rq);
  uint64_t oti_common = nanorq_oti_common(rq);
  uint32_t oti_scheme = nanorq_oti_scheme_specific(rq);
  FILE *oh = fopen("data.rq", "w+b");
  if (!oh ||
      fwrite(&oti_common, 1, sizeof(oti_common), oh) != sizeof(oti_common) ||
      fwrite(&oti_scheme, 1, sizeof(oti_scheme), oh) != sizeof(oti_scheme)) {
    fprintf(stderr, "could not create data.rq\n");
    if (oh)
      fclose(oh);
    nanorq_free(rq);
    myio->destroy(myio);
    return EXIT_FAILURE;
  }
  for (int b = 0; b < num_sbn; b++) {
    dump_block(rq, myio, oh, b, expected_loss, overhead);
  }
  fclose(oh);

  nanorq_free(rq);
  myio->destroy(myio);

  return 0;
}
