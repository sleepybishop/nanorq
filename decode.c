
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <nanorq.h>

void usage(char *prog) {
  fprintf(stdout, "usage:\n%s <filename> <packet_size>\n", prog);
  exit(1);
}

int main(int argc, char *argv[]) {

  if (argc < 2)
    usage(argv[0]);

  char *outfile = argv[1];
  struct ioctx *myio = ioctx_from_file(outfile, 0);
  if (!myio) {
    fprintf(stdout, "couldnt access file %s\n", outfile);
    return -1;
  }

  FILE *ih = fopen("data.rq", "r");
  uint32_t oti_scheme;
  uint64_t oti_common;

  fread(&oti_common, 1, sizeof(oti_common), ih);
  fread(&oti_scheme, 1, sizeof(oti_scheme), ih);

  nanorq *rq = nanorq_decoder_new(oti_common, oti_scheme);
  if (rq == NULL) {
    fprintf(stdout, "Could not initialize decoder.\n");
    return -1;
  }

  uint8_t num_sbn = nanorq_blocks(rq);
  uint32_t tag;
  size_t packet_size = nanorq_symbol_size(rq);
  uint8_t packet[packet_size];
  while (fread(&tag, 1, sizeof(tag), ih)) {
    fread(packet, packet_size, 1, ih);
    if (!nanorq_decoder_add_symbol(rq, (void *)packet, tag, myio)) {
      fprintf(stdout, "adding symbol %d failed.\n", tag);
      abort();
    }
  }
  for (int sbn = 0; sbn < num_sbn; sbn++) {
    fprintf(stdout, "block %d is %d packets, lost %d, have %d repair\n", sbn,
            (unsigned)nanorq_block_symbols(rq, sbn),
            (unsigned)nanorq_num_missing(rq, sbn),
            (unsigned)nanorq_num_repair(rq, sbn));
    if (!nanorq_repair_block(rq, myio, sbn)) {
      fprintf(stdout, "decode of sbn %d failed.\n", sbn);
    }
    nanorq_encoder_cleanup(rq, sbn);
  }
  fclose(ih);
  nanorq_free(rq);
  myio->destroy(myio);

  return 0;
}
