
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <nanorq.h>

void usage(char *prog) {
  fprintf(stdout, "usage:\n%s <output_filename>\n", prog);
  exit(1);
}

int main(int argc, char *argv[]) {

  if (argc < 2)
    usage(argv[0]);

  char *outfile = argv[1];
  FILE *ih = fopen("data.rq", "rb");
  if (!ih) {
    fprintf(stderr, "could not open data.rq\n");
    return EXIT_FAILURE;
  }
  uint32_t oti_scheme;
  uint64_t oti_common;

  if (fread(&oti_common, 1, sizeof(oti_common), ih) != sizeof(oti_common) ||
      fread(&oti_scheme, 1, sizeof(oti_scheme), ih) != sizeof(oti_scheme)) {
    fprintf(stderr, "data.rq has a truncated header\n");
    fclose(ih);
    return EXIT_FAILURE;
  }

  nanorq *rq = nanorq_decoder_new(oti_common, oti_scheme);
  if (rq == NULL) {
    fprintf(stdout, "Could not initialize decoder.\n");
    fclose(ih);
    return -1;
  }

  struct ioctx *myio = ioctx_from_file(outfile, IOCTX_MODE_WRITE);
  if (!myio) {
    fprintf(stderr, "could not create output file %s\n", outfile);
    nanorq_free(rq);
    fclose(ih);
    return EXIT_FAILURE;
  }

  int num_sbn = nanorq_blocks(rq);
  uint32_t tag;
  size_t packet_size = nanorq_symbol_size(rq);
  uint8_t packet[packet_size];
  bool failed = false;
  size_t got_tag;
  while ((got_tag = fread(&tag, 1, sizeof(tag), ih)) == sizeof(tag)) {
    if (fread(packet, 1, packet_size, ih) != packet_size) {
      fprintf(stderr, "data.rq has a truncated packet\n");
      failed = true;
      break;
    }
    if (NANORQ_SYM_ERR ==
        nanorq_decoder_add_symbol(rq, (void *)packet, tag, myio)) {
      fprintf(stdout, "adding symbol %d failed.\n", tag);
      failed = true;
      break;
    }
  }
  if (got_tag != 0) {
    fprintf(stderr, "data.rq has a truncated packet tag\n");
    failed = true;
  }
  if (!feof(ih) && ferror(ih)) {
    fprintf(stderr, "failed while reading data.rq\n");
    failed = true;
  }
  for (int sbn = 0; !failed && sbn < num_sbn; sbn++) {
    fprintf(stdout, "block %d is %d packets, lost %d, have %d repair\n", sbn,
            (unsigned)nanorq_block_symbols(rq, sbn),
            (unsigned)nanorq_num_missing(rq, sbn),
            (unsigned)nanorq_num_repair(rq, sbn));
    if (!nanorq_repair_block(rq, myio, sbn)) {
      fprintf(stdout, "decode of sbn %d failed.\n", sbn);
      failed = true;
    }
    nanorq_encoder_cleanup(rq, sbn);
  }
  fclose(ih);
  nanorq_free(rq);
  myio->destroy(myio);

  return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
