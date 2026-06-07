#include "../include/io.h"
#include "../include/nanorq.h"
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SYMBOL_SIZE 1024
#define NUM_SYMBOLS 1000
#define PAYLOAD_SIZE (NUM_SYMBOLS * SYMBOL_SIZE)
#define NUM_DEVICES 3

typedef struct {
  int id;
  int join_esi;
  float loss_rate;
  int received_count;
  nanorq *decoder;
  struct ioctx *dec_io;
  uint8_t *decoded_payload;
  bool done;
} Device;

int main() {
  srand(time(NULL));

  /* payload */
  uint8_t *firmware = malloc(PAYLOAD_SIZE);
  for (int i = 0; i < PAYLOAD_SIZE; i++)
    firmware[i] = rand() % 256;
  struct ioctx *enc_io = ioctx_from_mem(firmware, PAYLOAD_SIZE);

  /* decoder */
  nanorq *encoder =
      nanorq_encoder_new_ex(PAYLOAD_SIZE, SYMBOL_SIZE, NUM_SYMBOLS, 0, 8);
  nanorq_generate_symbols(encoder, 0, enc_io);
  nanorq_precalculate(encoder);

  uint64_t common = nanorq_oti_common(encoder);
  uint32_t specific = nanorq_oti_scheme_specific(encoder);

  /* our receivers */
  Device devices[NUM_DEVICES];
  devices[0] = (Device){.id = 1,
                        .join_esi = 0,
                        .loss_rate = 0.05,
                        .received_count = 0,
                        .done = false};
  devices[1] = (Device){.id = 2,
                        .join_esi = 500,
                        .loss_rate = 0.10,
                        .received_count = 0,
                        .done = false};
  devices[2] = (Device){.id = 3,
                        .join_esi = 2000,
                        .loss_rate = 0.20,
                        .received_count = 0,
                        .done = false};

  for (int i = 0; i < NUM_DEVICES; i++) {
    devices[i].decoder = nanorq_decoder_new(common, specific);
    devices[i].decoded_payload = calloc(1, PAYLOAD_SIZE);
    devices[i].dec_io =
        ioctx_from_mem(devices[i].decoded_payload, PAYLOAD_SIZE);
    printf("Device %d initialized: Joins at ESI=%d, Loss Rate=%.0f%%\n",
           devices[i].id, devices[i].join_esi, devices[i].loss_rate * 100);
  }
  printf("\nEncoder starting carousel broadcast...\n\n");

  uint32_t esi = 0;
  int devices_done = 0;

  while (devices_done < NUM_DEVICES) {
    uint8_t pkt[SYMBOL_SIZE];
    size_t encoded = nanorq_encode(encoder, pkt, esi, 0, enc_io);
    if (encoded == 0) {
      printf("Error generating ESI %u (Max ESI reached?)\n", esi);
      break;
    }

    uint32_t tag = nanorq_tag(0, esi);

    /* receivers enter at different points */
    for (int i = 0; i < NUM_DEVICES; i++) {
      if (devices[i].done)
        continue;
      if (esi < devices[i].join_esi)
        continue;

      double r = (double)rand() / RAND_MAX;
      if (r < devices[i].loss_rate)
        continue;
      int res = nanorq_decoder_add_symbol(devices[i].decoder, pkt, tag,
                                          devices[i].dec_io);
      if (res >= 0) {
        devices[i].received_count++;
      }

      /* attempt decode */
      if (devices[i].received_count >= NUM_SYMBOLS &&
          (devices[i].received_count - NUM_SYMBOLS) % 5 == 0) {

        if (nanorq_repair_block(devices[i].decoder, devices[i].dec_io, 0)) {
          if (memcmp(firmware, devices[i].decoded_payload, PAYLOAD_SIZE) == 0) {
            printf(">>> Device %d successfully verified firmware! (Joined: ESI "
                   "%d, Finished: ESI %d, Collected: %d "
                   "symbols)\n",
                   devices[i].id, devices[i].join_esi, esi,
                   devices[i].received_count);
            devices[i].done = true;
            devices_done++;
          }
        }
      }
    }

    esi++;
    if (esi % 1000 == 0 && devices_done < NUM_DEVICES) {
      printf("Carousel reached ESI %d... %d/%d devices updated.\n", esi,
             devices_done, NUM_DEVICES);
    }
  }

  printf("\nAll devices successfully updated firmware via carousel!\n");

  /* cleanup */
  for (int i = 0; i < NUM_DEVICES; i++) {
    nanorq_free(devices[i].decoder);
    devices[i].dec_io->destroy(devices[i].dec_io);
    free(devices[i].decoded_payload);
  }
  nanorq_free(encoder);
  enc_io->destroy(enc_io);
  free(firmware);

  return 0;
}
