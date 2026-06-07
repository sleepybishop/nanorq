#define _GNU_SOURCE
#define _POSIX_C_SOURCE 199309L
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "nanorq_core.h"
#include "nanorq_ops.h"

#define TEST_BYTES 256 * 1024 * 1024

uint64_t usecs() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (tv.tv_sec * (uint64_t)1000000 + tv.tv_usec);
}

void usage(char *prog) {
  fprintf(stderr,
          "usage:\n%s <packet_size> <num_packets> <overhead_pct> "
          "[<expected_loss>]\n",
          prog);
  exit(1);
}

int main(int argc, char *argv[]) {
  if (argc < 4)
    usage(argv[0]);

  srand((unsigned int)time(0));

  size_t T = strtol(argv[1], NULL, 10);       /* packet size */
  size_t K = strtol(argv[2], NULL, 10);       /* num packets */
  float overhead_pct = strtof(argv[3], NULL); /* overhead pct */
  float expected_loss = 6.0;
  if (argc >= 5)
    expected_loss = strtof(argv[4], NULL);

  nanorq_core rq;
  if (!nanorq_core_encoder_new(K, 0, &rq)) {
    fprintf(stderr, "failed to init encoder\n");
    return 1;
  }

  size_t prep_len = nanorq_core_calculate_prepare_memory(&rq);
  uint8_t *prep_mem = (uint8_t *)malloc(prep_len);
  size_t work_len = nanorq_core_calculate_work_memory(&rq);
  uint8_t *work_mem = (uint8_t *)malloc(work_len);

  schedule S_enc = {};
  size_t sched_bytes = ops_estimate_schedule_bytes(K);
  S_enc.ops.a = (sched_op *)malloc(sched_bytes);
  S_enc.ops.m = sched_bytes / sizeof(sched_op);
  nanorq_core_set_op_callback(&rq, &S_enc, ops_push);

  if (!nanorq_core_prepare(&rq, prep_mem, prep_len)) {
    fprintf(stderr, "nanorq_core_prepare failed\n");
    return 1;
  }

  u32 rows = nanorq_core_get_pc_rows(&rq);
  u32 SH = nanorq_core_get_pc_genc_offset(&rq);
  u32 stride = nanorq_core_recommended_stride(T);
  u32 mem = rows * stride;
  uint8_t *D = (uint8_t *)obl_alloc(rows, stride, nanorq_oblas.align_size);
  uint8_t *D_backup = (uint8_t *)malloc(mem);

  memset(D_backup, 0, mem);
  for (u32 row = SH; row < SH + K; row++) {
    for (u32 col = 0; col < T; col++) {
      D_backup[row * stride + col] = rand() & 0xff;
    }
  }

  /* encode without precalc */
  uint64_t bytes = 0;
  uint64_t t0 = usecs();
  while (bytes < TEST_BYTES) {
    if (!nanorq_core_prepare(&rq, prep_mem, prep_len)) {
      fprintf(stderr, "nanorq_core_prepare failed\n");
      exit(1);
    }
    nanorq_core_set_op_callback(&rq, &S_enc, ops_push);
    S_enc.ops.n = 0;
    S_enc.cpidx = 0;
    if (!nanorq_core_precalculate(&rq, work_mem, work_len)) {
      fprintf(stderr, "nanorq_core_precalculate failed\n");
      exit(1);
    }

    uint64_t t_mid0 = usecs();
    memcpy(D, D_backup, mem);
    uint64_t t_mid1 = usecs();
    t0 += (t_mid1 - t_mid0);

    ops_run(&rq, D, stride, &S_enc);
    bytes += T * K;
  }
  double elapsed_encode = (usecs() - t0) / 1000000.0;

  /* encode with precalc */
  if (!nanorq_core_prepare(&rq, prep_mem, prep_len)) {
    fprintf(stderr, "nanorq_core_prepare failed\n");
    exit(1);
  }
  nanorq_core_set_op_callback(&rq, &S_enc, ops_push);
  S_enc.ops.n = 0;
  S_enc.cpidx = 0;
  if (!nanorq_core_precalculate(&rq, work_mem, work_len)) {
    fprintf(stderr, "nanorq_core_precalculate failed\n");
    exit(1);
  }

  bytes = 0;
  t0 = usecs();
  while (bytes < TEST_BYTES) {
    uint64_t t_mid0 = usecs();
    memcpy(D, D_backup, mem);
    uint64_t t_mid1 = usecs();
    t0 += (t_mid1 - t_mid0);

    ops_run(&rq, D, stride, &S_enc);
    bytes += T * K;
  }
  double elapsed_precalc = (usecs() - t0) / 1000000.0;

  /* decode with no extra overhead */
  /* encode to compute intermediate symbols */
  memcpy(D, D_backup, mem);
  ops_run(&rq, D, stride, &S_enc);

  int drops = (int)(K * expected_loss / 100.0);
  if (drops <= 0)
    drops = 1;

  u32 *dropped_esi = (u32 *)malloc(sizeof(u32) * drops);
  for (u32 rp = 0; rp < (u32)drops; rp++) {
    dropped_esi[rp] = rp;
  }

  uint8_t *repair_pkts = (uint8_t *)malloc(drops * stride);
  for (u32 rp = 0; rp < (u32)drops; rp++) {
    ops_mix(&rq, D, stride, K + rp, repair_pkts + rp * stride);
  }

  uint8_t *orig_pkts = (uint8_t *)malloc(drops * stride);
  for (u32 rp = 0; rp < (u32)drops; rp++) {
    memcpy(orig_pkts + rp * stride, D_backup + (SH + dropped_esi[rp]) * stride,
           PAD(T));
  }

  nanorq_core dec;
  if (!nanorq_core_encoder_new(K, 0, &dec)) {
    fprintf(stderr, "failed to init decoder\n");
    return 1;
  }

  size_t dec_prep_len = nanorq_core_calculate_prepare_memory(&dec);
  uint8_t *dec_prep_mem = (uint8_t *)malloc(dec_prep_len);
  size_t dec_work_len = nanorq_core_calculate_work_memory(&dec);
  uint8_t *dec_work_mem = (uint8_t *)malloc(dec_work_len);

  schedule S_dec = {};
  S_dec.ops.a = (sched_op *)malloc(sched_bytes);
  S_dec.ops.m = sched_bytes / sizeof(sched_op);
  nanorq_core_set_op_callback(&dec, &S_dec, ops_push);

  uint8_t *D_dec = (uint8_t *)obl_alloc(rows, stride, nanorq_oblas.align_size);
  uint8_t *decoded_pkts = (uint8_t *)malloc(drops * stride);

  bytes = 0;
  t0 = usecs();
  double t_prep_dec = 0.0, t_precalc_dec = 0.0, t_ops_run_dec = 0.0,
         t_ops_mix_dec = 0.0;
  while (bytes < TEST_BYTES) {
    uint64_t t_start = usecs();
    if (!nanorq_core_prepare(&dec, dec_prep_mem, dec_prep_len)) {
      fprintf(stderr, "nanorq_core_prepare failed\n");
      exit(1);
    }
    nanorq_core_set_op_callback(&dec, &S_dec, ops_push);
    for (u32 rp = 0; rp < (u32)drops; rp++) {
      nanorq_core_replace_symbol(&dec, dropped_esi[rp], K + rp);
    }
    nanorq_core_patch_matrix(&dec);
    t_prep_dec += (usecs() - t_start) / 1000000.0;

    t_start = usecs();
    S_dec.ops.n = 0;
    S_dec.cpidx = 0;
    if (!nanorq_core_precalculate(&dec, dec_work_mem, dec_work_len)) {
      fprintf(stderr, "decoder precalculate failed\n");
      exit(1);
    }
    t_precalc_dec += (usecs() - t_start) / 1000000.0;

    uint64_t t_mid0 = usecs();
    memcpy(D_dec, D_backup, mem);
    for (u32 rp = 0; rp < (u32)drops; rp++) {
      memcpy(D_dec + (SH + dropped_esi[rp]) * stride, repair_pkts + rp * stride,
             PAD(T));
    }
    uint64_t t_mid1 = usecs();
    t0 += (t_mid1 - t_mid0);

    t_start = usecs();
    ops_run(&dec, D_dec, stride, &S_dec);
    t_ops_run_dec += (usecs() - t_start) / 1000000.0;

    t_start = usecs();
    for (u32 rp = 0; rp < (u32)drops; rp++) {
      ops_mix(&dec, D_dec, stride, dropped_esi[rp], decoded_pkts + rp * stride);
    }
    t_ops_mix_dec += (usecs() - t_start) / 1000000.0;

    bytes += T * K;
  }
  double elapsed_decode = (usecs() - t0) / 1000000.0;

  if (K == 50000 && 0) {
    fprintf(stderr,
            "\n[K=50000 decode Breakdown] prep: %.4fs, precalc: %.4fs, "
            "ops_run: %.4fs, ops_mix: %.4fs\n",
            t_prep_dec, t_precalc_dec, t_ops_run_dec, t_ops_mix_dec);
  }

  for (u32 rp = 0; rp < (u32)drops; rp++) {
    if (memcmp(decoded_pkts + rp * stride, orig_pkts + rp * stride, T) != 0) {
      fprintf(stderr, "verification failed for drop %u!\n", rp);
      fprintf(stderr, "Expected: ");
      for (u32 i = 0; i < 16 && i < T; i++)
        fprintf(stderr, "%02x ", orig_pkts[rp * stride + i]);
      fprintf(stderr, "\nDecoded:  ");
      for (u32 i = 0; i < 16 && i < T; i++)
        fprintf(stderr, "%02x ", decoded_pkts[rp * stride + i]);
      fprintf(stderr, "\n");
      exit(1);
    }
  }

  /* decode with extra overhead */
  int overhead = (int)(K * overhead_pct) / 100;
  if (overhead <= 0)
    overhead = 1;

  free(repair_pkts);
  repair_pkts = (uint8_t *)malloc((drops + overhead) * stride);
  for (u32 rp = 0; rp < (u32)(drops + overhead); rp++) {
    ops_mix(&rq, D, stride, K + rp, repair_pkts + rp * stride);
  }

  nanorq_core dec_oh;
  if (!nanorq_core_encoder_new(K, overhead, &dec_oh)) {
    fprintf(stderr, "failed to init decoder_oh\n");
    return 1;
  }

  size_t dec_oh_prep_len = nanorq_core_calculate_prepare_memory(&dec_oh);
  uint8_t *dec_oh_prep_mem = (uint8_t *)malloc(dec_oh_prep_len);
  size_t dec_oh_work_len = nanorq_core_calculate_work_memory(&dec_oh);
  uint8_t *dec_oh_work_mem = (uint8_t *)malloc(dec_oh_work_len);

  schedule S_dec_oh = {};
  size_t dec_oh_sched_bytes = ops_estimate_schedule_bytes(K + overhead);
  S_dec_oh.ops.a = (sched_op *)malloc(dec_oh_sched_bytes);
  S_dec_oh.ops.m = dec_oh_sched_bytes / sizeof(sched_op);
  nanorq_core_set_op_callback(&dec_oh, &S_dec_oh, ops_push);

  if (!nanorq_core_prepare(&dec_oh, dec_oh_prep_mem, dec_oh_prep_len)) {
    fprintf(stderr, "nanorq_core_prepare failed\n");
    return 1;
  }

  u32 rows_oh = nanorq_core_get_pc_rows(&dec_oh);
  u32 mem_oh = rows_oh * stride;
  uint8_t *D_dec_oh =
      (uint8_t *)obl_alloc(rows_oh, stride, nanorq_oblas.align_size);

  double t_prep = 0.0, t_precalc = 0.0, t_ops_run = 0.0, t_ops_mix = 0.0;
  bytes = 0;
  t0 = usecs();
  while (bytes < TEST_BYTES) {
    uint64_t t_start = usecs();
    if (!nanorq_core_prepare(&dec_oh, dec_oh_prep_mem, dec_oh_prep_len)) {
      fprintf(stderr, "nanorq_core_prepare failed\n");
      exit(1);
    }
    nanorq_core_set_op_callback(&dec_oh, &S_dec_oh, ops_push);
    for (u32 rp = 0; rp < (u32)drops; rp++) {
      nanorq_core_replace_symbol(&dec_oh, dropped_esi[rp], K + overhead + rp);
    }
    for (u32 i = 0; i < (u32)overhead; i++) {
      nanorq_core_replace_symbol(&dec_oh, K + i, K + i);
    }
    nanorq_core_patch_matrix(&dec_oh);
    t_prep += (usecs() - t_start) / 1000000.0;

    t_start = usecs();
    S_dec_oh.ops.n = 0;
    S_dec_oh.cpidx = 0;
    if (!nanorq_core_precalculate(&dec_oh, dec_oh_work_mem, dec_oh_work_len)) {
      fprintf(stderr, "decoder precalculate failed\n");
      exit(1);
    }
    t_precalc += (usecs() - t_start) / 1000000.0;

    uint64_t t_mid0 = usecs();
    memcpy(D_dec_oh, D_backup, mem);
    memset(D_dec_oh + mem, 0, mem_oh - mem);
    for (u32 rp = 0; rp < (u32)drops; rp++) {
      memcpy(D_dec_oh + (SH + dropped_esi[rp]) * stride,
             repair_pkts + (overhead + rp) * stride, PAD(T));
    }
    for (u32 i = 0; i < (u32)overhead; i++) {
      memcpy(D_dec_oh + (SH + K + i) * stride, repair_pkts + i * stride,
             PAD(T));
    }
    uint64_t t_mid1 = usecs();
    t0 += (t_mid1 - t_mid0);

    t_start = usecs();
    ops_run(&dec_oh, D_dec_oh, stride, &S_dec_oh);
    t_ops_run += (usecs() - t_start) / 1000000.0;

    t_start = usecs();
    for (u32 rp = 0; rp < (u32)drops; rp++) {
      ops_mix(&dec_oh, D_dec_oh, stride, dropped_esi[rp],
              decoded_pkts + rp * stride);
    }
    t_ops_mix += (usecs() - t_start) / 1000000.0;

    bytes += T * (K + overhead);
  }
  double elapsed_decode_oh = (usecs() - t0) / 1000000.0;

  if (K == 50000 && 0) {
    fprintf(stderr,
            "\n[K=50000 oh-5 Breakdown] prep: %.4fs, precalc: %.4fs, ops_run: "
            "%.4fs, ops_mix: %.4fs\n",
            t_prep, t_precalc, t_ops_run, t_ops_mix);
    size_t gf256_enc = 0, gf256_dec = 0, gf256_dec_oh = 0;
    for (size_t i = 0; i < S_enc.ops.n; i++)
      if (S_enc.ops.a[i].u > 1)
        gf256_enc++;
    for (size_t i = 0; i < S_dec.ops.n; i++)
      if (S_dec.ops.a[i].u > 1)
        gf256_dec++;
    for (size_t i = 0; i < S_dec_oh.ops.n; i++)
      if (S_dec_oh.ops.a[i].u > 1)
        gf256_dec_oh++;
    fprintf(stderr,
            "Schedule sizes: S_enc=%zu (GF256:%zu), S_dec=%zu (GF256:%zu), "
            "S_dec_oh=%zu (GF256:%zu)\n",
            S_enc.ops.n, gf256_enc, S_dec.ops.n, gf256_dec, S_dec_oh.ops.n,
            gf256_dec_oh);
  }

  for (u32 rp = 0; rp < (u32)drops; rp++) {
    if (memcmp(decoded_pkts + rp * stride, orig_pkts + rp * stride, T) != 0) {
      fprintf(stderr, "verification for decode-oh failed for drop %u!\n", rp);
      exit(1);
    }
  }

  fprintf(stdout, "%10d %10.1f %10.1f %10.1f %10.1f\n", (int)K,
          (8.0 * TEST_BYTES / (1024 * 1024 * elapsed_encode)),
          (8.0 * TEST_BYTES / (1024 * 1024 * elapsed_precalc)),
          (8.0 * TEST_BYTES / (1024 * 1024 * elapsed_decode)),
          (8.0 * TEST_BYTES / (1024 * 1024 * elapsed_decode_oh)));

  /* cleanup */
  free(prep_mem);
  free(work_mem);
  free(S_enc.ops.a);
  free(D);
  free(D_backup);

  free(dropped_esi);
  free(repair_pkts);
  free(orig_pkts);

  free(dec_prep_mem);
  free(dec_work_mem);
  free(S_dec.ops.a);
  free(D_dec);
  free(decoded_pkts);

  free(dec_oh_prep_mem);
  free(dec_oh_work_mem);
  free(S_dec_oh.ops.a);
  free(D_dec_oh);

  return 0;
}
