#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <err.h>

#include <nanorq.h>
#include "nanorq_ops.h"

#define NSECS 1000000000L
static int64_t timespec_diff(const struct timespec *t1, const struct timespec *t0)
{
    int64_t sec = (int64_t)t1->tv_sec - (int64_t)t0->tv_sec;
    int64_t nsec = (int64_t)t1->tv_nsec - (int64_t)t0->tv_nsec;
    return (int64_t)((int64_t)sec * NSECS) + nsec;
}

static struct timespec now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts;
}

static void prepare_data_mat(uint8_t *D, const char *f, u32 rows, u32 T, u32 K, u32 SH, u32 stride)
{
    FILE *h = fopen(f, "r");
    if (!h)
        return;
    for (u32 row = 0; row < SH; row++)
        for (u32 col = 0; col < T; col++)
            D[row * stride + col] = 0;
    for (u32 row = SH, off = 0; row < K + SH; row++, off += T) {
        u8 *ptr = D + row * stride;
        size_t got = 0;
        if (!feof(h) && !ferror(h))
            got = fread(ptr, 1, T, h);
        for (u32 col = got; col < T; col++)
            D[row * stride + col] = 0x5a;
    }
    fclose(h);
    for (u32 row = K + SH; row < rows; row++)
        for (u32 col = 0; col < T; col++)
            D[row * stride + col] = 0;
}

int main(int argc, char *argv[])
{
    nanorq rq;
    uint8_t *D;

    struct timespec calc_start, calc_end, ops_start, ops_end;
    if (argc < 5)
        errx(EXIT_FAILURE, "usage: %s <K> <T> <drops> <file>\n", argv[0]);

    int K = strtol(argv[1], NULL, 10);
    if (K < 10 || K > 56403)
        errx(EXIT_FAILURE, "K [packet_count] out of range (10,56403)\n");

    int T = strtol(argv[2], NULL, 10);
    if (T < 8 || T > 65536)
        errx(EXIT_FAILURE, "T [packet_size] out of range (8,65536)\n");

    int drops = strtol(argv[3], NULL, 10);
    if (drops < 1 || drops > K)
        errx(EXIT_FAILURE, "drops should be > 0 and <= K\n");

    if (!nanorq_encoder_new(K, 0, &rq))
        errx(EXIT_FAILURE, "failed to init codec\n");

    size_t sched_bytes = ops_estimate_schedule_bytes(K);
    schedule S_enc, S_dec;
    schedule_init(&S_enc, malloc(sched_bytes), sched_bytes);
    schedule_init(&S_dec, malloc(sched_bytes), sched_bytes);

    calc_start = now();
    size_t prep_len = nanorq_calculate_prepare_memory(&rq);
    uint8_t *prep_mem = malloc(prep_len);
    if (!nanorq_prepare(&rq, prep_mem, prep_len))
        errx(1, "OOM prepare");

    size_t work_len = nanorq_calculate_work_memory(&rq);
    uint8_t *work_mem = malloc(work_len);
    nanorq_set_op_callback(&rq, &S_enc, ops_push);
    if (!nanorq_precalculate(&rq, work_mem, work_len))
        errx(EXIT_FAILURE, "encoder precalculate failed\n");

    u32 rows = nanorq_get_pc_rows(&rq);
    u32 SH = nanorq_get_pc_genc_offset(&rq);
    uint32_t stride = nanorq_recommended_stride(T);
    u32 mem = rows * stride;

    D = obl_alloc(rows, stride, nanorq_oblas.align_size);
    prepare_data_mat(D, argv[4], rows, T, K, SH, stride);

    /* backup d matrix before ops_run */
    u8 *D_backup = malloc((size_t)mem);
    for (u32 row = 0; row < rows; row++) {
        u8 *src = (D + row * stride);
        u8 *dst = D_backup + row * stride;
        for (u32 col = 0; col < PAD(T); col++) {
            dst[col] = src[col];
        }
    }

    u8 *orig_pkts = NULL;
    if (drops > 0)
        orig_pkts = calloc(drops, (size_t)stride);
    for (u32 rp = 0; rp < drops; rp++) {
        u8 *src = (D + (SH + rp) * stride);
        for (u32 i = 0; i < PAD(T); i++) {
            orig_pkts[rp * stride + i] = src[i];
        }
    }

    /* encode to get repair packets */
    ops_run(&rq, D, stride, &S_enc);

    /* save repair packets */
    u8 *repair_pkts = NULL;
    if (drops > 0)
        repair_pkts = calloc(drops, (size_t)stride);
    for (u32 rp = 0; rp < drops; rp++) {
        ops_mix(&rq, D, stride, K + rp, repair_pkts + rp * stride);
    }

    calc_start = now();
    /* prepare fresh decoder matrix */
    if (!nanorq_prepare(&rq, prep_mem, prep_len))
        errx(1, "OOM prepare");

    /* simulate drop and substitute repair data */
    for (u32 rp = 0; rp < drops; rp++) {
        /* drop first 'drops' source packets */
        u32 dropped_esi = rp;
        u32 repair_esi = K + rp;
        nanorq_replace_symbol(&rq, dropped_esi, repair_esi);

        /* insert repair packet data into d matrix at dropped row */
        nanorq_place_symbol(&rq, D, stride, dropped_esi, repair_pkts + rp * stride, PAD(T));
    }

    nanorq_patch_matrix(&rq);

    free(work_mem);
    work_len = nanorq_calculate_work_memory(&rq);
    work_mem = malloc(work_len);

    nanorq_set_op_callback(&rq, &S_dec, ops_push);
    if (!nanorq_precalculate(&rq, work_mem, work_len))
        errx(EXIT_FAILURE, "decoder precalculate failed\n");
    calc_end = now();

    /* restore d matrix before decoding */
    for (u32 row = 0; row < rows; row++) {
        u8 *dst = (D + row * stride);
        u8 *src = D_backup + row * stride;
        for (u32 i = 0; i < PAD(T); i++) {
            dst[i] = src[i];
        }
    }
    free(D_backup);

    /* re-inject repair packet data into d matrix at dropped rows */
    for (u32 rp = 0; rp < drops; rp++) {
        u32 dropped_esi = rp;
        nanorq_place_symbol(&rq, D, stride, dropped_esi, repair_pkts + rp * stride, PAD(T));
    }

    ops_start = now();
    ops_run(&rq, D, stride, &S_dec);
    ops_end = now();

    int diffs = 0;
    u8 *reppkt = calloc(1, (size_t)stride);
    for (u32 rp = 0; rp < drops; rp++) {
        ops_mix(&rq, D, stride, rp, reppkt);
        u8 *expected = orig_pkts + rp * stride;
        for (u32 i = 0; i < T; i++) {
            if (reppkt[i] != expected[i]) {
                diffs++;
            }
        }
    }
    free(reppkt);

    if (diffs > 0) {
        fprintf(stderr, "Decode failed with %d byte mismatches!\n", diffs);
        return 1;
    }

    double calc_time = timespec_diff(&calc_end, &calc_start) / (double)NSECS;
    double ops_time = timespec_diff(&ops_end, &ops_start) / (double)NSECS;

    fprintf(stderr, "calc: %.6fs ops: %.6fs bytes: %u\n", calc_time, ops_time, T * K);

    free(prep_mem);
    free(work_mem);
    free(D);
    if (orig_pkts)
        free(orig_pkts);
    if (repair_pkts)
        free(repair_pkts);
    free(S_enc.ops.a);
    free(S_dec.ops.a);

    return 0;
}
