#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <err.h>

#include <nanorq.h>
#include "operations.h"

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

static void prepare_data_mat(u8_vec *D, const char *f, u32 rows, u32 T, u32 K, u32 SH)
{
    FILE *h = fopen(f, "r");
    if (!h)
        return;
    for (u32 row = 0; row < SH; row++)
        for (u32 col = 0; col < T; col++)
            uv_E(*D, row, col) = 0;
    for (u32 row = SH, off = 0; row < K + SH; row++, off += T) {
        u8 *ptr = uv_R(*D, row);
        size_t got = fread(ptr, 1, T, h);
        for (u32 col = got; col < T; col++)
            uv_E(*D, row, col) = 0x5a;
    }
    fclose(h);
    for (u32 row = K + SH; row < rows; row++)
        for (u32 col = 0; col < T; col++)
            uv_E(*D, row, col) = 0;
}

int main(int argc, char *argv[])
{
    int ok = 0;
    nanorq rq;
    schedule S = {};
    u8_vec D;

    struct timespec calc_start, calc_end, ops_start, ops_end;
    if (argc < 5)
        errx(EXIT_FAILURE, "usage: %s <K> <T> <R> <file>\n", argv[0]);

    int K = strtol(argv[1], NULL, 10);
    if (K < 10 || K > 56403)
        errx(EXIT_FAILURE, "K [packet_count] out of range (10,56403)\n");

    int T = strtol(argv[2], NULL, 10);
    if (T < 8 || T > 65536)
        errx(EXIT_FAILURE, "T [packet_size] out of range (8,65536)\n");

    int R = strtol(argv[3], NULL, 10);
    if (R < 1)
        errx(EXIT_FAILURE, "R [repair_count] should be > 0\n");

    if (0 != nanorq_encoder_new(K, 0, &rq))
        errx(EXIT_FAILURE, "failed to init codec\n");

    calc_start = now();
    size_t prep_len = nanorq_calculate_prepare_memory(&rq);
    uint8_t *prep_mem = malloc(prep_len);
    nanorq_prepare(&rq, prep_mem, prep_len);

    size_t work_len = nanorq_calculate_work_memory(&rq);
    uint8_t *work_mem = malloc(work_len);
    nanorq_set_op_callback(&rq, &S, ops_push);
    ok = nanorq_precalculate(&rq, work_mem, work_len);
    assert(ok);
    calc_end = now();

    u32 rows = nanorq_get_pc_rows(&rq);
    u32 SH = nanorq_get_pc_genc_offset(&rq);
    u32 mem = rows * PAD(T);
    u8 *base = malloc(mem);
    u8_vec_init(&D, base, mem, mem, PAD(T));
    prepare_data_mat(&D, argv[4], rows, T, K, SH);

    ops_start = now();
    ops_run(&rq, &D, &S);
    ops_end = now();

    u8 reppkt[PAD(T)];
    for (u32 rp = 0; rp < R; rp++) {
        ops_mix(&rq, &D, K + rp, reppkt);
        fprintf(stdout, "RP %3d:", rp);
        for (u32 i = 0; i < T; i++)
            fprintf(stdout, "%02x", reppkt[i]);
        fprintf(stdout, "\n");
    }

    double calc_time = timespec_diff(&calc_end, &calc_start) / (double)NSECS;
    double ops_time = timespec_diff(&ops_end, &ops_start) / (double)NSECS;

    fprintf(stderr, "calc: %.6fs ops: %.6fs bytes: %u\n", calc_time, ops_time, T * K);

    free(prep_mem);
    free(work_mem);
    free(base);
    kv_destroy(S.ops);

    return 0;
}
