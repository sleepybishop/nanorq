#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <err.h>

#include <nanorq.h>
#include "operations.h"

static size_t get_filesize(const char *f)
{
    size_t sz = 0;
    FILE *h = fopen(f, "r");
    if (!h)
        return 0;
    fseek(h, 0, SEEK_END);
    sz = ftell(h);
    fclose(h);
    return sz;
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

    int F = get_filesize(argv[4]);
    if (F < 0 || F > 56403 * T)
        errx(EXIT_FAILURE, "Size of data file out of range (%d, %d)\n", 0, 56403 * T);

    if (0 != nanorq_encoder_new(K, 0, &rq))
        errx(EXIT_FAILURE, "failed to init codec\n");

    size_t prep_len = nanorq_calculate_prepare_memory(&rq);
    uint8_t *prep_mem = calloc(1, prep_len);
    ok = nanorq_prepare(&rq, prep_mem, prep_len);
    assert(ok);

    u32 rows = nanorq_get_pc_rows(&rq);
    u32 SH = nanorq_get_pc_genc_offset(&rq);
    u32 mem = rows * PAD(T);
    u8 *base = malloc(mem);
    u8_vec_init(&D, base, mem, mem, PAD(T));
    prepare_data_mat(&D, argv[4], rows, T, K, SH);

    size_t work_len = nanorq_calculate_work_memory(&rq);
    uint8_t *work_mem = calloc(1, work_len);
    nanorq_set_op_callback(&rq, &S, ops_push);
    ok = nanorq_precalculate(&rq, work_mem, work_len);
    assert(ok);

    ops_run(&rq, &D, &S);

    u8 reppkt[PAD(T)];
    for (u32 rp = 0; rp < R; rp++) {
        ops_mix(&rq, &D, K + rp, reppkt);
        fprintf(stdout, "RP %3d:", rp);
        for (u32 i = 0; i < T; i++) {
            fprintf(stdout, "%02x", reppkt[i]);
        }
        fprintf(stdout, "\n");
    }

    free(prep_mem);
    free(work_mem);
    free(base);
    kv_destroy(S.ops);

    return 0;
}
