#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <err.h>

#include <util.h>

#define NSECS 1000000000L

static u64 timespec_subtract(const struct timespec *s, const struct timespec *e)
{
    u64 sec = (u64)e->tv_sec - (u64)s->tv_sec;
    i64 nsec = (i64)e->tv_nsec - (i64)s->tv_nsec;
    return (u64)((u64)sec * NSECS) + nsec;
}

static inline struct timespec now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts;
}

int main(int argc, char *argv[])
{
    void *sched = NULL;

    if (argc < 2) {
        fprintf(stderr, "usage: %s <T>\n", argv[0]);
        return -1;
    }

    u32 T = strtol(argv[1], NULL, 10);
    if (T < 20 || T > 65536) {
        fprintf(stderr, "bad block size\n");
        return -1;
    }

    u32 mem = 256 * 1024 * 1024;
    u8 *base = calloc(1, mem);
    u8_vec data;
    u8_vec_init(&data, base, mem, mem, T);

    for (u32 x = 1; x < mem; x++)
        base[x] = rand() % 0xff;

    struct timespec axpy_start = now();
    for (u32 row = 1; row < mem / T; row++) {
        u8 u = rand() % 0xff;
        obl_axpy(uv_R(data, 0), uv_R(data, row), u, T);
    }
    struct timespec axpy_end = now();

    struct timespec xor_start = now();
    for (u32 row = 1; row < mem / T; row++)
        obl_axpy(uv_R(data, 0), uv_R(data, row), 1, T);
    struct timespec xor_end = now();

    fprintf(stdout, "AXPY/s: %.5f\n", 1.0 * NSECS * (mem / T) / timespec_subtract(&axpy_start, &axpy_end));
    fprintf(stdout, "1XPY/s: %.5f\n", 1.0 * NSECS * (mem / T) / timespec_subtract(&xor_start, &xor_end));

    return 0;
}
