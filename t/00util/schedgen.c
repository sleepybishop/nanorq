#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include <nanorq.h>

void sched_op(void *arg, u32 i, u16 j, u8 u)
{
    FILE *stream = (FILE *)arg;
    fprintf(stream, "%u %u %u\n", i, j, u);
}

int main(int argc, char *argv[])
{
    nanorq rq;

    if (argc < 2) {
        fprintf(stderr, "usage: %s <K>\n", argv[0]);
        return -1;
    }

    int K = strtol(argv[1], NULL, 10);
    if (K < 5 || K > 56403 || 0 != nanorq_encoder_new(K, 0, &rq)) {
        fprintf(stderr, "failed to init codec\n");
        return -1;
    }

    size_t prep_len = nanorq_calculate_prepare_memory(&rq);
    uint8_t *prep_mem = malloc(prep_len);
    nanorq_prepare(&rq, prep_mem, prep_len);

    size_t work_len = nanorq_calculate_work_memory(&rq);
    uint8_t *work_mem = malloc(work_len);
    nanorq_set_op_callback(&rq, stdout, sched_op);
    nanorq_precalculate(&rq, work_mem, work_len);
    fflush(stdout);
    free(prep_mem);
    free(work_mem);

    return 0;
}
