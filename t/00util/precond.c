#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include <nanorq.h>

#define hf_at(W, i, j) (uv_A((W)->F.type, (i)) ? uv_E((W)->UL, uv_A((W)->F.rowmap, (i)), (j)) : bm_at(&(W)->U, (i), (j)))
void work_matrix_print(params *P, pc *W, FILE *stream)
{
    u32 m = W->rows, n = W->u;
    u32_vec *d = &W->d;
    fprintf(stream, "U[%ux%u]\n", m, n);
    for (u32 row = 0; row < m; row++) {
        int drow = d ? uv_A(*d, row) : row;
        for (u32 col = 0; col < n; col++) {
            u8 val = hf_at(W, drow, col);
            if (val > 1)
                fprintf(stream, "x");
            else
                fprintf(stream, "%d", val);
        }
        fprintf(stream, "\n");
    }
    fflush(stream);
}

int main(int argc, char *argv[])
{
    nanorq rq;
    void *sched = NULL;

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
    uint8_t *prep_mem = calloc(1, prep_len);
    nanorq_prepare(&rq, prep_mem, prep_len);
    size_t work_len = nanorq_calculate_work_memory(&rq);
    uint8_t *work_mem = calloc(1, work_len);
    nanorq_precalculate(&rq, work_mem, work_len);
    work_matrix_print(&rq.P, &rq.W, stdout);

    free(prep_mem);
    free(work_mem);

    return 0;
}
