#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include <precode.h>

void hdpc_matrix_print(params *P, pc *W, FILE *stream)
{
    u32 m = P->H, n = P->Kprime + P->S;
    fprintf(stream, "HDPC[%ux%u]\n", m, n);
    for (u32 row = 0; row < m; row++) {
        for (u32 col = 0; col < n; col++)
            fprintf(stream, "%02x", uv_E(W->HDPC, row, col));
        fprintf(stream, "\n");
    }
    fflush(stream);
}

int main(int argc, char *argv[])
{
    int ok = 0;
    void *sched = NULL;

    if (argc < 2) {
        fprintf(stderr, "usage: %s <K>\n", argv[0]);
        return -1;
    }

    int K = strtol(argv[1], NULL, 10);
    if (K < 5 || K > 56403) {
        fprintf(stderr, "failed to init codec\n");
        return -1;
    }

    params P = params_init(K);
    pc W = {};

    u32 mem = PAD(P.H * (P.Kprime + P.S));
    u32 tmp = P.H * (P.Kprime + P.S);
    u8 hdpc_base[mem];
    u8_vec_init(&W.HDPC, hdpc_base, tmp, tmp, P.Kprime + P.S);

    precode_matrix_make_HDPC(&P, &W);
    hdpc_matrix_print(&P, &W, stdout);
    return 0;
}
