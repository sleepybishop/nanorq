#include "params.h"
#include "rand.h"
#include "tuple.h"

static bool is_prime(uint16_t n) {
  if (n <= 1)
    return false;
  if (n <= 3)
    return true;

  if (n % 2 == 0 || n % 3 == 0)
    return false;

  for (int i = 5; i * i <= n; i = i + 6)
    if (n % i == 0 || n % (i + 2) == 0)
      return false;

  return true;
}

params params_init(uint16_t K) {
  params P;

  for (unsigned i = 0; i < K_padded_size; i++) {
    if (K <= K_padded[i]) {
      P.Kprime = K_padded[i];
      P.J = J_K_padded[i];
      P.S = S_H_W[i][0];
      P.H = S_H_W[i][1];
      P.W = S_H_W[i][2];
      break;
    }
  }

  P.L = P.Kprime + P.S + P.H;
  P.P = P.L - P.W;
  P.U = P.P - P.H;
  P.B = P.W - P.S;
  P.P1 = P.P;

  while (!is_prime(P.P1))
    P.P1++;

  return P;
}

void params_set_idxs(uint32_t X, params *P, uint_vec *dst) {
  tuple t = gen_tuple(X, P);

  kv_push(unsigned, *dst, t.b);
  for (unsigned j = 1; j < t.d; j++) {
    t.b = (t.b + t.a) % P->W;
    kv_push(unsigned, *dst, t.b);
  }
  while (t.b1 >= P->P)
    t.b1 = (t.b1 + t.a1) % P->P1;

  kv_push(unsigned, *dst, P->W + t.b1);
  for (unsigned j = 1; j < t.d1; j++) {
    t.b1 = (t.b1 + t.a1) % P->P1;
    while (t.b1 >= P->P)
      t.b1 = (t.b1 + t.a1) % P->P1;
    kv_push(unsigned, *dst, P->W + t.b1);
  }
}
