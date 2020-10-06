#include "params.h"
#include "rand.h"
#include "tuple.h"

static bool is_prime(uint16_t n) {
  if (n <= 3)
    return true;
  if (n % 2 == 0 || n % 3 == 0)
    return false;

  uint32_t i = 5;
  uint32_t w = 2;
  while (i * i <= n) {
    if (n % i == 0)
      return false;
    i += w;
    w = 6 - w;
  }
  return true;
}

params params_init(uint16_t K) {
  params P;

  for (int i = 0; i < K_padded_size; i++) {
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

int_vec params_get_idxs(uint32_t X, params *P) {
  int_vec ret;
  tuple t = gen_tuple(X, P);

  kv_init(ret);
  kv_resize(int, ret, t.d + t.d1);
  kv_push(int, ret, t.b);

  for (int j = 1; j < t.d; j++) {
    t.b = (t.b + t.a) % P->W;
    kv_push(int, ret, t.b);
  }
  while (t.b1 >= P->P)
    t.b1 = (t.b1 + t.a1) % P->P1;

  kv_push(int, ret, P->W + t.b1);
  for (int j = 1; j < t.d1; j++) {
    t.b1 = (t.b1 + t.a1) % P->P1;
    while (t.b1 >= P->P)
      t.b1 = (t.b1 + t.a1) % P->P1;
    kv_push(int, ret, P->W + t.b1);
  }
  return ret;
}
