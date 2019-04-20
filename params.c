#include "params.h"
#include "rand.h"

struct ptuple {
  uint16_t d;
  uint16_t a;
  uint16_t b;
  uint16_t d1;
  uint16_t a1;
  uint16_t b1;
};

static const uint32_t degree_dist[] = {
    0,       5243,    529531,  704294,  791675,  844104,  879057,  904023,
    922747,  937311,  948962,  958494,  966438,  973160,  978921,  983914,
    988283,  992138,  995565,  998631,  1001391, 1003887, 1006157, 1008229,
    1010129, 1011876, 1013490, 1014983, 1016370, 1017662, 1048576};

static const uint16_t degree_dist_size =
    sizeof(degree_dist) / sizeof(degree_dist[0]);

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

static uint16_t deg(uint32_t v, uint16_t W) {
  for (int d = 0; d < degree_dist_size; d++) {
    if (v < degree_dist[d])
      return (d < (W - 2)) ? d : (W - 2);
  }
  return 0;
}

static struct ptuple gen_tuple(uint32_t X, uint16_t J, uint16_t W,
                               uint16_t P1) {

  struct ptuple ret;

  size_t A = 53591 + J * 997;

  if (A % 2 == 0)
    A++;
  size_t B1 = 10267 * (J + 1);
  uint32_t y = (uint32_t)((B1 + X * A));
  uint32_t v = rnd_get(y, 0, (uint32_t)((1 << 20)));
  ret.d = deg(v, W);
  ret.a = 1 + (uint16_t)(rnd_get(y, 1, W - 1));
  ret.b = (uint16_t)(rnd_get(y, 2, W));
  if (ret.d < 4) {
    ret.d1 = 2 + (uint16_t)(rnd_get(X, 3, 2));
  } else {
    ret.d1 = 2;
  }
  ret.a1 = 1 + (uint16_t)(rnd_get(X, 4, P1 - 1));
  ret.b1 = (uint16_t)(rnd_get(X, 5, P1));

  return ret;
}

struct pparams params_init(uint16_t symbols) {
  uint16_t idx;
  struct pparams prm = {0};

  for (idx = 0; idx < K_padded_size; idx++) {
    if (K_padded[idx] >= symbols) {
      prm.K_padded = K_padded[idx];
      break;
    }
  }

  prm.J = J_K_padded[idx];
  prm.S = S_H_W[idx][0];
  prm.H = S_H_W[idx][1];
  prm.W = S_H_W[idx][2];

  prm.L = prm.K_padded + prm.S + prm.H;
  prm.P = prm.L - prm.W;
  prm.U = prm.P - prm.H;
  prm.B = prm.W - prm.S;
  prm.P1 = prm.P + 1;

  while (!is_prime(prm.P1))
    prm.P1++;

  return prm;
}

uint16_vec params_get_idxs(struct pparams *prm, uint32_t X) {
  uint16_vec ret;
  struct ptuple t = gen_tuple(X, prm->J, prm->W, prm->P1);

  kv_init(ret);
  kv_resize(uint16_t, ret, t.d + t.d1);
  kv_push(uint16_t, ret, t.b);

  for (int j = 1; j < t.d; j++) {
    t.b = (t.b + t.a) % prm->W;
    kv_push(uint16_t, ret, t.b);
  }
  while (t.b1 >= prm->P)
    t.b1 = (t.b1 + t.a1) % prm->P1;

  kv_push(uint16_t, ret, prm->W + t.b1);
  for (int j = 1; j < t.d1; j++) {
    t.b1 = (t.b1 + t.a1) % prm->P1;
    while (t.b1 >= prm->P)
      t.b1 = (t.b1 + t.a1) % prm->P1;
    kv_push(uint16_t, ret, prm->W + t.b1);
  }
  return ret;
}
