#include "tuple.h"
#include "rand.h"

static const uint32_t degree_dist[] = {
    0,       5243,    529531,  704294,  791675,  844104,  879057,  904023,
    922747,  937311,  948962,  958494,  966438,  973160,  978921,  983914,
    988283,  992138,  995565,  998631,  1001391, 1003887, 1006157, 1008229,
    1010129, 1011876, 1013490, 1014983, 1016370, 1017662, 1048576};

static const uint16_t degree_dist_size =
    sizeof(degree_dist) / sizeof(degree_dist[0]);

static uint16_t deg(uint32_t v, uint16_t W) {
  for (int d = 0; d < degree_dist_size; d++) {
    if (v < degree_dist[d])
      return (d < (W - 2)) ? d : (W - 2);
  }
  return 0;
}

tuple gen_tuple(uint32_t X, params *P) {
  tuple ret;

  size_t A = 53591 + P->J * 997;

  if (A % 2 == 0)
    A++;
  size_t B1 = 10267 * (P->J + 1);
  uint32_t y = (uint32_t)((B1 + X * A));
  uint32_t v = rnd_get(y, 0, (uint32_t)((1 << 20)));
  ret.d = deg(v, P->W);
  ret.a = 1 + rnd_get(y, 1, P->W - 1);
  ret.b = rnd_get(y, 2, P->W);
  if (ret.d < 4) {
    ret.d1 = 2 + rnd_get(X, 3, 2);
  } else {
    ret.d1 = 2;
  }
  ret.a1 = 1 + rnd_get(X, 4, P->P1 - 1);
  ret.b1 = rnd_get(X, 5, P->P1);

  return ret;
}
