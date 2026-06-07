#include "partition.h"
#include "sopi.h"
#include <assert.h>
#include <stdio.h>

void test_partition() {
  printf("Testing partition...\n");
  struct partition p = partition_fill(1000, 3);
  assert(p.IL == 334);
  assert(p.IS == 333);
  assert(p.JL == 1);
  assert(p.JS == 2);

  p = partition_fill(100, 10);
  assert(p.IL == 0);
  assert(p.IS == 10);
  assert(p.JL == 0);
  assert(p.JS == 10);

  printf("Partition tests passed.\n");
}

void test_sopi() {
  printf("Testing SOPI...\n");
  rq_sopi sopi = {100, 50, 200, 10};

  // Single ESI mapping test
  uint32_t esi = sopi_get_esi_single(&sopi, 5);
  uint64_t expected_esi = (100ULL + 5ULL * 50ULL) % SOPI_N;
  assert(esi == expected_esi);

  // Large index that overflows 32-bit but stays within limits
  uint64_t large_i = 100000000ULL;
  uint64_t large_b = 300ULL;
  uint32_t val = sopi_mersenne_mod(sopi.a, large_i, large_b);
  uint64_t expected_val = (sopi.a + large_i * large_b) % SOPI_N;
  assert(val == expected_val);

  // Multi block mapping test
  uint32_t Z = 5;
  uint32_t out_sbn, out_esi;
  sopi_get_mapping_multi(&sopi, 12, Z, &out_sbn, &out_esi);

  uint64_t r = 12 / Z;
  uint64_t expected_multi_esi = (sopi.a + r * sopi.b) % SOPI_N;
  uint32_t expected_multi_sbn = (12 + sopi.c + r * sopi.d) % Z;
  assert(out_esi == expected_multi_esi);
  assert(out_sbn == expected_multi_sbn);

  printf("SOPI tests passed.\n");
}

int main() {
  test_partition();
  test_sopi();
  printf("All utility tests passed successfully!\n");
  return 0;
}
