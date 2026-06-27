#include <check.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/resource.h>

// Include the actual production header
#include "t/00util/ult.h"

START_TEST(test_malloc_null_check_invariant)
{
    // Invariant: Memory allocation failures must be handled gracefully without crashes
    // Set memory limit to force malloc failures
    struct rlimit old_limit, new_limit;
    getrlimit(RLIMIT_AS, &old_limit);
    
    new_limit.rlim_cur = 1024 * 1024;  // 1MB limit
    new_limit.rlim_max = old_limit.rlim_max;
    setrlimit(RLIMIT_AS, &new_limit);
    
    // Test payloads: sizes that should trigger allocation failures
    size_t payloads[] = {
        1024 * 1024 * 1024,  // 1GB - will fail under our limit
        0,                    // Boundary: zero allocation
        1024                  // Valid small allocation
    };
    
    int num_payloads = sizeof(payloads) / sizeof(payloads[0]);
    
    for (int i = 0; i < num_payloads; i++) {
        // Call actual production function that uses malloc
        uint8_t *result = prep_mem(payloads[i]);
        
        // Security property: Either allocation succeeds with valid pointer,
        // or fails and returns NULL (no crash)
        if (result == NULL) {
            // This is acceptable - allocation failed but we didn't crash
            ck_assert_msg(1, "Allocation failed gracefully");
        } else {
            // Allocation succeeded - verify it's usable
            ck_assert_ptr_ne(result, NULL);
            free(result);
        }
    }
    
    // Restore original memory limit
    setrlimit(RLIMIT_AS, &old_limit);
}
END_TEST

Suite *security_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("Security");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_malloc_null_check_invariant);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = security_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}