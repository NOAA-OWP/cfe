#ifndef GENERAL_TEST_UTILS
#define GENERAL_TEST_UTILS

#include <stdbool.h>

# define TEST_RETURN_CODE_PASS 0
# define TEST_RETURN_CODE_FAIL 1

/**
 * Helper func to dynamically allocate memory for an array of strings.
 *
 * @param str_count The number of strings
 * @param max_str_length The maximum size of each string
 * @return Pointer to the allocated memory
 */
char** allocate_array_of_strings(int str_count, int max_str_length);

/**
 * Free memory from an array of strings dynamically allocated using ``allocate_array_of_strings``.
 *
 * @param array_of_strings The (outer) array to free
 * @param str_count The number of strings (i.e., inner arrays) free individually
 */
void free_array_of_strings(char** array_of_strings, int str_count);

/**
 * Confirm the actual matches the expected for arrays of strings, ignoring the order within the outer array.
 *
 * @param expected
 * @param actual
 * @param array_size
 * @return
 */
bool confirm_matches_expected_str_arrays(const char** expected, const char** actual, const int array_size);

/**
 * Reusable convenience method for checking that ``char*`` results of a BMI function call were as expected.
 *
 * @param expected The expected value.
 * @param actual The actual value set by the function call.
 * @return
 */
bool confirm_matches_expected_strs(const char* expected, const char* actual);

/**
 * Reusable convenience method for checking that ``int`` results of a BMI function call were as expected.
 *
 * @param expected The expected value.
 * @param actual The actual value set by the function call.
 * @return
 */
bool confirm_matches_expected_ints(const int expected, const int actual);

/**
 * Convenience implementation to check if an int matches an ``expected`` value provided as a string.
 *
 * @param expected Expected int value, but represented as a string
 * @param actual Actual int value returned as part of a test.
 * @return
 */
bool confirm_matches_expected_ints_parse_str(const char* expected, const int actual);

/**
 * Search the array to see if it contains the given string.
 *
 * @param str The string to search for.
 * @param array_of_strings  The array to search.
 * @param array_size The size of the array to search.
 * @return The index with the array to search where the sought string can be found, or -1 if it is not found.
 */
int find_in_array(const char* str, const char** array_of_strings, const int array_size);

#endif