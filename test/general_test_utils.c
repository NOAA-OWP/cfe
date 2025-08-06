#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "general_test_utils.h"

/**
 * Helper func to dynamically allocate memory for an array of strings.
 *
 * @param str_count The number of strings
 * @param max_str_length The maximum size of each string
 * @return Pointer to the allocated memory
 */
char** allocate_array_of_strings(int str_count, int max_str_length)
{
    char** string_array = malloc(sizeof(char *) * str_count);
    for (int i = 0; i < str_count; i++)
        string_array[i] = (char*) malloc(sizeof(char) * max_str_length);
    return string_array;
}

/**
 * Free memory from an array of strings dynamically allocated using ``allocate_array_of_strings``.
 *
 * @param array_of_strings The (outer) array to free
 * @param str_count The number of strings (i.e., inner arrays) free individually
 */
void free_array_of_strings(char** array_of_strings, int str_count)
{
    for (int i = 0; i < str_count; i++)
        free(array_of_strings[i]);
    free(array_of_strings);
}

/**
 * Reusable convenience method for checking that ``double`` results of a BMI function call were as expected.
 *
 * @param expected The expected value.
 * @param actual The actual value set by the function call.
 * @return
 */
bool confirm_matches_expected_doubles(const double expected, const double actual)
{
    if (expected != actual) {
        printf("Expected:  %f  |  Value:  %f\n", expected, actual);
        return false;
    }
    return true;
}

/**
 * Confirm the actual matches the expected for arrays of strings, ignoring the order within the outer array.
 *
 * @param expected
 * @param actual
 * @param array_size
 * @return
 */
bool confirm_matches_expected_str_arrays(char** expected, char** actual, const int array_size)
{
    for (int i = 0; i < array_size; i++) {
        int match_found = find_in_array(expected[i], actual, array_size);
        if (match_found < 0) {
            printf("\nFailed to find expected string %s in actual resutls", expected[i]);
            return false;
        }
    }
    return true;
}

/**
 * Reusable convenience method for checking that ``char*`` results of a BMI function call were as expected.
 *
 * @param expected The expected value.
 * @param actual The actual value set by the function call.
 * @return
 */
bool confirm_matches_expected_strs(const char* expected, const char* actual)
{
    if (strcmp(expected, actual) != 0) {
        printf("\nExpected:  %s  |  Actual:  %s", expected, actual);
        return false;
    }
    return true;
}

/**
 * Reusable convenience method for checking that ``int`` results of a BMI function call were as expected.
 *
 * @param expected The expected value.
 * @param actual The actual value set by the function call.
 * @return
 */
bool confirm_matches_expected_ints(const int expected, const int actual)
{
    if (expected != actual) {
        printf("Expected:  %i  |  Value:  %i\n", expected, actual);
        return false;
    }
    return true;
}

/**
 * Convenience implementation to check if an int matches an ``expected`` value provided as a string.
 *
 * @param expected Expected int value, but represented as a string
 * @param actual Actual int value returned as part of a test.
 * @return
 */
bool confirm_matches_expected_ints_parse_str(const char* expected, const int actual)
{
    return confirm_matches_expected_ints(atoi(expected), actual);
}

/**
 * Search the array to see if it contains the given string.
 *
 * @param str The string to search for.
 * @param array_of_strings  The array to search.
 * @param array_size The size of the array to search.
 * @return The index with the array to search where the sought string can be found, or -1 if it is not found.
 */
int find_in_array(const char* str, char** array_of_strings, const int array_size)
{
    for (int i = 0; i < array_size; i++) {
        const char* val = array_of_strings[i];

        if (strcmp(val, str) == 0)
            return i;
    }
    return -1;
}