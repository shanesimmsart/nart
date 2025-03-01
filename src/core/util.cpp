#include "../../include/nart/core/util.h"

//  Note: This function assumes values between 0 and 1 - epsilon
uint32_t BinarySearch(float value, const std::vector<float>& v, uint32_t start,
    uint32_t end) {
    uint32_t i = start;

    while (start < end) {
        i = start + ((end - start) / 2);

        if (v[i] > value) {
            end = i;
            i -= 1;
        }

        else {
            start = i + 1;
        }
    }

    return i;
}


