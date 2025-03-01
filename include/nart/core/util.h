#pragma once

#include <cstdint>
#include <vector>

//  Note: This function assumes values between 0 and 1 - epsilon
uint32_t BinarySearch(float value, const std::vector<float>& v, uint32_t start,
                      uint32_t end);
