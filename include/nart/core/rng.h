#pragma once

#include <cstdint>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <iostream>

class RNG {
public:
    void Seed(uint32_t seed) {
        // This is the seed Marsaglia uses
        y = seed + 2463534242;
    }

    float UniformFloat() {
        // Marsaglia's Xorshift RNG
        // This is equivalent to multiplying by a binary matrix T with a period
        // of 2^32-1 where T is created by multiplying a left shift matrix
        // several times, then a right shift matrix, and then another left one,
        // to shuffle the values Many different values have been found that
        // work, but these are Marsaglia's favourites
        y ^= (y << 13);
        y ^= (y >> 17);
        y ^= (y << 5);

        // When seeding from pixel coords that are spatially correlated, the
        // outputs are correlated due to linearity of the RNG Multiplying by
        // 0x9E3779BB introduces non-linearity, and mixes the bits up more
        // thoroughly 0x9E3779BB is fractional part of the golden ratio * 2^32,
        // and helps distribute values Small changes in the seed result in
        // significant, non-linear changes
        float f = glm::min(1.f - glm::epsilon<float>(),
                           float((y * 0x9E3779BB) * 2.3283064365386963e-10f));

        return f;
    }

    uint32_t UniformInt32(uint32_t max)  // 32-bit integer between 0 and max
    {
        y ^= (y << 13);
        y ^= (y >> 17);
        y ^= (y << 5);

        // Multiply-high method
        // This is essentially the same as taking y % (max + 1) to remap to [0,
        // max), but a bit more efficient and slightly less biased when (x + 1)
        // does not evenly divide 2^32 Treat y as a fraction, as though it's
        // divided by 2^32 Multiply by (max + 1) to remap to [0, max + 1) Take
        // the integer part by bit-shifting to the right by 32 and casting to 32
        // bits
        return (static_cast<uint64_t>(y * 0x9E3779B9) *
                (static_cast<uint64_t>(max) + 1)) >>
               32;
    }

private:
    uint32_t y = 2463534242;
};