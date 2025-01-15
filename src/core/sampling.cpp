#include "../../include/nart/core/sampling.h"

#include "../../include/nart/core/rng.h"

glm::vec2 UniformSampleDisk(glm::vec2 sample) {
    float r = glm::sqrt(sample.x);
    float theta = sample.y * glm::two_pi<float>();

    float cosTheta = glm::cos(theta);
    float sinTheta = glm::sin(theta);

    float x = r * cosTheta;
    float y = r * sinTheta;

    return glm::vec2(x, y);
}

glm::vec2 UniformSampleRing(glm::vec2 sample, float& pdf, float innerRadius) {
    float r = glm::sqrt(glm::mix(innerRadius, 1.f, sample.x));
    float theta = sample.y * glm::two_pi<float>();

    float cosTheta = glm::cos(theta);
    float sinTheta = glm::sin(theta);

    float x = r * cosTheta;
    float y = r * sinTheta;

    pdf = 1.f / (glm::pi<float>() * (1.f - innerRadius));

    return glm::vec2(x, y);
}

glm::vec3 CosineSampleHemisphere(glm::vec2 sample, float& pdf) {
    // Using Malley's method, we can uniformly sample a disk, then project
    // the sample points upwards onto the hemisphere to achieve a
    // cosine-weighted distribution
    glm::vec2 diskSample = UniformSampleDisk(sample);
    float z = glm::sqrt(
        1.f - (diskSample.x * diskSample.x + diskSample.y * diskSample.y));

    pdf = z * glm::one_over_pi<float>();

    return glm::vec3(diskSample, z);
}

float StratifiedSample1D(RNG& rng, uint32_t n, uint32_t nSamples) {
    float invNSamples = 1.f / static_cast<float>(nSamples);
    return (static_cast<float>(n) + rng.UniformFloat()) * invNSamples;
}

// Using this technique instead of stratifying in x and y independently lets
// have any number of spp including primes, e.g 2, 3, 5, 7... and have them be
// evenly distributed
void LatinSquare(RNG& rng, uint32_t nSamples, std::vector<glm::vec2>& samples) {
    // Create stratified samples along diagonal
    for (uint32_t i = 0; i < nSamples; ++i) {
        samples[i] = glm::vec2(StratifiedSample1D(rng, i, nSamples),
                               StratifiedSample1D(rng, i, nSamples));
    }

    // Shuffle dimensions
    for (uint32_t i = 0; i < nSamples; ++i) {
        uint32_t choice = rng.UniformInt32(nSamples - 1 - i);
        std::swap(samples[i].x, samples[choice].x);
        choice = rng.UniformInt32(nSamples - 1 - i);
        std::swap(samples[i].y, samples[choice].y);
    }
}
