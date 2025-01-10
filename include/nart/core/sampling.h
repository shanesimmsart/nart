#pragma once

#include "rng.h"

#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <random>

glm::vec2 UniformSampleDisk(glm::vec2 sample);

glm::vec2 UniformSampleRing(glm::vec2 sample, float& pdf, float innerRadius);

glm::vec3 CosineSampleHemisphere(glm::vec2 sample, float& pdf);

float StratifiedSample1D(RNG& rng, uint32_t n, uint32_t nSamples);

void LatinSquare(RNG& rng, uint32_t nSamples, std::vector<glm::vec2>& samples);


