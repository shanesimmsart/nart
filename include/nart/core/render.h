#pragma once

#include <tbb/task_group.h>
#include <OpenEXR/ImfRgbaFile.h>
#include <OpenEXR/ImfArray.h>

#include "scene.h"
#include "geometry.h"
#include "light.h"

struct RenderParams
{
    // Default render parameters
    // (Overridden by scene file if left at default)
    uint32_t imageWidth = 0;
    uint32_t imageHeight = 0;
    uint32_t bucketSize = 0;
    uint32_t spp = 0;
    float filterWidth = -1.f;
    float rougheningFactor = -1.f;
};

struct Pixel
{
    glm::vec4 contribution = glm::vec4(0.f);
    float filterWeightSum = 0.f;
};

inline float Gaussian(float width, float x)
{
    if (x >= width) return 0.f;

    // In a Gaussian distribution, any value beyond 3 standard deviations is negligible (<0.3%)
    float sigma = width / 3.f;

    return (1.f / glm::sqrt(2.f * glm::pi<float>() * sigma * sigma)) * glm::exp(-(x * x) / (2.f * sigma * sigma));
}

class RenderSession
{
public:
    RenderSession(const Scene& scene, RenderParams params);

    void AddSample(const float* filterTable, const glm::vec2& sampleCoords, const glm::vec4& L, std::vector<Pixel>& pixels) const;

    std::vector<Pixel> RenderTile(const float* filterTable, uint32_t x0, uint32_t x1, uint32_t y0, uint32_t y1) const;

    glm::vec3 EstimateDirect(const glm::vec3 wo, BSDF& bsdf, const Intersection& isect, const Ray& ray, RNG& rng, uint8_t& flags) const;

    std::vector<Pixel> Render() const;

    void WriteImageToEXR(const std::vector<Pixel>& pixels, const char* filePath) const;

private:
    const Scene& scene;

    const RenderParams params;

    // Discrete bounds of filter
    uint32_t filterBounds;
    // Size of each tile, including filter width border
    uint32_t tileSize;
    // Total image dimensions, including filter width border
    uint32_t totalWidth;
    uint32_t totalHeight;
 
    const float shadowBias = 0.01f;
};

using RenderSessionPtr = std::unique_ptr<RenderSession>;

// Parses render parameter arguments, returns true if successful
bool ParseRenderParamArguments(int argc, char* argv[], RenderParams& params);

std::vector<std::unique_ptr<RenderSession>> LoadSessions(const std::string& scenePath, const Scene& scene, const RenderParams& _params);


