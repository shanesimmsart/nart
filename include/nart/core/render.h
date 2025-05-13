#pragma once

// TODO: GET LOOKUP TABLE WORKING!!! >:d
//       OR F16C SSE
#define IMATH_HALF_USE_LOOKUP_TABLE 0
#define IMATH_HALF_NO_LOOKUP_TABLE
#include <OpenEXR/ImfArray.h>
#include <OpenEXR/ImfRgbaFile.h>
#include <tbb/task_group.h>

#include "geometry.h"
#include "light.h"
#include "scene.h"

struct RenderParams {
    // Default render parameters
    // (Overridden by scene file if left at default)
    uint32_t imageWidth = 0;
    uint32_t imageHeight = 0;
    uint32_t bucketSize = 0;
    uint32_t spp = 0;
    uint32_t bounces = 0;
    float filterWidth = -1.f;
    float rougheningFactor = -1.f;
};

struct Pixel {
    glm::vec4 contribution = glm::vec4(0.f);
    float filterWeightSum = 0.f;
};

inline float Gaussian(float width, float x) {
    if (x >= width) return 0.f;

    // In a Gaussian distribution, any value beyond 3 standard deviations is
    // negligible (<0.3%)
    float sigma = width / 3.f;

    return (1.f / glm::sqrt(2.f * glm::pi<float>() * sigma * sigma)) *
           glm::exp(-(x * x) / (2.f * sigma * sigma));
}

struct IntersectionInfo {
    IntersectionInfo(uint32_t meshID, uint8_t priority, float eta)
        : meshID(meshID), priority(priority), eta(eta) {}

    uint32_t meshID = -1;
    uint8_t priority = 0;
    float eta = 1.f;
};

class RenderSession {
public:
    RenderSession(const Scene& scene, RenderParams params);

    void AddSample(const float* filterTable, const glm::vec2& sampleCoords,
                   const glm::vec4& L, std::vector<Pixel>& pixels) const;

    std::vector<Pixel> RenderTile(const float* filterTable, uint32_t x0,
                                  uint32_t x1, uint32_t y0, uint32_t y1) const;

    using IsectInfoList =
        std::vector<IntersectionInfo, ArenaAllocator<IntersectionInfo>>;

    // Check if we should ignore the intersection and update IOR of current
    // medium
    bool IsectIsValid(const Intersection& isect, const IsectInfoList& isectList,
                      float& eta_outer) const;

    glm::vec3 EstimateDirect(const glm::vec3 wo, BSDF& bsdf,
                             const Intersection& isect, const Ray& ray,
                             RNG& rng, uint8_t& flags, float eta_outer) const;

    void UpdateIsectList(IsectInfoList& isectList, const Intersection& isect,
                         float eta_sampled) const {
        // Add intersection info to list
        bool inList = false;
        uint32_t matchingIDIndex = 0;

        // Check if intersected mesh already in isectList
        for (uint32_t k = isectList.size(); k-- > 0;) {
            if (isectList[k].meshID == isect.meshID) {
                inList = true;
                matchingIDIndex = k;
                break;
            }
        }

        // If it isn't, we are entering said mesh, and add
        // it to the list
        if (!inList) {
            isectList.emplace_back(
                IntersectionInfo(isect.meshID, isect.priority, eta_sampled));
        }

        // If it is, we are exiting it, and it must be
        // removed
        else {
            isectList.erase(isectList.begin() + matchingIDIndex);
        }
    }

    std::vector<Pixel> Render() const;

    void WriteImageToEXR(const std::vector<Pixel>& pixels,
                         const char* filePath) const;

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

    const float shadowBias = 0.001f;
};

using RenderSessionPtr = std::unique_ptr<RenderSession>;

// Parses render parameter arguments, returns true if successful
bool ParseRenderParamArguments(int argc, char* argv[], RenderParams& params);

std::vector<std::unique_ptr<RenderSession>> LoadSessions(
    const std::string& scenePath, const Scene& scene,
    const RenderParams& _params);
