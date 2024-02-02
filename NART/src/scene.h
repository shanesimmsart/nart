#pragma once

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

#include "geometry.h"
#include "lights.h"
#include "bvh.h"

std::shared_ptr<TriMesh> LoadMeshFromFile(std::string filePath, glm::mat4& objectToWorld, std::shared_ptr<Material> material);

struct Camera
{
    Camera(float fov, glm::mat4 cameraToWorld);

    const float fov = 10.f;
    const glm::mat4 cameraToWorld = glm::mat4(1.f);
};

// Information about the render settings to be used
struct RenderInfo
{
    RenderInfo(uint32_t imageWidth, uint32_t imageHeight, uint32_t bucketSize, uint32_t spp, float filterWidth);

    const uint32_t imageWidth;
    const uint32_t imageHeight;
    const uint32_t bucketSize;
    const uint32_t spp;
    const float filterWidth;

    // Discrete bounds of filter
    uint32_t filterBounds;
    // Size of each tile, including border added to include filter width
    uint32_t tileSize;
    // Total image dimensions, including border added to include filter width
    uint32_t totalWidth;
    uint32_t totalHeight;
};

class Scene
{
public:
    Scene(RenderInfo info, Camera camera, std::vector<std::shared_ptr<Light>> lights, std::shared_ptr<BVH> bvh);

    bool Intersect(const Ray& ray, Intersection* isect) const;

    RenderInfo info;
    Camera camera;
    std::vector<std::shared_ptr<Light>> lights;
    std::shared_ptr<BVH> bvh;
};


Scene LoadScene(std::string scenePath);


