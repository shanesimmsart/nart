#pragma once

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

#include "../cameras/pinholecamera.h"
#include "../materials/diffusematerial.h"
#include "../materials/glossydielectricmaterial.h"
#include "../materials/plasticmaterial.h"
#include "../materials/specularmaterial.h"
#include "../materials/glassmaterial.h"
#include "geometry.h"
#include "../lights/distantlight.h"
#include "../lights/disklight.h"
#include "../lights/ringlight.h"
#include "bvh.h"

class Scene
{
public:
    Scene(std::string scenePath);

    bool Intersect(const Ray& ray, Intersection* isect) const;

    const Light& GetLight(uint8_t index) const;

    uint8_t GetNumLights() const;

    const Camera& GetCamera() const;

    const BVH& GetBVH() const;

private:
    glm::mat4 MatrixFromVector(std::vector<float> vector);

    TriMeshPtr LoadMeshFromFile(std::string filePath, glm::mat4& objectToWorld, std::shared_ptr<Material> material);

    void LoadMeshes(const nlohmann::json json);
    
    void LoadCamera(const nlohmann::json json);

    void LoadLights(const nlohmann::json json);

    std::vector<LightPtr> lights;
    CameraPtr camera;
    BVHPtr bvh;
};


