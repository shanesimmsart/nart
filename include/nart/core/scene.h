#pragma once

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

#include "../cameras/pinholecamera.h"
#include "../lights/disklight.h"
#include "../lights/distantlight.h"
#include "../lights/ringlight.h"
#include "../materials/diffusematerial.h"
#include "../materials/glassmaterial.h"
#include "../materials/glossydielectricmaterial.h"
#include "../materials/plasticmaterial.h"
#include "../materials/specularmaterial.h"
#include "bvh.h"
#include "geometry.h"

class Scene {
public:
    Scene(std::string scenePath);

    bool Intersect(const Ray& ray, Intersection& isect) const;

    const Light& GetLight(uint8_t index) const;

    uint8_t GetNumLights() const;

    const Camera& GetCamera() const;

    const BVH& GetBVH() const;

private:
    glm::mat4 MatrixFromVector(const std::vector<float>& vector) const;

    TriMeshPtr LoadMeshFromFile(const std::string& filePath,
                                glm::mat4& objectToWorld,
                                std::unique_ptr<Material>&& material) const;

    void LoadMeshes(const nlohmann::json& json);

    void LoadCamera(const nlohmann::json& json);

    void LoadLights(const nlohmann::json& json);

    std::vector<LightPtr> lights;
    CameraPtr camera;
    BVHPtr bvh;
};
