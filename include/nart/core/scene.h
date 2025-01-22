#pragma once

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

#include "../cameras/pinholecamera.h"
#include "../lights/disklight.h"
#include "../lights/distantlight.h"
#include "../lights/environmentlight.h"
#include "../lights/ringlight.h"
#include "../materials/diffusematerial.h"
#include "../materials/glassmaterial.h"
#include "../materials/glossydielectricmaterial.h"
#include "../materials/plasticmaterial.h"
#include "../materials/specularmaterial.h"
#include "../patterns/constantpattern.h"
#include "../patterns/texturepattern.h"
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
                                MaterialPtr&& material) const;

    void LoadMeshes(const nlohmann::json& json);

    void LoadCamera(const nlohmann::json& json);

    void LoadLights(const nlohmann::json& json);

    PatternPtr GetRho_d(const nlohmann::json& material);

    PatternPtr GetRho_s(const nlohmann::json& material);

    PatternPtr GetEta(const nlohmann::json& material);

    PatternPtr GetTau(const nlohmann::json& material);

    PatternPtr GetAlpha(const nlohmann::json& material);

    PatternPtr GetLe(const nlohmann::json& material);

    std::vector<LightPtr> lights;
    CameraPtr camera;
    BVHPtr bvh;
};
