#pragma once

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

#include "../cameras/pinholecamera.h"
#include "../materials/diffusematerial.h"
#include "../materials/glossydielectricmaterial.h"
#include "../materials/plasticmaterial.h"
#include "../materials/specularmaterial.h"
#include "../materials/speculardielectricmaterial.h"
#include "geometry.h"
#include "../lights/distantlight.h"
#include "../lights/disklight.h"
#include "../lights/ringlight.h"
#include "bvh.h"

class Scene
{
public:
    Scene(const Camera* camera, std::vector<const Light*> lights, const BVH* bvh);

    ~Scene()
    {
        delete camera;

        for (const Light* light : lights)
        {
            delete light;
        }

        lights.clear();

        delete bvh;
    }

    bool Intersect(const Ray& ray, Intersection* isect) const;

    const Camera* camera;
    std::vector<const Light*> lights;
    const BVH* bvh;
};

Scene LoadScene(std::string scenePath);
glm::mat4 SceneMatrixFromVector(std::vector<float> vector);
std::shared_ptr<TriMesh> LoadMeshFromFile(std::string filePath, glm::mat4& objectToWorld, std::shared_ptr<Material> material);
std::vector<std::shared_ptr<TriMesh>> LoadMeshes(const nlohmann::json json);


