#pragma once

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

#include "../cameras/pinholecamera.h"
#include "../materials/diffusematerial.h"
#include "../materials/glossydielectricmaterial.h"
#include "../materials/plasticmaterial.h"
#include "../materials/specularmaterial.h"
#include "geometry.h"
#include "../lights/distantlight.h"
#include "../lights/disklight.h"
#include "../lights/ringlight.h"
#include "bvh.h"

std::shared_ptr<TriMesh> LoadMeshFromFile(std::string filePath, glm::mat4& objectToWorld, std::shared_ptr<Material> material);

class Scene
{
public:
    Scene(std::shared_ptr<Camera> camera, std::vector<const Light*> lights, std::shared_ptr<BVH> bvh);

    ~Scene()
    {
        for (const Light* light : lights)
        {
            delete light;
        }

        lights.clear();
    }

    bool Intersect(const Ray& ray, Intersection* isect) const;

    std::shared_ptr<Camera> camera;
    std::vector<const Light*> lights;
    std::shared_ptr<BVH> bvh;
};


Scene LoadScene(std::string scenePath);


