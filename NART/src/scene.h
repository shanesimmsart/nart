#pragma once

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

#include "geometry.h"

std::shared_ptr<TriMesh> LoadMeshFromFile(std::string filePath, glm::mat4& objectToWorld, std::shared_ptr<Material> material)
{
    std::ifstream file;
    uint32_t numFaces;
    file.open(filePath);

    if (!file)
    {
        std::cerr << "Error: Mesh file " << filePath << " could not be opened.\n";
        return nullptr;
    }

    file >> numFaces;

    std::vector<uint32_t> faces(numFaces);
    uint32_t maxNormIndex = 0;

    uint32_t numVertIndices = 0;

    for (uint32_t i = 0; i < numFaces; ++i)
    {
        if (!(file >> faces[i]))
        {
            std::cerr << "Error: Mesh file could not be read.\n";
            return nullptr;
        }
        numVertIndices += faces[i];
    }

    std::vector<uint32_t> vertIndices(numVertIndices);
    uint32_t maxVertIndex = 0;

    uint32_t k = 0;

    for (uint32_t i = 0; i < numFaces; ++i)
    {
        for (uint32_t j = 0; j < faces[i]; ++j)
        {
            uint32_t vertIndex;

            if (!(file >> vertIndex))
            {
                std::cerr << "Error: Mesh file could not be read.\n";
                return nullptr;
            }

            vertIndices[k] = vertIndex;

            k += 1;

            maxVertIndex = glm::max(maxVertIndex, vertIndex);
        }
    }

    uint32_t numVertCoords = (maxVertIndex + 1) * 3;
    std::vector<float> vertCoords(numVertCoords);

    for (uint32_t i = 0; i < numVertCoords; ++i)
    {
        float vertCoord;
        if (!(file >> vertCoord))
        {
            std::cerr << "Error: Mesh file could not be read.\n";
            return nullptr;
        }

        vertCoords[i] = vertCoord;
    }

    std::vector<uint32_t> normIndices(numVertIndices);
    k = 0;

    for (uint32_t i = 0; i < numFaces; ++i)
    {
        for (uint32_t j = 0; j < faces[i]; ++j)
        {
            uint32_t normIndex;

            if (!(file >> normIndex))
            {
                std::cerr << "Error: Mesh file could not be read.\n";
                return nullptr;
            }

            normIndices[k] = normIndex;
            k += 1;

            maxNormIndex = glm::max(maxNormIndex, normIndex);
        }
    }

    uint32_t numNormCoords = (maxNormIndex + 1) * 3;
    std::vector<float> normCoords(numNormCoords);

    for (uint32_t i = 0; i < numNormCoords; ++i)
    {
        float normCoord;
        if (!(file >> normCoord))
        {
            std::cerr << "Error: Mesh file could not be read.\n";
            return nullptr;
        }

        normCoords[i] = normCoord;
    }



    // Calculate number of tris
    uint32_t numTris = 0;
    for (uint32_t i = 0; i < numFaces; ++i)
    {
        numTris += faces[i] - 2;
    }

    // Create vector of vertices
    std::vector<glm::vec3> verts(maxVertIndex + 1);
    uint32_t j = 0;
    for (uint32_t i = 0; i < maxVertIndex + 1; ++i)
    {
        verts[i] = glm::transpose(objectToWorld) * glm::vec4(vertCoords[j], vertCoords[j + 1], vertCoords[j + 2], 1.f);
        j += 3;
    }

    // Create vector of normals
    std::vector<glm::vec3> norms(maxNormIndex + 1);
    j = 0;
    for (uint32_t i = 0; i < maxNormIndex + 1; ++i)
    {
        // Need to transform normals by inverse transpose - otherwise the normal will get scaled by any non-uniform scaling
        norms[i] = glm::normalize(glm::vec3(glm::inverse(((objectToWorld))) * glm::vec4(normCoords[j], normCoords[j + 1], normCoords[j + 2], 0.f)));
        j += 3;
    }

    // Create tri indices
    std::vector<uint32_t> triIndices(numTris * 3);
    k = 0;
    uint32_t l = 0;
    for (uint32_t i = 0; i < numFaces; ++i)
    {
        for (uint32_t j = 0; j < faces[i] - 2; ++j)
        {
            triIndices[k] = vertIndices[l];
            triIndices[k + 1] = vertIndices[l + j + 1];
            triIndices[k + 2] = vertIndices[l + j + 2];
            k += 3;
        }
        l += faces[i];
    }

    // Create tri normal indices
    std::vector<uint32_t> triNormIndices(numTris * 3);
    k = 0;
    l = 0;
    for (uint32_t i = 0; i < numFaces; ++i)
    {
        for (uint32_t j = 0; j < faces[i] - 2; ++j)
        {
            triNormIndices[k] = normIndices[l];
            triNormIndices[k + 1] = normIndices[l + j + 1];
            triNormIndices[k + 2] = normIndices[l + j + 2];
            k += 3;
        }
        l += faces[i];
    }

    std::shared_ptr<TriMesh> mesh = std::make_unique<TriMesh>(material);
    mesh->triangles.reserve(numTris);

    // Now we can finally build our triangles
    j = 0;
    for (uint32_t i = 0; i < numTris; ++i)
    {
        Triangle tri(verts[triIndices[j]], verts[triIndices[j + 1]], verts[triIndices[j + 2]], norms[triNormIndices[j]], norms[triNormIndices[j + 1]], norms[triNormIndices[j + 2]]);
        mesh->triangles.emplace_back(tri);

        j += 3;
    }

    return mesh;
}

struct Camera
{
    Camera(float fov, glm::mat4 cameraToWorld) : fov(fov), cameraToWorld(cameraToWorld) {}

    const float fov = 10.f;
    const glm::mat4 cameraToWorld = glm::mat4(1.f);
};

// Information about the render settings to be used
struct RenderInfo
{
    RenderInfo(uint32_t imageWidth, uint32_t imageHeight, uint32_t bucketSize, uint32_t spp, float filterWidth) :
        imageWidth(imageWidth), imageHeight(imageHeight), bucketSize(bucketSize), spp(spp), filterWidth(filterWidth)
    {
        filterBounds = static_cast<uint32_t>(glm::ceil(filterWidth));
        tileSize = bucketSize + (filterBounds * 2);
        totalWidth = imageWidth + (filterBounds * 2);
        totalHeight = imageHeight + (filterBounds * 2);
    }

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
    Scene(RenderInfo info, Camera camera, std::vector<std::shared_ptr<Light>> lights, std::shared_ptr<BVH> bvh) : info(info), camera(camera), lights(lights), bvh(bvh) {}

    bool Intersect(const Ray& ray, Intersection* isect) const
    {
        return bvh->Intersect(ray, isect);
    }

    RenderInfo info;
    Camera camera;
    std::vector<std::shared_ptr<Light>> lights;
    std::shared_ptr<BVH> bvh;
};


Scene LoadScene(std::string scenePath)
{
    nlohmann::json json;
    std::ifstream ifs;
    ifs.open(scenePath);

    if (!ifs)
    {
        std::cerr << "Error: Scene file could not be opened.\n";
    }

    try
    {
        ifs >> json;
    }

    catch (nlohmann::json::exception& e)
    {
        std::cerr << "Error parsing scene: " << e.what() << "\nAborting.\n";
        std::abort();
    }

    nlohmann::json jsonInfo = json["renderInfo"];

    // Default values if not found in json
    uint32_t imageWidth = 64;
    uint32_t imageHeight = 64;
    uint32_t bucketSize = 32;
    uint32_t spp = 1;
    float filterWidth = 1.f;

    if (!jsonInfo.is_null())
    {
        try
        {
            imageWidth = jsonInfo["imageWidth"].get<uint32_t>();
            imageHeight = jsonInfo["imageHeight"].get<uint32_t>();
            bucketSize = jsonInfo["bucketSize"].get<uint32_t>();
            spp = jsonInfo["spp"].get<uint32_t>();
            filterWidth = jsonInfo["filterWidth"].get<float>();
        }

        catch (nlohmann::json::exception& e)
        {
            std::cerr << "Error in renderInfo: " << e.what() << "\n";
        }
    }

    else
    {
        std::cerr << "Warning: No render info found.\n";
    }

    RenderInfo info(imageWidth, imageHeight, bucketSize, spp, filterWidth);

    // Default values if not found in json
    glm::mat4 cameraToWorld(1.f);
    float fov = 1.f;
    std::vector<float> camTransform = { 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f };

    try {
        fov = json["camera"]["fov"].get<float>();
        camTransform = json["camera"]["transform"].get<std::vector<float>>();
    }

    catch (nlohmann::json::exception& e)
    {
        std::cerr << "Error in camera: " << e.what() << "\n";
    }

    uint8_t j = 0;
    for (uint8_t i = 0; i < 4; ++i)
    {
        cameraToWorld[i] = glm::vec4(camTransform[j], camTransform[j + 1], camTransform[j + 2], camTransform[j + 3]);
        j += 4;
    }

    Camera camera(fov, cameraToWorld);

    std::vector<std::string> mesh_filePaths;
    std::vector<std::shared_ptr<Material>> mesh_materials;
    if (!json["meshes"].is_null()) {
        for (auto& elem : json["meshes"])
        {
            std::string filePath = "";
            try
            {
                filePath = elem["filePath"].get<std::string>();
            }
            catch (nlohmann::json::exception& e)
            {
                std::cerr << "Error in mesh file path: " << e.what() << "\n";
            }
            mesh_filePaths.push_back(filePath);

            for (auto& mat : elem["material"])
            {
                std::shared_ptr<Material> material;
                try
                {
                    std::string type = mat["type"].get<std::string>();

                    if (type == "lambert")
                    {
                        std::vector<float> rhoGet = mat["rho"].get<std::vector<float>>();
                        glm::vec3 rho = glm::vec3(rhoGet[0], rhoGet[1], rhoGet[2]);
                        material = std::make_shared<DiffuseMaterial>(rho);
                    }

                    else if (type == "specular")
                    {
                        std::vector<float> RGet = mat["R"].get<std::vector<float>>();
                        glm::vec3 R = glm::vec3(glm::min(RGet[0], 1.f - glm::epsilon<float>()), glm::min(RGet[1], 1.f - glm::epsilon<float>()), glm::min(RGet[2], 1.f - glm::epsilon<float>()));
                        float eta = mat["eta"].get<float>();
                        material = std::make_shared<SpecularMaterial>(R, eta);
                    }

                    else if (type == "glossy")
                    {
                        std::vector<float> RGet = mat["R"].get<std::vector<float>>();
                        glm::vec3 R = glm::vec3(glm::min(RGet[0], 1.f - glm::epsilon<float>()), glm::min(RGet[1], 1.f - glm::epsilon<float>()), glm::min(RGet[2], 1.f - glm::epsilon<float>()));
                        float eta = mat["eta"].get<float>();
                        float alpha = mat["alpha"].get<float>();
                        material = std::make_shared<GlossyDielectricMaterial>(R, eta, alpha);
                    }

                    else if (type == "plastic")
                    {
                        std::vector<float> rhoGet = mat["rho"].get<std::vector<float>>();
                        glm::vec3 rho = glm::vec3(rhoGet[0], rhoGet[1], rhoGet[2]);
                        std::vector<float> RGet = mat["R"].get<std::vector<float>>();
                        glm::vec3 R = glm::vec3(glm::min(RGet[0], 1.f - glm::epsilon<float>()), glm::min(RGet[1], 1.f - glm::epsilon<float>()), glm::min(RGet[2], 1.f - glm::epsilon<float>()));
                        float eta = mat["eta"].get<float>();
                        float alpha = mat["alpha"].get<float>();
                        material = std::make_shared<PlasticMaterial>(rho, R, eta, alpha);
                    }

                    else
                    {
                        std::cerr << "Error: '" << type << "' is not a material type.\nAborting.\n";
                        std::abort();
                    }
                }
                catch (nlohmann::json::exception& e)
                {
                    std::cerr << "Error in mesh material type: " << e.what() << "\n";
                }
                mesh_materials.push_back(material);
            }
        }
    }

    else
    {
        std::cerr << "Error: No meshes found.\n";
    }

    std::vector<glm::mat4> meshTransforms;
    if (!json["meshes"].is_null()) {
        for (auto& elem : json["meshes"])
        {
            std::vector<float> meshTransform = { 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f };
            try
            {
                meshTransform = elem["transform"].get<std::vector<float>>();
            }
            catch (nlohmann::json::exception& e)
            {
                std::cerr << "Error in mesh transform: " << e.what() << "\n";
            }

            glm::mat4 objectToWorld(1.f);
            uint8_t j = 0;
            for (uint8_t i = 0; i < 4; ++i)
            {
                objectToWorld[i] = glm::vec4(meshTransform[j], meshTransform[j + 1], meshTransform[j + 2], meshTransform[j + 3]);
                j += 4;
            }
            meshTransforms.push_back(objectToWorld);
        }
    }

    std::vector<std::shared_ptr<TriMesh>> meshes;

    for (uint32_t i = 0; i < mesh_filePaths.size(); ++i) {
        meshes.push_back(LoadMeshFromFile(mesh_filePaths[i], meshTransforms[i], mesh_materials[i]));
    }

    std::cout << "Building BVH...\n";
    std::shared_ptr<BVH> bvh = std::make_shared<BVH>(BVH(meshes));

    std::vector<std::shared_ptr<Light>> lights;

    if (!json["lights"].is_null()) {
        for (auto& elem : json["lights"])
        {
            std::vector<float> lightTransform = { 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f };
            try
            {
                lightTransform = elem["transform"].get<std::vector<float>>();
            }
            catch (nlohmann::json::exception& e)
            {
                std::cerr << "Error in light transform: " << e.what() << "\n";
            }

            glm::mat4 lightToWorld(1.f);
            uint8_t j = 0;
            for (uint8_t i = 0; i < 4; ++i)
            {
                lightToWorld[i] = glm::vec4(lightTransform[j], lightTransform[j + 1], lightTransform[j + 2], lightTransform[j + 3]);
                j += 4;
            }

            try
            {
                std::string type = elem["type"].get<std::string>();

                if (type == "distant")
                {
                    std::vector<float> LeGet = elem["Le"].get<std::vector<float>>();
                    float intensity = elem["intensity"].get<float>();
                    glm::vec3 Le = glm::vec3(LeGet[0], LeGet[1], LeGet[2]);
                    std::shared_ptr<Light> light = std::make_shared<DistantLight>(Le, intensity, lightToWorld);
                    lights.push_back(light);
                }

                if (type == "disk")
                {
                    float radius = elem["radius"].get<float>();
                    std::vector<float> LeGet = elem["Le"].get<std::vector<float>>();
                    float intensity = elem["intensity"].get<float>();
                    glm::vec3 Le = glm::vec3(LeGet[0], LeGet[1], LeGet[2]);
                    std::shared_ptr<Light> light = std::make_shared<DiskLight>(radius, Le, intensity, lightToWorld);
                    lights.push_back(light);
                }

                if (type == "ring")
                {
                    float radius = elem["radius"].get<float>();
                    float innerRadius = elem["innerRadius"].get<float>();
                    std::vector<float> LeGet = elem["Le"].get<std::vector<float>>();
                    float intensity = elem["intensity"].get<float>();
                    glm::vec3 Le = glm::vec3(LeGet[0], LeGet[1], LeGet[2]);
                    std::shared_ptr<Light> light = std::make_shared<RingLight>(radius, innerRadius, Le, intensity, lightToWorld);
                    lights.push_back(light);
                }
            }

            catch (nlohmann::json::exception& e)
            {
                std::cerr << "Error in mesh material type: " << e.what() << "\n";
            }
        }
    }

    return Scene(info, camera, lights, bvh);
}


