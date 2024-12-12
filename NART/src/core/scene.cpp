#include "scene.h"

Scene::Scene(std::string scenePath)
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

    LoadCamera(json);

    LoadMeshes(json);

    LoadLights(json);
}

bool Scene::Intersect(const Ray& ray, Intersection* isect) const
{
    return bvh->Intersect(ray, isect);
}

const Light& Scene::GetLight(uint8_t index) const
{
    if (!(lights[index]))
    {
        throw std::runtime_error("Light is null");
    }

    return *(lights[index]);
}

uint8_t Scene::GetNumLights() const
{
    return lights.size();
}

const Camera& Scene::GetCamera() const
{
    if (!camera)
    {
        throw std::runtime_error("Camera is null");
    }

    return *(camera);
}

const BVH& Scene::GetBVH() const
{
    if (!bvh)
    {
        throw std::runtime_error("BVH is null");
    }

    return *(bvh);
}

glm::mat4 Scene::MatrixFromVector(std::vector<float> vector)
{
    glm::mat4 matrix(1.f);

    uint8_t j = 0;
    for (uint8_t i = 0; i < 4; ++i)
    {
        matrix[i] = glm::vec4(vector[j], vector[j + 1], vector[j + 2], vector[j + 3]);
        j += 4;
    }

    return matrix;
}

TriMeshPtr Scene::LoadMeshFromFile(std::string filePath, glm::mat4& objectToWorld, std::shared_ptr<Material> material)
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

    TriMeshPtr mesh = std::make_unique<TriMesh>(material);
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

void Scene::LoadMeshes(const nlohmann::json json)
{
    std::vector<std::string> meshFilePaths;
    std::vector<std::shared_ptr<Material>> meshMaterials;
    std::vector<glm::mat4> meshTransforms;
    if (!json["meshes"].is_null()) {
        uint8_t numMeshes = json["meshes"].size();

        meshFilePaths.reserve(numMeshes);
        meshMaterials.reserve(numMeshes);
        meshTransforms.reserve(numMeshes);

        for (auto& elem : json["meshes"])
        {
            std::string filePath = "";
            try
            {
                meshFilePaths.push_back(elem["filePath"].get<std::string>());
            }
            catch (nlohmann::json::exception& e)
            {
                std::cerr << "Error in mesh file path: " << e.what() << "\n";
            }

            for (auto& mat : elem["material"])
            {
                std::shared_ptr<Material> material;
                try
                {
                    std::string type = mat["type"].get<std::string>();

                    if (type == "lambert")
                    {
                        std::vector<float> rho_dGet = mat["rho_d"].get<std::vector<float>>();
                        glm::vec3 rho_d = glm::vec3(rho_dGet[0], rho_dGet[1], rho_dGet[2]);
                        meshMaterials.push_back(std::make_shared<DiffuseMaterial>(rho_d));
                    }

                    else if (type == "specular")
                    {
                        std::vector<float> rho_sGet = mat["rho_s"].get<std::vector<float>>();
                        glm::vec3 rho_s = glm::vec3(glm::min(rho_sGet[0], 1.f - glm::epsilon<float>()), glm::min(rho_sGet[1], 1.f - glm::epsilon<float>()), glm::min(rho_sGet[2], 1.f - glm::epsilon<float>()));
                        float eta = mat["eta"].get<float>();
                        meshMaterials.push_back(std::make_shared<SpecularMaterial>(rho_s, eta));
                    }

                    else if (type == "glass")
                    {
                        std::vector<float> rho_sGet = mat["rho_s"].get<std::vector<float>>();
                        glm::vec3 rho_s = glm::vec3(glm::min(rho_sGet[0], 1.f - glm::epsilon<float>()), glm::min(rho_sGet[1], 1.f - glm::epsilon<float>()), glm::min(rho_sGet[2], 1.f - glm::epsilon<float>()));
                        std::vector<float> tauGet = mat["tau"].get<std::vector<float>>();
                        glm::vec3 tau = glm::vec3(glm::min(tauGet[0], 1.f - glm::epsilon<float>()), glm::min(tauGet[1], 1.f - glm::epsilon<float>()), glm::min(tauGet[2], 1.f - glm::epsilon<float>()));
                        float eta = mat["eta"].get<float>();
                        float roughness = mat["roughness"].get<float>();
                        float alpha = roughness * roughness;
                        meshMaterials.push_back(std::make_shared<GlassMaterial>(rho_s, tau, eta, alpha));
                    }

                    else if (type == "glossy")
                    {
                        std::vector<float> rho_sGet = mat["rho_s"].get<std::vector<float>>();
                        glm::vec3 rho_s = glm::vec3(glm::min(rho_sGet[0], 1.f - glm::epsilon<float>()), glm::min(rho_sGet[1], 1.f - glm::epsilon<float>()), glm::min(rho_sGet[2], 1.f - glm::epsilon<float>()));
                        float eta = mat["eta"].get<float>();
                        float roughness = mat["roughness"].get<float>();
                        float alpha = roughness * roughness;
                        meshMaterials.push_back(std::make_shared<GlossyDielectricMaterial>(rho_s, eta, alpha));
                    }

                    else if (type == "plastic")
                    {
                        std::vector<float> rho_dGet = mat["rho_d"].get<std::vector<float>>();
                        glm::vec3 rho_d = glm::vec3(rho_dGet[0], rho_dGet[1], rho_dGet[2]);
                        std::vector<float> rho_sGet = mat["rho_s"].get<std::vector<float>>();
                        glm::vec3 rho_s = glm::vec3(glm::min(rho_sGet[0], 1.f - glm::epsilon<float>()), glm::min(rho_sGet[1], 1.f - glm::epsilon<float>()), glm::min(rho_sGet[2], 1.f - glm::epsilon<float>()));
                        float eta = mat["eta"].get<float>();
                        float roughness = mat["roughness"].get<float>();
                        float alpha = roughness * roughness;
                        meshMaterials.push_back(std::make_shared<PlasticMaterial>(rho_d, rho_s, eta, alpha));
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
            }
        }

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

            glm::mat4 objectToWorld = MatrixFromVector(meshTransform);
            meshTransforms.push_back(objectToWorld);
        }
    }

    else
    {
        std::cerr << "Error: No meshes found.\n";
    }

    std::vector<TriMeshPtr> meshes;
    meshes.reserve(meshFilePaths.size());

    for (uint32_t i = 0; i < meshFilePaths.size(); ++i) {
        meshes.push_back(LoadMeshFromFile(meshFilePaths[i], meshTransforms[i], meshMaterials[i]));
    }

    bvh = std::make_unique<BVH>(std::move(meshes));
}

void Scene::LoadCamera(const nlohmann::json json)
{
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

    cameraToWorld = MatrixFromVector(camTransform);

    camera = std::make_unique<PinholeCamera>(fov, cameraToWorld);
}

void Scene::LoadLights(const nlohmann::json json)
{
    if (!json["lights"].is_null()) {

        uint8_t numLights = json["lights"].size();
        lights.reserve(numLights);

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

            glm::mat4 lightToWorld = MatrixFromVector(lightTransform);

            try
            {
                std::string type = elem["type"].get<std::string>();

                if (type == "distant")
                {
                    std::vector<float> LeGet = elem["Le"].get<std::vector<float>>();
                    float intensity = elem["intensity"].get<float>();
                    glm::vec3 Le = glm::vec3(LeGet[0], LeGet[1], LeGet[2]);
                    lights.emplace_back(std::make_unique<DistantLight>(Le, intensity, lightToWorld));
                }

                if (type == "disk")
                {
                    float radius = elem["radius"].get<float>();
                    std::vector<float> LeGet = elem["Le"].get<std::vector<float>>();
                    float intensity = elem["intensity"].get<float>();
                    glm::vec3 Le = glm::vec3(LeGet[0], LeGet[1], LeGet[2]);
                    lights.emplace_back(std::make_unique<DiskLight>(radius, Le, intensity, lightToWorld));
                }

                if (type == "ring")
                {
                    float radius = elem["radius"].get<float>();
                    float innerRadius = elem["innerRadius"].get<float>();
                    std::vector<float> LeGet = elem["Le"].get<std::vector<float>>();
                    float intensity = elem["intensity"].get<float>();
                    glm::vec3 Le = glm::vec3(LeGet[0], LeGet[1], LeGet[2]);
                    lights.emplace_back(std::make_unique<RingLight>(radius, innerRadius, Le, intensity, lightToWorld));
                }
            }

            catch (nlohmann::json::exception& e)
            {
                std::cerr << "Error in mesh material type: " << e.what() << "\n";
            }
        }
    }
}


