#include "../../include/nart/core/scene.h"

Scene::Scene(std::string scenePath) {
    nlohmann::json json;
    std::ifstream ifs;
    ifs.open(scenePath);

    if (!ifs) {
        std::cerr << "Error: Scene file could not be opened.\n";
    }

    try {
        ifs >> json;
    }

    catch (nlohmann::json::exception& e) {
        std::cerr << "Error parsing scene: " << e.what() << "\nAborting.\n";
        std::abort();
    }

    LoadCamera(json);

    LoadMeshes(json);

    LoadLights(json);
}

bool Scene::Intersect(const Ray& ray, Intersection& isect) const {
    return bvh->Intersect(ray, isect);
}

const Light& Scene::GetLight(uint8_t index) const {
    if (index >= lights.size()) {
        throw std::out_of_range("Light index is out of range");
    }

    if (!(lights[index])) {
        throw std::runtime_error("Light is null");
    }

    return *lights[index];
}

uint8_t Scene::GetNumLights() const { return lights.size(); }

const Camera& Scene::GetCamera() const {
    if (!camera) {
        throw std::runtime_error("Camera is null");
    }

    return *camera;
}

const BVH& Scene::GetBVH() const {
    if (!bvh) {
        throw std::runtime_error("BVH is null");
    }

    return *bvh;
}

// JSON library returns array of floats as std::vector,
// this function just reformats that as glm::mat4
glm::mat4 Scene::MatrixFromVector(const std::vector<float>& vector) const {
    glm::mat4 matrix(1.f);

    uint8_t j = 0;
    for (uint8_t i = 0; i < 4; ++i) {
        matrix[i] =
            glm::vec4(vector[j], vector[j + 1], vector[j + 2], vector[j + 3]);
        j += 4;
    }

    return matrix;
}

TriMeshPtr Scene::LoadMeshFromFile(const std::string& filePath,
                                   glm::mat4& objectToWorld,
                                   MaterialPtr&& material, uint32_t meshID, uint8_t priority) const {
    std::ifstream file;
    uint32_t numFaces;
    file.open(filePath);

    if (!file) {
        std::cerr << "Error: Mesh file " << filePath
                  << " could not be opened.\n";
        return nullptr;
    }

    file >> numFaces;

    std::vector<uint32_t> faces(numFaces);
    uint32_t maxNormIndex = 0;
    uint32_t maxUVIndex = 0;

    uint32_t numVertIndices = 0;

    for (uint32_t i = 0; i < numFaces; ++i) {
        if (!(file >> faces[i])) {
            std::cerr << "Error: Mesh file could not be read.\n";
            return nullptr;
        }
        numVertIndices += faces[i];
    }

    std::vector<uint32_t> vertIndices(numVertIndices);
    uint32_t maxVertIndex = 0;

    uint32_t k = 0;

    for (uint32_t i = 0; i < numFaces; ++i) {
        for (uint32_t j = 0; j < faces[i]; ++j) {
            uint32_t vertIndex;

            if (!(file >> vertIndex)) {
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

    for (uint32_t i = 0; i < numVertCoords; ++i) {
        float vertCoord;
        if (!(file >> vertCoord)) {
            std::cerr << "Error: Mesh file could not be read.\n";
            return nullptr;
        }

        vertCoords[i] = vertCoord;
    }

    std::vector<uint32_t> normIndices(numVertIndices);
    k = 0;

    for (uint32_t i = 0; i < numFaces; ++i) {
        for (uint32_t j = 0; j < faces[i]; ++j) {
            uint32_t normIndex;

            if (!(file >> normIndex)) {
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

    for (uint32_t i = 0; i < numNormCoords; ++i) {
        float normCoord;
        if (!(file >> normCoord)) {
            std::cerr << "Error: Mesh file could not be read.\n";
            return nullptr;
        }

        normCoords[i] = normCoord;
    }

    // UVs (optional)
    std::vector<uint32_t> UVIndices(numVertIndices);
    k = 0;
    bool noUVs = false;

    for (uint32_t i = 0; i < numFaces; ++i) {
        for (uint32_t j = 0; j < faces[i]; ++j) {
            uint32_t UVIndex;

            if (!(file >> UVIndex)) {
                if (i == 0) {
                    // No UVs present
                    noUVs = true;
                    break;
                }

                else {
                    std::cerr << "Error: Mesh file could not be read.\n";
                    return nullptr;
                }
            }

            UVIndices[k] = UVIndex;
            k += 1;

            maxUVIndex = glm::max(maxUVIndex, UVIndex);
        }

        if (noUVs) {
            break;
        }
    }

    uint32_t numUVCoords = 0;
    std::vector<float> UVCoords(0);

    if (!noUVs) {
        numUVCoords = (maxUVIndex + 1) * 2;
        UVCoords.reserve(numUVCoords);

        for (uint32_t i = 0; i < numUVCoords; ++i) {
            float UVCoord;
            if (!(file >> UVCoord)) {
                std::cerr << "Error: Mesh file could not be read.\n";
                return nullptr;
            }

            UVCoords.push_back(UVCoord);
        }
    }

    // Finished reading mesh file

    // Start building mesh

    // Calculate number of tris
    uint32_t numTris = 0;
    for (uint32_t i = 0; i < numFaces; ++i) {
        numTris += faces[i] - 2;
    }

    // Create vector of vertices
    std::vector<glm::vec3> verts(maxVertIndex + 1);
    uint32_t j = 0;
    for (uint32_t i = 0; i < maxVertIndex + 1; ++i) {
        verts[i] =
            glm::transpose(objectToWorld) *
            glm::vec4(vertCoords[j], vertCoords[j + 1], vertCoords[j + 2], 1.f);
        j += 3;
    }

    // Create vector of normals
    std::vector<glm::vec3> norms(maxNormIndex + 1);
    j = 0;
    for (uint32_t i = 0; i < maxNormIndex + 1; ++i) {
        // Need to transform normals by inverse transpose - otherwise the normal
        // will get scaled by any non-uniform scaling
        norms[i] =
            glm::normalize(glm::vec3(glm::inverse(((objectToWorld))) *
                                     glm::vec4(normCoords[j], normCoords[j + 1],
                                               normCoords[j + 2], 0.f)));
        j += 3;
    }

    // Create vector of UVs
    std::vector<glm::vec2> UVs(maxUVIndex + 1);
    if (!noUVs) {
        j = 0;
        for (uint32_t i = 0; i < maxUVIndex + 1; ++i) {
            // Need to transform normals by inverse transpose - otherwise the
            // normal will get scaled by any non-uniform scaling
            UVs[i] = glm::vec2(UVCoords[j], UVCoords[j + 1]);
            j += 2;
        }
    }

    // Create tri indices
    std::vector<uint32_t> triIndices(numTris * 3);
    k = 0;
    uint32_t l = 0;
    for (uint32_t i = 0; i < numFaces; ++i) {
        for (uint32_t j = 0; j < faces[i] - 2; ++j) {
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
    for (uint32_t i = 0; i < numFaces; ++i) {
        for (uint32_t j = 0; j < faces[i] - 2; ++j) {
            triNormIndices[k] = normIndices[l];
            triNormIndices[k + 1] = normIndices[l + j + 1];
            triNormIndices[k + 2] = normIndices[l + j + 2];
            k += 3;
        }
        l += faces[i];
    }

    // Create tri UV indices
    std::vector<uint32_t> triUVIndices(numTris * 3);
    if (!noUVs) {
        k = 0;
        l = 0;
        for (uint32_t i = 0; i < numFaces; ++i) {
            for (uint32_t j = 0; j < faces[i] - 2; ++j) {
                triUVIndices[k] = UVIndices[l];
                triUVIndices[k + 1] = UVIndices[l + j + 1];
                triUVIndices[k + 2] = UVIndices[l + j + 2];
                k += 3;
            }
            l += faces[i];
        }
    }

    std::vector<Triangle> triangles;
    triangles.reserve(numTris);

    // Now we can finally build our triangles
    j = 0;
    k = 0;
    for (uint32_t i = 0; i < numTris; ++i) {
        if (noUVs) {
            triangles.emplace_back(
                verts[triIndices[j]], verts[triIndices[j + 1]],
                verts[triIndices[j + 2]], norms[triNormIndices[j]],
                norms[triNormIndices[j + 1]], norms[triNormIndices[j + 2]]);
        }

        else {
            triangles.emplace_back(
                verts[triIndices[j]], verts[triIndices[j + 1]],
                verts[triIndices[j + 2]], norms[triNormIndices[j]],
                norms[triNormIndices[j + 1]], norms[triNormIndices[j + 2]],
                UVs[triUVIndices[k + 0]], UVs[triUVIndices[k + 1]],
                UVs[triUVIndices[k + 2]]);
        }

        j += 3;
        k += 3;
    }

    return std::make_unique<TriMesh>(std::move(triangles), std::move(material), meshID, priority);
}

PatternPtr Scene::GetRho_d(const nlohmann::json& material) {
    PatternPtr rho_dPtn;

    // rho_d
    if (material["rho_d"].is_object()) {
        std::string ptnType = material["rho_d"]["type"].get<std::string>();

        if (ptnType == "constant") {
            std::vector<float> rho_dGet =
                material["rho_d"]["value"].get<std::vector<float>>();
            glm::vec3 rho_d =
                glm::vec3(glm::min(rho_dGet[0], 1.f - glm::epsilon<float>()),
                          glm::min(rho_dGet[1], 1.f - glm::epsilon<float>()),
                          glm::min(rho_dGet[2], 1.f - glm::epsilon<float>()));

            rho_dPtn = std::make_unique<ConstantPattern>(rho_d);
        }

        if (ptnType == "texture") {
            std::string rho_dPath =
                material["rho_d"]["filePath"].get<std::string>();

            rho_dPtn = std::make_unique<TexturePattern>(rho_dPath);
        }

        else {
            std::cerr << "Error: '" << ptnType
                      << "' is not a pattern type.\nAborting.\n";
            std::abort();
        };
    }

    else {
        std::vector<float> rho_dGet =
            material["rho_d"].get<std::vector<float>>();
        glm::vec3 rho_d = glm::vec3(rho_dGet[0], rho_dGet[1], rho_dGet[2]);

        rho_dPtn = std::make_unique<ConstantPattern>(rho_d);
    }

    return rho_dPtn;
}

PatternPtr Scene::GetRho_s(const nlohmann::json& material) {
    PatternPtr rho_sPtn;

    if (material["rho_s"].is_object()) {
        std::string ptnType = material["rho_s"]["type"].get<std::string>();

        if (ptnType == "constant") {
            std::vector<float> rho_sGet =
                material["rho_s"]["value"].get<std::vector<float>>();
            glm::vec3 rho_s =
                glm::vec3(glm::min(rho_sGet[0], 1.f - glm::epsilon<float>()),
                          glm::min(rho_sGet[1], 1.f - glm::epsilon<float>()),
                          glm::min(rho_sGet[2], 1.f - glm::epsilon<float>()));

            rho_sPtn = std::make_unique<ConstantPattern>(rho_s);
        }

        if (ptnType == "texture") {
            std::string rho_sPath =
                material["rho_s"]["filePath"].get<std::string>();

            rho_sPtn = std::make_unique<TexturePattern>(rho_sPath);
        }

        else {
            std::cerr << "Error: '" << ptnType
                      << "' is not a pattern type.\nAborting.\n";
            std::abort();
        };
    }

    else {
        std::vector<float> rho_sGet =
            material["rho_s"].get<std::vector<float>>();
        glm::vec3 rho_s =
            glm::vec3(glm::min(rho_sGet[0], 1.f - glm::epsilon<float>()),
                      glm::min(rho_sGet[1], 1.f - glm::epsilon<float>()),
                      glm::min(rho_sGet[2], 1.f - glm::epsilon<float>()));

        rho_sPtn = std::make_unique<ConstantPattern>(rho_s);
    }

    return rho_sPtn;
}

PatternPtr Scene::GetEta(const nlohmann::json& material) {
    PatternPtr etaPtn;

    if (material["eta"].is_object()) {
        std::string ptnType = material["eta"]["type"].get<std::string>();

        if (ptnType == "constant") {
            float eta = material["eta"]["value"].get<float>();

            etaPtn = std::make_unique<ConstantPattern>(glm::vec3(eta));
        }

        if (ptnType == "texture") {
            std::string etaPath =
                material["eta"]["filePath"].get<std::string>();

            etaPtn = std::make_unique<TexturePattern>(etaPath);
        }

        else {
            std::cerr << "Error: '" << ptnType
                      << "' is not a pattern type.\nAborting.\n";
            std::abort();
        };
    }

    else {
        float eta = material["eta"].get<float>();

        etaPtn = std::make_unique<ConstantPattern>(glm::vec3(eta));
    }

    return etaPtn;
}

PatternPtr Scene::GetTau(const nlohmann::json& material) {
    PatternPtr tauPtn;

    if (material["tau"].is_object()) {
        std::string ptnType = material["tau"]["type"].get<std::string>();

        if (ptnType == "constant") {
            std::vector<float> tauGet =
                material["tau"]["value"].get<std::vector<float>>();
            glm::vec3 tau =
                glm::vec3(glm::min(tauGet[0], 1.f - glm::epsilon<float>()),
                          glm::min(tauGet[1], 1.f - glm::epsilon<float>()),
                          glm::min(tauGet[2], 1.f - glm::epsilon<float>()));

            tauPtn = std::make_unique<ConstantPattern>(tau);
        }

        if (ptnType == "texture") {
            std::string tauPath =
                material["tau"]["filePath"].get<std::string>();

            tauPtn = std::make_unique<TexturePattern>(tauPath);
        }

        else {
            std::cerr << "Error: '" << ptnType
                      << "' is not a pattern type.\nAborting.\n";
            std::abort();
        };
    }

    else {
        std::vector<float> tauGet = material["tau"].get<std::vector<float>>();
        glm::vec3 tau =
            glm::vec3(glm::min(tauGet[0], 1.f - glm::epsilon<float>()),
                      glm::min(tauGet[1], 1.f - glm::epsilon<float>()),
                      glm::min(tauGet[2], 1.f - glm::epsilon<float>()));

        tauPtn = std::make_unique<ConstantPattern>(tau);
    }

    return tauPtn;
}

PatternPtr Scene::GetAlpha(const nlohmann::json& material) {
    PatternPtr alphaPtn;

    if (material["roughness"].is_object()) {
        std::string ptnType = material["roughness"]["type"].get<std::string>();

        if (ptnType == "constant") {
            float roughness = material["roughness"]["value"].get<float>();

            alphaPtn = std::make_unique<ConstantPattern>(
                glm::vec3(roughness * roughness));
        }

        if (ptnType == "texture") {
            std::string roughnessPath =
                material["roughness"]["filePath"].get<std::string>();

            alphaPtn = std::make_unique<TexturePattern>(roughnessPath, true);
        }

        else {
            std::cerr << "Error: '" << ptnType
                      << "' is not a pattern type.\nAborting.\n";
            std::abort();
        };
    }

    else {
        float alpha = material["roughness"].get<float>();

        alphaPtn = std::make_unique<ConstantPattern>(glm::vec3(alpha * alpha));
    }

    return alphaPtn;
}

// TODO: Add option to turn off PDF?
PatternPtr Scene::GetLe(const nlohmann::json& light) {
    PatternPtr LePtn;

    if (light["Le"].is_object()) {
        std::string ptnType = light["Le"]["type"].get<std::string>();

        if (ptnType == "constant") {
            std::vector<float> LeGet =
                light["Le"]["value"].get<std::vector<float>>();
            glm::vec3 Le =
                glm::vec3(glm::min(LeGet[0], 1.f - glm::epsilon<float>()),
                          glm::min(LeGet[1], 1.f - glm::epsilon<float>()),
                          glm::min(LeGet[2], 1.f - glm::epsilon<float>()));

            LePtn = std::make_unique<ConstantPattern>(Le);
        }

        if (ptnType == "texture") {
            std::string filePath = light["Le"]["filePath"].get<std::string>();

            LePtn = std::make_unique<TexturePattern>(filePath, 0, 1);
        }

        else {
            std::cerr << "Error: '" << ptnType
                      << "' is not a pattern type.\nAborting.\n";
            std::abort();
        };
    }

    else {
        std::vector<float> LeGet = light["Le"].get<std::vector<float>>();
        glm::vec3 Le =
            glm::vec3(glm::min(LeGet[0], 1.f - glm::epsilon<float>()),
                      glm::min(LeGet[1], 1.f - glm::epsilon<float>()),
                      glm::min(LeGet[2], 1.f - glm::epsilon<float>()));

        LePtn = std::make_unique<ConstantPattern>(Le);
    }

    return LePtn;
}

PatternPtr Scene::GetNormal(const nlohmann::json& material) {
    PatternPtr normalPtn;

    if (material.contains("normal")) {
        if (material["normal"].is_object()) {
            std::string ptnType = material["normal"]["type"].get<std::string>();

            if (ptnType == "constant") {
                std::vector<float> normalGet =
                    material["normal"]["value"].get<std::vector<float>>();
                glm::vec3 normal = glm::vec3(
                    glm::min(normalGet[0], 1.f - glm::epsilon<float>()),
                    glm::min(normalGet[1], 1.f - glm::epsilon<float>()),
                    glm::min(normalGet[2], 1.f - glm::epsilon<float>()));

                normalPtn = std::make_unique<ConstantPattern>(normal);
            }

            if (ptnType == "texture") {
                std::string filePath =
                    material["normal"]["filePath"].get<std::string>();

                normalPtn = std::make_unique<TexturePattern>(filePath);
            }

            else {
                std::cerr << "Error: '" << ptnType
                          << "' is not a pattern type.\nAborting.\n";
                std::abort();
            };
        }

        else if (material["normal"].is_array()) {
            std::vector<float> normalGet =
                material["normal"].get<std::vector<float>>();
            glm::vec3 normal =
                glm::vec3(glm::min(normalGet[0], 1.f - glm::epsilon<float>()),
                          glm::min(normalGet[1], 1.f - glm::epsilon<float>()),
                          glm::min(normalGet[2], 1.f - glm::epsilon<float>()));

            normalPtn = std::make_unique<ConstantPattern>(normal);
        }

        else {
            return nullptr;
        }

        return normalPtn;
    }
    return nullptr;
}

void Scene::LoadMeshes(nlohmann::json& json) {
    std::vector<std::string> meshFilePaths;
    std::vector<std::unique_ptr<Material>> meshMaterials;
    std::vector<uint8_t> meshPriorities;
    std::vector<glm::mat4> meshTransforms;

    if (!json["meshes"].is_null()) {
        uint8_t numMeshes = json["meshes"].size();

        meshFilePaths.reserve(numMeshes);
        meshMaterials.reserve(numMeshes);
        meshPriorities.reserve(numMeshes);
        meshTransforms.reserve(numMeshes);

        // Load materials
        for (auto& elem : json["meshes"]) {
            std::string filePath = "";
            try {
                meshFilePaths.push_back(elem["filePath"].get<std::string>());
            } catch (nlohmann::json::exception& e) {
                std::cerr << "Error in mesh file path: " << e.what() << "\n";
            }

            auto& mat = elem["material"];
            try {
                std::string type = mat["type"].get<std::string>();

                if (type == "lambert") {
                    PatternPtr rho_dPtn = GetRho_d(mat);
                    PatternPtr normalPtn = GetNormal(mat);

                    meshMaterials.push_back(std::make_unique<DiffuseMaterial>(
                        std::move(rho_dPtn), std::move(normalPtn)));
                }

                else if (type == "specular") {
                    PatternPtr rho_sPtn = GetRho_s(mat);
                    PatternPtr etaPtn = GetEta(mat);
                    PatternPtr normalPtn = GetNormal(mat);

                    meshMaterials.push_back(std::make_unique<SpecularMaterial>(
                        std::move(rho_sPtn), std::move(etaPtn),
                        std::move(normalPtn)));
                }

                else if (type == "glass") {
                    PatternPtr rho_sPtn = GetRho_s(mat);
                    PatternPtr tauPtn = GetTau(mat);
                    PatternPtr etaPtn = GetEta(mat);
                    PatternPtr alphaPtn = GetAlpha(mat);
                    PatternPtr normalPtn = GetNormal(mat);

                    meshMaterials.push_back(std::make_unique<GlassMaterial>(
                        std::move(rho_sPtn), std::move(tauPtn),
                        std::move(etaPtn), std::move(alphaPtn),
                        std::move(normalPtn)));
                }

                else if (type == "glossy") {
                    PatternPtr rho_sPtn = GetRho_s(mat);
                    PatternPtr etaPtn = GetEta(mat);
                    PatternPtr alphaPtn = GetAlpha(mat);
                    PatternPtr normalPtn = GetNormal(mat);

                    meshMaterials.push_back(
                        std::make_unique<GlossyDielectricMaterial>(
                            std::move(rho_sPtn), std::move(etaPtn),
                            std::move(alphaPtn), std::move(normalPtn)));
                }

                else if (type == "plastic") {
                    PatternPtr rho_dPtn = GetRho_d(mat);
                    PatternPtr rho_sPtn = GetRho_s(mat);
                    PatternPtr etaPtn = GetEta(mat);
                    PatternPtr alphaPtn = GetAlpha(mat);
                    PatternPtr normalPtn = GetNormal(mat);

                    meshMaterials.push_back(std::make_unique<PlasticMaterial>(
                        std::move(rho_dPtn), std::move(rho_sPtn),
                        std::move(etaPtn), std::move(alphaPtn),
                        std::move(normalPtn)));
                }

                else {
                    std::cerr << "Error: '" << type
                              << "' is not a material type.\nAborting.\n";
                    std::abort();
                }
            } catch (nlohmann::json::exception& e) {
                std::cerr << "Error in mesh material type: " << e.what()
                          << "\n";
            }
        }

        // Load priorities
        for (auto& elem : json["meshes"]) {
            uint8_t priority = 0;

            if (!elem["priority"].is_null()) {
                try {
                    priority = elem["priority"].get<uint8_t>();
                } catch (nlohmann::json::exception& e) {
                    std::cerr << "Error in mesh transform: " << e.what()
                              << "\n";
                }
            }

            meshPriorities.push_back(priority);
        }

        // Load transforms
        for (auto& elem : json["meshes"]) {
            std::vector<float> meshTransform = {1.f, 0.f, 0.f, 0.f, 0.f, 1.f,
                                                0.f, 0.f, 0.f, 0.f, 1.f, 0.f,
                                                0.f, 0.f, 0.f, 1.f};
            try {
                meshTransform = elem["transform"].get<std::vector<float>>();
            } catch (nlohmann::json::exception& e) {
                std::cerr << "Error in mesh transform: " << e.what() << "\n";
            }

            glm::mat4 objectToWorld = MatrixFromVector(meshTransform);
            meshTransforms.push_back(objectToWorld);
        }
    }

    else {
        std::cerr << "Error: No meshes found.\n";
    }

    std::vector<TriMeshPtr> meshes;
    meshes.reserve(meshFilePaths.size());

    for (uint32_t i = 0; i < meshFilePaths.size(); ++i) {
        meshes.push_back(LoadMeshFromFile(meshFilePaths[i], meshTransforms[i],
                                          std::move(meshMaterials[i]), i, meshPriorities[i]));
    }

    bvh = std::make_unique<BVH>(std::move(meshes));
}

void Scene::LoadCamera(const nlohmann::json& json) {
    glm::mat4 cameraToWorld(1.f);
    float fov = 1.f;
    std::vector<float> camTransform = {1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f,
                                       0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f};

    try {
        fov = json["camera"]["fov"].get<float>();
        camTransform = json["camera"]["transform"].get<std::vector<float>>();
    }

    catch (nlohmann::json::exception& e) {
        std::cerr << "Error in camera: " << e.what() << "\n";
    }

    cameraToWorld = MatrixFromVector(camTransform);

    camera = std::make_unique<PinholeCamera>(fov, cameraToWorld);
}

void Scene::LoadLights(const nlohmann::json& json) {
    if (!json["lights"].is_null()) {
        uint8_t numLights = json["lights"].size();
        lights.reserve(numLights);

        for (auto& elem : json["lights"]) {
            std::vector<float> lightTransform = {1.f, 0.f, 0.f, 0.f, 0.f, 1.f,
                                                 0.f, 0.f, 0.f, 0.f, 1.f, 0.f,
                                                 0.f, 0.f, 0.f, 1.f};
            try {
                lightTransform = elem["transform"].get<std::vector<float>>();
            } catch (nlohmann::json::exception& e) {
                std::cerr << "Error in light transform: " << e.what() << "\n";
            }

            glm::mat4 lightToWorld = MatrixFromVector(lightTransform);

            try {
                std::string type = elem["type"].get<std::string>();

                if (type == "disk") {
                    float radius = elem["radius"].get<float>();
                    PatternPtr LePtn = GetLe(elem);
                    float intensity = elem["intensity"].get<float>();

                    lights.emplace_back(std::make_unique<DiskLight>(
                        radius, std::move(LePtn), intensity, lightToWorld));
                }

                if (type == "ring") {
                    float radius = elem["radius"].get<float>();
                    float innerRadius = elem["innerRadius"].get<float>();
                    PatternPtr LePtn = GetLe(elem);
                    float intensity = elem["intensity"].get<float>();

                    lights.emplace_back(std::make_unique<RingLight>(
                        radius, innerRadius, std::move(LePtn), intensity,
                        lightToWorld));
                }

                if (type == "environment") {
                    PatternPtr LePtn = GetLe(elem);
                    float intensity = elem["intensity"].get<float>();

                    lights.emplace_back(std::make_unique<EnvironmentLight>(
                        std::move(LePtn), intensity, lightToWorld));
                }
            }

            catch (nlohmann::json::exception& e) {
                std::cerr << "Error in mesh material type: " << e.what()
                          << "\n";
            }
        }
    }
}
