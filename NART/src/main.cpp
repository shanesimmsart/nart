#include <iostream>
#include <random>
#include <chrono>
#include <fstream>
#include <queue>
#include <tbb/task_group.h>
#include <OpenEXR/ImfRgbaFile.h>
#include <OpenEXR/ImfArray.h>
#include <glm/common.hpp>
#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/trigonometric.hpp>
#include <glm/geometric.hpp>
#include <nlohmann/json.hpp>

static const float Pi = 3.14159265358979323846;
static const float OneOverPi = 0.31830988618379067154;
static const float OneMinusEpsilon = 0x1.fffffep-1;

#define Infinity std::numeric_limits<float>::infinity()
#define FilterTableResolution 64

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

class Ray
{
public:
    Ray(glm::vec3 o, glm::vec3 d) : o(o), d(d) {}

    glm::vec3 o = glm::vec3(0.f, 0.f, 0.f);
    glm::vec3 d = glm::vec3(0.f, 1.f, 0.f);
};

class Light
{
public:
    Light(glm::vec3 Le, glm::mat4 LightToWorld) : Le(Le), LightToWorld(LightToWorld) {}

    // Return incident radiance from light, set wi
    virtual glm::vec3 Li(glm::vec3 p, glm::vec3* wi) = 0;

    // Does light have a delta distribution?
    bool isDelta = false;

protected:
    // Radiance emitted
    glm::vec3 Le;
    glm::mat4 LightToWorld;
};

class DistantLight : public Light
{
public:
    DistantLight(glm::vec3 Le, glm::mat4 LightToWorld) : Light(Le, LightToWorld)
    {
        isDelta = true;
        // Pointing down by default
        direction = glm::vec4(0.f, 0.f, -1.f, 0.f) * LightToWorld;
    }

    glm::vec3 Li(glm::vec3 p, glm::vec3* wi)
    {
        *wi = -direction;
        return Le;
    }
    
    glm::vec3 direction;
};

class BxDF
{
public:
    // BxDF - wo and wi must always be in local coord sys
    virtual glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi) = 0;

    // Sample BRDF, setting wi and pdf
    // virtual glm::vec3 sample_f(glm::vec3 wo, glm::vec3* wi, float *pdf) = 0;

protected:
    BxDF() {}
};

class BSDF
{
public:
    BSDF(glm::vec3 n, uint8_t numBxDFs) : n(n), numBxDFs(numBxDFs)
    {
        // Build local coord sys from normal
        if (n.x > n.y)
        {
            nt = glm::normalize(glm::normalize(glm::vec3(n.z, 0.f, -n.x)));
        }

        else
        {
            nt = glm::normalize(glm::normalize(glm::vec3(0.f, n.z, -n.y)));
        }
        nb = glm::cross(nt, n);

        bxdfs.reserve(numBxDFs);
    }

    glm::vec3 ToLocal(glm::vec3 v)
    {
        return glm::vec3(glm::dot(v, nt), glm::dot(v, nb), glm::dot(v, n));
    }

    void AddBxDF(std::shared_ptr<BxDF> bxdf)
    {
        bxdfs.emplace_back(bxdf);
    }

    // Sum BxDFs
    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi)
    {
        glm::vec3 f = glm::vec3(0.f);
        for (uint8_t i = 0; i < numBxDFs; ++i)
        {
            f += bxdfs[i]->f(wo, wi);
        }
        return f;
    }

    // Sample BSDF, setting wi and pdf
    // virtual glm::vec3 sample_f(glm::vec3 wo, glm::vec3* wi, float *pdf) = 0;

protected:
    // Coord sys in worldspace
    glm::vec3 nt, nb, n;
    uint8_t numBxDFs = 0;
    std::vector<std::shared_ptr<BxDF>> bxdfs;
};

class LambertBRDF : public BxDF
{
public:
    LambertBRDF(glm::vec3 rho) : rho(rho) {}

    glm::vec3 f(const glm::vec3& wo, const glm::vec3& wi)
    {
        return rho * OneOverPi;
    }

private:
    // Albedo
    glm::vec3 rho;
};

class Material
{
public:
    virtual BSDF CreateBSDF(glm::vec3 n) = 0;

protected:
    Material() {};
};

class DiffuseMaterial : public Material
{
public:
    DiffuseMaterial(glm::vec3 rho) : rho(rho) // your boat...
    {}

    BSDF CreateBSDF(glm::vec3 n)
    {
        BSDF bsdf(n, 1);
        bsdf.AddBxDF(std::make_shared<LambertBRDF>(rho));
        return bsdf;
    }

private:
    glm::vec3 rho;
};

struct Intersection
{
    std::shared_ptr<Material> material;
    // Barycentric coords
    float u, v;
    // Point, geometric normal, shading normal
    glm::vec3 p, gn, sn;

    // Upper and lower bounds of ray to be considered
    float tMin = 0.f;
    float tMax = Infinity;
};

class Triangle
{
public:
    Triangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2) :
        v0(v0), v1(v1), v2(v2), n0(n0), n1(n1), n2(n2) {}
    bool Intersect(const Ray& ray, Intersection* isect)
    {
        // Un-normalised direction of surface
        glm::vec3 n = glm::cross((v2 - v1), (v0 - v1));
        // Length of cross product == area of triangle bounded by vectors * 2
        float area = glm::length(n) / 2.f;

        // Find t of intersection with plane
        float t = (glm::dot(v0, n) - glm::dot(ray.o, n)) / glm::dot(ray.d, n);

        if (t <= isect->tMin || t >= isect->tMax) return false;

        glm::vec3 p = ray.o + t * ray.d;

        // Check if p is on the correct side of each vector, and compute barycentric coords
        glm::vec3 v0v1 = v1 - v0;
        glm::vec3 v0p = p - v0;
        if (glm::dot((glm::cross(v0v1, v0p)), n) < 0.f) return false;

        glm::vec3 v1v2 = v2 - v1;
        glm::vec3 v1p = p - v1;
        glm::vec3 c1 = glm::cross(v1v2, v1p);
        if (glm::dot(c1, n) < 0.f) return false;
        float u = glm::length(c1) / 2.f / area;

        glm::vec3 v2v0 = v0 - v2;
        glm::vec3 v2p = p - v2;
        glm::vec3 c2 = glm::cross(v2v0, v2p);
        if (glm::dot(c2, n) < 0.f) return false;

        isect->tMax = t;
        isect->u = u;
        isect->v = glm::length(c2) / 2.f / area;
        isect->gn = glm::normalize(n);
        isect->sn = n0 * isect->u + n1 * isect->v + n2 * (1 - isect->u - isect->v);
        isect->p = p;

        return true;
    }

    const glm::vec3 v0, v1, v2;
    const glm::vec3 n0, n1, n2;
};

class TriMesh
{
public:
    TriMesh(std::shared_ptr<Material> material) : material(material) {}

    // Triangles are passed directly, outside of the ctor, to reduce memory allocations
    std::vector<Triangle> triangles;
    std::shared_ptr<Material> material;
};

struct BoundingVolume
{
    BoundingVolume()
    {
        // Bounding volumes are defined by slabs along the cardinal direction vectors
        // and diagonally oriented vectors
        float OneOverSqrt3 = 1.f / glm::sqrt(3);
        normals[0] = glm::vec3(1, 0, 0);
        normals[1] = glm::vec3(0, 1, 0);
        normals[2] = glm::vec3(0, 0, 1);
        normals[3] = glm::vec3(OneOverSqrt3, OneOverSqrt3, OneOverSqrt3);
        normals[4] = glm::vec3(OneOverSqrt3, -OneOverSqrt3, OneOverSqrt3);
        normals[5] = glm::vec3(-OneOverSqrt3, OneOverSqrt3, OneOverSqrt3);
        normals[6] = glm::vec3(-OneOverSqrt3, -OneOverSqrt3, OneOverSqrt3);

        for (uint8_t i = 0; i < 7; ++i)
        {
            boundsMin[i] = Infinity;
            boundsMax[i] = -Infinity;
        }
    }

    // Extend bounding volume by another to ensure it is fully contained
    void ExtendBy(const BoundingVolume& bv)
    {
        for (uint8_t i = 0; i < 7; ++i)
        {
            boundsMin[i] = glm::min(boundsMin[i], bv.boundsMin[i]);
            boundsMax[i] = glm::max(boundsMax[i], bv.boundsMax[i]);
        }
    }

    bool Intersect(const Ray& ray, Intersection* isect)
    {
        float tMin = -Infinity;
        float tMax = Infinity;
        float dotLow = Infinity;

        for (uint8_t i = 0; i < 7; ++i)
        {
            // Equation of a plane:
            // (o+td)dotN - bound = 0
            // OdotN + tDdotN = bound
            // tDdotN = bound - OdotN
            // t = (bound - OdotN) / DdotN

            float slabMin = (boundsMin[i] - glm::dot(ray.o, normals[i])) / glm::dot(ray.d, normals[i]);
            float slabMax = (boundsMax[i] - glm::dot(ray.o, normals[i])) / glm::dot(ray.d, normals[i]);

            if (slabMin > slabMax)
            {
                std::swap(slabMin, slabMax);
            }

            if (slabMin > tMax || tMin > slabMax)
            {
                return false;
            }

            if (glm::dot(ray.o + (slabMin * ray.d), normals[i]) < dotLow)
            {
                dotLow = glm::dot(ray.o + (slabMin * ray.d), normals[i]);
            }

            tMin = glm::max(tMin, slabMin);
            tMax = glm::min(tMax, slabMax);
        }

        isect->tMin = tMin;
        isect->tMax = tMax;

        return true;
    }

    glm::vec3 normals[7];
    float boundsMin[7];
    float boundsMax[7];
};

// A chunk of triangles to be placed inside of a bounding volume
class Chunk
{
    // Stores a mesh pointer and an index into one of it's triangles
    struct TriangleIndex
    {
        TriangleIndex(std::shared_ptr<TriMesh> mesh, const uint32_t& index) : mesh(mesh), index(index) {}
        std::shared_ptr<TriMesh> mesh;
        const uint32_t index;
    };

public:
    void Insert(std::shared_ptr<TriMesh> mesh, const uint32_t& index)
    {
        triangleIndices.push_back(TriangleIndex(mesh, index));
    }

    bool Intersect(const Ray& ray, Intersection* isect) const
    {
        bool hit = false;
        for (uint32_t i = 0; i < triangleIndices.size(); ++i)
        {
            if (triangleIndices[i].mesh->triangles[triangleIndices[i].index].Intersect(ray, isect))
            {
                hit = true;
                isect->material = triangleIndices[i].mesh->material;
            }
        }
        return hit;
    }

    void CalculateBounds()
    {
        float OneOverSqrt3 = 1.f / glm::sqrt(3);

        glm::vec3 normals[7];
        normals[0] = glm::vec3(1, 0, 0);
        normals[1] = glm::vec3(0, 1, 0);
        normals[2] = glm::vec3(0, 0, 1);
        normals[3] = glm::vec3(OneOverSqrt3, OneOverSqrt3, OneOverSqrt3);
        normals[4] = glm::vec3(OneOverSqrt3, -OneOverSqrt3, OneOverSqrt3);
        normals[5] = glm::vec3(-OneOverSqrt3, OneOverSqrt3, OneOverSqrt3);
        normals[6] = glm::vec3(-OneOverSqrt3, -OneOverSqrt3, OneOverSqrt3);

        for (TriangleIndex triIndex : triangleIndices)
        {
            Triangle triangle = triIndex.mesh->triangles[triIndex.index];

            bboxMin = glm::min(bboxMin, triangle.v0);
            bboxMax = glm::max(bboxMax, triangle.v0);

            bboxMin = glm::min(bboxMin, triangle.v1);
            bboxMax = glm::max(bboxMax, triangle.v1);

            bboxMin = glm::min(bboxMin, triangle.v2);
            bboxMax = glm::max(bboxMax, triangle.v2);

            for (uint8_t i = 0; i < 7; ++i)
            {
                float v0DotNorm = glm::dot(triangle.v0, normals[i]);
                bv.boundsMin[i] = glm::min(v0DotNorm, bv.boundsMin[i]);
                bv.boundsMax[i] = glm::max(v0DotNorm, bv.boundsMax[i]);

                float v1DotNorm = glm::dot(triangle.v1, normals[i]);
                bv.boundsMin[i] = glm::min(v1DotNorm, bv.boundsMin[i]);
                bv.boundsMax[i] = glm::max(v1DotNorm, bv.boundsMax[i]);

                float v2DotNorm = glm::dot(triangle.v2, normals[i]);
                bv.boundsMin[i] = glm::min(v2DotNorm, bv.boundsMin[i]);
                bv.boundsMax[i] = glm::max(v2DotNorm, bv.boundsMax[i]);
            }
        }
    }

    // Bounding box (used for octree position)
    glm::vec3 bboxMin = glm::vec3(Infinity);
    glm::vec3 bboxMax = glm::vec3(-Infinity);
    // Bounding volume (used for intersecting)
    BoundingVolume bv = BoundingVolume();

private:
    std::vector<TriangleIndex> triangleIndices;
};

// BVH uses an octree structure for containing bounding volumes
class Octree
{
    struct OctreeNode
    {
        OctreeNode(glm::vec3 nodeMin, glm::vec3 nodeMax) : nodeMin(nodeMin), nodeMax(nodeMax) {}

        std::vector<std::shared_ptr<Chunk>> chunks;
        std::shared_ptr<OctreeNode> children[8] = { nullptr };
        bool isLeaf = true;

        const glm::vec3 nodeMin;
        const glm::vec3 nodeMax;

        BoundingVolume bv = BoundingVolume();
    };

public:
    Octree(glm::vec3 sceneMin, glm::vec3 sceneMax) : sceneMin(sceneMin), sceneMax(sceneMax)
    {
        root = std::make_shared<OctreeNode>(sceneMin, sceneMax);
    }

    void Insert(std::shared_ptr<Chunk> chunk)
    {
        chunk->CalculateBounds();
        InsertChunkIntoNode(chunk, root);
    }

    void Build()
    {
        BuildBoundingVolumes(root);
    }

    bool Intersect(const Ray& ray, Intersection* isect)
    {
        // Using priority queue to keep track of closest node intersection
        std::priority_queue<std::pair<float, std::shared_ptr<OctreeNode>>, std::vector<std::pair<float, std::shared_ptr<OctreeNode>>>, std::greater<std::pair<float, std::shared_ptr<OctreeNode>>>> queue;
        queue.push(std::make_pair(Infinity, root));

        bool hit = false;

        while (!queue.empty())
        {
            std::shared_ptr<OctreeNode> topNode = queue.top().second;
            queue.pop();
            Intersection triIsect;
            for (uint8_t i = 0; i < 8; ++i)
            {
                std::shared_ptr<OctreeNode> node = topNode->children[i];
                if (node)
                {
                    Intersection bvIsect;
                    if (node->bv.Intersect(ray, &bvIsect))
                    {
                        queue.push(std::make_pair(bvIsect.tMin, topNode->children[i]));

                        if (node->isLeaf)
                        {
                            for (uint32_t j = 0; j < node->chunks.size(); ++j)
                            {
                                if (node->chunks[j]->Intersect(ray, &triIsect))
                                {
                                    hit = true;
                                }
                            }
                        }
                    }
                }
            }
            // If the nearest triangle is closer than the nearest bounding volume in the queue, we don't need to continue
            if (!queue.empty() && isect->tMax < queue.top().first) break;
            if (isect->tMax > triIsect.tMax) *isect = triIsect;
        }

        return hit;
    }

private:
    void InsertChunkIntoNode(std::shared_ptr<Chunk> chunk, std::shared_ptr<OctreeNode> node, uint8_t depth = 0)
    {
        if (node->isLeaf)
        {
            if (node->chunks.empty() || depth >= maxDepth)
            {
                node->chunks.push_back(chunk);
            }

            else
            {
                node->isLeaf = false;

                for (uint32_t i = 0; i < node->chunks.size(); ++i)
                {
                    InsertChunkIntoNode(node->chunks.back(), node, depth++);
                    node->chunks.pop_back();
                }

                InsertChunkIntoNode(chunk, node, depth++);
            }
        }

        else
        {
            // Figure out child node index based on chunk centroid and node centroid
            uint8_t nodeIndex = 0;

            glm::vec3 chunkCentroid = chunk->bboxMin + ((chunk->bboxMax - chunk->bboxMin) * 0.5f);
            glm::vec3 nodeCentroid = node->nodeMin + ((node->nodeMax - node->nodeMin) * 0.5f);

            if (chunkCentroid.x > nodeCentroid.x) nodeIndex = nodeIndex | 1;
            if (chunkCentroid.y > nodeCentroid.y) nodeIndex = nodeIndex | 2;
            if (chunkCentroid.z > nodeCentroid.z) nodeIndex = nodeIndex | 4;

            if (!node->children[nodeIndex])
            {
                glm::vec3 childNodeSize = (node->nodeMax - node->nodeMin) * 0.5f;
                glm::vec3 childNodeMin = nodeCentroid;
                glm::vec3 childNodeMax = nodeCentroid;

                // Compute child bounds
                if (nodeIndex & 1) childNodeMax.x += childNodeSize.x;
                else childNodeMin.x -= childNodeSize.x;
                if (nodeIndex & 2) childNodeMax.y += childNodeSize.y;
                else childNodeMin.y -= childNodeSize.y;
                if (nodeIndex & 4) childNodeMax.z += childNodeSize.z;
                else childNodeMin.z -= childNodeSize.z;

                node->children[nodeIndex] = std::make_shared<OctreeNode>(childNodeMin, childNodeMax);
            }

            InsertChunkIntoNode(chunk, node->children[nodeIndex], depth++);
        }
    }

    void BuildBoundingVolumes(std::shared_ptr<OctreeNode> node)
    {
        if (node->isLeaf)
        {
            for (uint32_t i = 0; i < node->chunks.size(); ++i)
            {
                node->bv.ExtendBy(node->chunks[i]->bv);
            }
        }

        else
        {
            for (uint8_t i = 0; i < 8; ++i)
            {
                if (node->children[i])
                {
                    BuildBoundingVolumes(node->children[i]);
                    node->bv.ExtendBy(node->children[i]->bv);
                }
            }
        }
    }

    uint8_t maxDepth = 6;
    std::shared_ptr<OctreeNode> root;
    const glm::vec3 sceneMin;
    const glm::vec3 sceneMax;
};

class BVH
{
public:
    BVH(std::vector<std::shared_ptr<TriMesh>> meshes) : meshes(meshes)
    {
        // Calculate scene extents and number of triangles
        uint32_t numTriangles = 0;
        glm::vec3 sceneMax(-Infinity);
        glm::vec3 sceneMin(Infinity);
        for (std::shared_ptr<TriMesh> mesh : meshes)
        {
            numTriangles += mesh->triangles.size();
            for (Triangle triangle : mesh->triangles)
            {
                sceneMax = glm::max(triangle.v0, sceneMax);
                sceneMin = glm::min(triangle.v0, sceneMin);

                sceneMax = glm::max(triangle.v1, sceneMax);
                sceneMin = glm::min(triangle.v1, sceneMin);

                sceneMax = glm::max(triangle.v2, sceneMax);
                sceneMin = glm::min(triangle.v2, sceneMin);
            }
        }

        glm::vec3 sceneSize = sceneMax - sceneMin;

        // Calculate grid resolution
        // Cleary et al. 1983. Design and analysis of a parallel ray tracing computer.
        // In their formula, lambda controls the granularity of the grid;
        // values between 3 and 5 are recommended
        float lambda = 3.f;
        float sceneVolume = sceneSize.x * sceneSize.y * sceneSize.z;
        // Calculate resolution of grid, and create a chunk per grid cell
        glm::vec3 resolution = glm::floor(sceneSize * glm::pow((glm::vec3(static_cast<float>(numTriangles)) / sceneVolume) * lambda, glm::vec3(1.f / 3.f)));

        resolution = glm::clamp(resolution, glm::vec3(1.f), glm::vec3(128.f));
        uint32_t numChunks = static_cast<uint32_t>(resolution.x * resolution.y * resolution.z);
        chunks = std::vector<std::shared_ptr<Chunk>>(numChunks + 1);

        // Add triangles to chunks
        for (uint32_t j = 0; j < meshes.size(); ++j)
        {
            for (uint32_t i = 0; i < meshes[j]->triangles.size(); ++i)
            {
                glm::vec3 triangleMin = glm::min(glm::min(meshes[j]->triangles[i].v0, meshes[j]->triangles[i].v1), meshes[j]->triangles[i].v2) - sceneMin;
                glm::vec3 chunkCoords = glm::floor((triangleMin / sceneSize) * resolution);
                uint32_t chunkIndex = static_cast<uint32_t>(glm::floor(chunkCoords.x * resolution.y * resolution.z + chunkCoords.y * resolution.z + chunkCoords.x));
                if (!chunks[chunkIndex])
                {
                    std::shared_ptr<Chunk> chunk = std::make_shared<Chunk>();
                    chunks[chunkIndex] = chunk;
                }
                chunks[chunkIndex]->Insert(meshes[j], i);
            }
        }

        // Create octree structure from chunks
        octree = std::make_unique<Octree>(sceneMin, sceneMax);

        for (uint32_t i = 0; i < chunks.size(); ++i)
        {
            if (chunks[i])
            {
                octree->Insert(chunks[i]);
            }
        }

        // Create bounding volume for each octree node, from leaves to root
        octree->Build();
    }

    bool Intersect(const Ray& ray, Intersection* isect) const
    {
        return octree->Intersect(ray, isect);
    }

    std::vector<std::shared_ptr<Chunk>> chunks;
    std::vector<std::shared_ptr<TriMesh>> meshes;
    std::shared_ptr<Octree> octree;
};

// TODO: Add material to mesh
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
        norms[i] = glm::normalize(glm::inverse(objectToWorld) * glm::vec4(normCoords[j], normCoords[j + 1], normCoords[j + 2], 0.f));
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

struct Pixel
{
    glm::vec4 contribution = glm::vec4(0.f);
    float filterWeightSum = 0.f;
};

float StratifiedSample1D(std::default_random_engine& rng, uint32_t n, uint32_t nSamples)
{
    float invNSamples = 1.f / static_cast<float>(nSamples);
    std::uniform_real_distribution<float> distribution(0.f, OneMinusEpsilon);
    return (static_cast<float>(n) + distribution(rng)) * invNSamples;
}

// Using this technique instead of stratifying in x and y independently lets have any number of spp
// including primes, e.g 2, 3, 5, 7... and have them be evenly distributed
void LatinSquare(std::default_random_engine& rng, uint32_t nSamples, std::vector<glm::vec2>& samples)
{
    // Create stratified samples along diagonal
    for (uint32_t i = 0; i < nSamples; ++i)
    {
        samples[i] = glm::vec2(StratifiedSample1D(rng, i, nSamples), StratifiedSample1D(rng, i, nSamples));
    }

    // Shuffle dimensions
    for (uint32_t i = 0; i < nSamples; ++i)
    {
        std::uniform_int_distribution<uint32_t> distribution(0, nSamples - 1 - i);
        uint32_t choice = distribution(rng);
        std::swap(samples[i].x, samples[choice].x);
        choice = distribution(rng);
        std::swap(samples[i].y, samples[choice].y);
    }
}

float Gaussian(float width, float x)
{
    if (x >= width) return 0.f;

    // In a Gaussian distribution, any value beyond 3 standard deviations is negligible (<0.3%)
    float sigma = width / 3.f;

    return (1.f / glm::sqrt(2.f * Pi * sigma * sigma)) * glm::exp(-(x * x) / (2.f * sigma * sigma));
}

void AddSample(const RenderInfo& info, const float* filterTable, glm::vec2 sampleCoords, glm::vec4 L, std::vector<Pixel>& pixels)
{
    // Calculate discrete x and y bounds of filter
    uint32_t x0 = static_cast<uint32_t>(glm::floor(sampleCoords.x - info.filterWidth));
    uint32_t x1 = static_cast<uint32_t>(glm::ceil(sampleCoords.x + info.filterWidth));
    uint32_t y0 = static_cast<uint32_t>(glm::floor(sampleCoords.y - info.filterWidth));
    uint32_t y1 = static_cast<uint32_t>(glm::ceil(sampleCoords.y + info.filterWidth));

    for (uint32_t y = y0; y < y1; ++y)
    {
        for (uint32_t x = x0; x < x1; ++x)
        {
            // Calculate distance from sampleCoords
            float distX = (static_cast<float>(x) + 0.5f) - sampleCoords.x;
            float distY = (static_cast<float>(y) + 0.5f) - sampleCoords.y;
            float dist = glm::sqrt(distX * distX + distY * distY);

            // Get index into precomputed filter table based on distance from sample, as opposed to:
            // float filterWeight = Gaussian(scene.info.filterWidth, dist);
            uint8_t filterIndex = glm::min(static_cast<uint8_t>((dist / info.filterWidth) * FilterTableResolution), static_cast<uint8_t>(FilterTableResolution - 1));
            float filterWeight = filterTable[filterIndex];

            // Transform image coords into discrete tile coords
            uint32_t tileX = static_cast<uint32_t>(glm::floor(distX + glm::mod(sampleCoords.x - static_cast<float>(info.filterBounds), static_cast<float>(info.bucketSize)) + static_cast<float>(info.filterBounds)));
            uint32_t tileY = static_cast<uint32_t>(glm::floor(distY + glm::mod(sampleCoords.y - static_cast<float>(info.filterBounds), static_cast<float>(info.bucketSize)) + static_cast<float>(info.filterBounds)));

            uint32_t tileIndex = (tileY * info.tileSize) + tileX;

            pixels[tileIndex].contribution += L * filterWeight;
            pixels[tileIndex].filterWeightSum += filterWeight;
        }
    }
}

struct Camera
{
    Camera(float fov, glm::mat4 cameraToWorld) : fov(fov), cameraToWorld(cameraToWorld) {}

    const float fov = 10.f;
    const glm::mat4 cameraToWorld = glm::mat4(1.f);
};

struct Scene
{
    Scene(RenderInfo info, Camera camera, std::vector<std::shared_ptr<Light>> lights, std::shared_ptr<BVH> bvh) : info(info), camera(camera), lights(lights), bvh(bvh) {}

    RenderInfo info;
    Camera camera;
    std::vector<std::shared_ptr<Light>> lights;
    std::shared_ptr<BVH> bvh;
};


// TODO: Add materials to meshes from JSON
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
                    glm::vec3 Le = glm::vec3(LeGet[0], LeGet[1], LeGet[2]);
                    std::shared_ptr<Light> light = std::make_shared<DistantLight>(Le, lightToWorld);
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

std::vector<Pixel> RenderTile(const Scene& scene, const float* filterTable, uint32_t x0, uint32_t x1, uint32_t y0, uint32_t y1)
{
    std::vector<Pixel> pixels((scene.info.tileSize) * (scene.info.tileSize));

    for (uint32_t y = y0; y < y1; ++y)
    {
        for (uint32_t x = x0; x < x1; ++x)
        {
            // Don't render beyond total image extents when the current tile goes beyond them
            if (x < scene.info.totalWidth && y < scene.info.totalHeight)
            {
                // Create stratified image samples in a Latin square pattern
                std::default_random_engine rng;
                rng.seed(y * scene.info.totalWidth + x);
                std::vector<glm::vec2> imageSamples(scene.info.spp);
                LatinSquare(rng, scene.info.spp, imageSamples);

                for (uint32_t i = 0; i < scene.info.spp; ++i)
                {
                    glm::vec2 imageSample = imageSamples[i];

                    // Create ray in camera-space
                    float aspectRatio = static_cast<float>(scene.info.imageWidth) / static_cast<float>(scene.info.imageHeight);
                    float px = (((static_cast<float>(x) + imageSample.x) / static_cast<float>(scene.info.imageWidth)) * 2.f - 1.f) * glm::tan(glm::radians(scene.camera.fov)) * aspectRatio;
                    float py = (((static_cast<float>(y) + imageSample.y) / static_cast<float>(scene.info.imageHeight)) * -2.f + 1.f) * glm::tan(glm::radians(scene.camera.fov));

                    glm::vec4 o = glm::vec4(0.f, 0.f, 0.f, 1.f);
                    glm::vec4 d(px, py, -1.f, 0.f);
                    glm::normalize(d);

                    // Transform ray into world-space using camera transform
                    glm::mat4 camToWorld = scene.camera.cameraToWorld;

                    o = o * camToWorld;
                    d = d * camToWorld;

                    Ray ray(o, d);

                    // Get sample for L
                    glm::vec4 L(0.f, 0.f, 0.f, 0.f);
                    Intersection isect;
                    if (scene.bvh->Intersect(ray, &isect))
                    {
                        // std::unique_ptr<BxDF> brdf = std::make_unique<LambertBRDF>(isect.sn, glm::vec3(0.8f));
                        BSDF bsdf = isect.material->CreateBSDF(isect.sn);
                        glm::vec3 wo = bsdf.ToLocal(-ray.o);

                        for (const auto& light : scene.lights)
                        {
                            glm::vec3 wi;
                            glm::vec3 Li = light->Li(isect.p, &wi);
                            wi = bsdf.ToLocal(wi);
                            glm::vec3 f = bsdf.f(wo, wi);
                            L += glm::vec4(f * Li * glm::max(wi.z, 0.f), 1.f);
                        }
                        // float Lf = glm::max(glm::dot(glm::normalize(glm::vec3(1, 0.5, 1)), isect.sn), 0.f);
                        // L = glm::vec4(Lf, Lf, Lf, 1.f);
                    }

                    // Transform image sample to "total" image coords (image including filter bounds)
                    glm::vec2 sampleCoords = glm::vec2(static_cast<float>(x + scene.info.filterBounds) + imageSample.x, static_cast<float>(y + scene.info.filterBounds) + imageSample.y);
                    // Add sample to pixels within filter width
                    // Sample coords are respective to total image size
                    AddSample(scene.info, filterTable, sampleCoords, L, pixels);
                }
            }
        }
    }

    return pixels;
}

std::vector<Pixel> Render(const Scene& scene)
{
    std::vector<Pixel> pixels(scene.info.totalHeight * scene.info.totalWidth);

    uint32_t nBucketsX = uint32_t(glm::ceil(static_cast<float>(scene.info.imageWidth) / static_cast<float>(scene.info.bucketSize)));
    uint32_t nBucketsY = uint32_t(glm::ceil(static_cast<float>(scene.info.imageHeight) / static_cast<float>(scene.info.bucketSize)));
    uint32_t nBuckets = nBucketsX * nBucketsY;

    // Pre-compute filter values (only need to do this in 1D as filter is currently only isotropic)
    float filterTable[FilterTableResolution];
    for (uint8_t i = 0; i < FilterTableResolution; ++i)
    {
        filterTable[i] = Gaussian(FilterTableResolution - 1, i);
    }

    // Create tiles and render each one in parallel
    std::vector<Pixel> p(scene.info.tileSize * scene.info.tileSize);
    std::vector<std::vector<Pixel>> tiles(nBuckets, p);

    tbb::task_group tg;

    for (uint32_t y = 0; y < nBucketsY; ++y)
    {
        for (uint32_t x = 0; x < nBucketsX; ++x)
        {
            uint32_t index = y * nBucketsX + x;
            tg.run([&scene, &filterTable, &tiles, index, x, y]
                {
                    std::vector<Pixel> tile = RenderTile(scene, filterTable, scene.info.bucketSize * x, scene.info.bucketSize * (x + 1), scene.info.bucketSize * y, scene.info.bucketSize * (y + 1));
                    tiles[index] = tile;
                });
        }
    }

    tg.wait();

    // Combine tiles into image
    std::cout << "Combining tiles into image...\n";

    for (uint32_t j = 0; j < nBucketsY; ++j)
    {
        for (uint32_t i = 0; i < nBucketsX; ++i)
        {
            for (uint32_t y = 0; y < scene.info.tileSize; ++y)
            {
                for (uint32_t x = 0; x < scene.info.tileSize; ++x)
                {
                    std::vector<Pixel> v = tiles[j * nBucketsX + i];
                    uint32_t pX = x + (i * scene.info.bucketSize);
                    uint32_t pY = y + (j * scene.info.bucketSize);
                    // Ignore tile pixels beyond image bounds
                    if (pX < (scene.info.imageWidth + scene.info.filterBounds) && pY < (scene.info.imageHeight + scene.info.filterBounds))
                    {
                        uint32_t pIndex = (pY * scene.info.totalWidth) + pX;
                        pixels[pIndex].contribution += v[(y * scene.info.tileSize) + x].contribution;
                        pixels[pIndex].filterWeightSum += v[(y * scene.info.tileSize) + x].filterWeightSum;
                    }
                }
            }
        }
    }

    return pixels;
}

void WriteImageToEXR(const RenderInfo& info, const std::vector<Pixel>& pixels, const char* filePath)
{
    // Write image to EXR
    Imf::Array2D<Imf::Rgba> imagePixels(info.imageHeight, info.imageWidth);

    // Ignore pixels beyond image bounds
    for (uint32_t y = info.filterBounds; y < (info.imageHeight + info.filterBounds); ++y)
    {
        for (uint32_t x = info.filterBounds; x < (info.imageWidth + info.filterBounds); ++x)
        {
            glm::vec4 result(0.f);

            result = pixels[y * info.totalWidth + x].contribution / pixels[y * info.totalWidth + x].filterWeightSum;

            imagePixels[y - info.filterBounds][x - info.filterBounds].r = result.r;
            imagePixels[y - info.filterBounds][x - info.filterBounds].g = result.g;
            imagePixels[y - info.filterBounds][x - info.filterBounds].b = result.b;
            imagePixels[y - info.filterBounds][x - info.filterBounds].a = result.a;
        }
    }

    Imf::RgbaOutputFile file(filePath, info.imageWidth, info.imageHeight, Imf::WRITE_RGBA);
    file.setFrameBuffer(&imagePixels[0][0], 1, info.imageWidth);
    file.writePixels(info.imageHeight);
}

int main(int argc, char* argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "Loading " << argv[1] << "...\n";
    Scene scene = LoadScene(argv[1]);

    std::cout << "Rendering...\n";
    std::vector<Pixel> image = Render(scene);

    std::cout << "Writing to " << argv[2] << "...\n";
    WriteImageToEXR(scene.info, image, argv[2]);
        
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = end - start;
    std::cout << "Completed in " << duration.count() << "s\n";

    return 0;
}