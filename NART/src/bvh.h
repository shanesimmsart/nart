#pragma once

#include <queue>

#include "geometry.h"

struct BoundingVolume
{
    BoundingVolume();

    // Extend bounding volume by another to ensure it is fully contained
    void ExtendBy(const BoundingVolume& bv);

    bool Intersect(const Ray& ray, Intersection* isect);

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
    void Insert(std::shared_ptr<TriMesh> mesh, const uint32_t& index);

    bool Intersect(const Ray& ray, Intersection* isect) const;

    void CalculateBounds();

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
        OctreeNode(glm::vec3 nodeMin, glm::vec3 nodeMax);

        std::vector<std::shared_ptr<Chunk>> chunks;
        std::shared_ptr<OctreeNode> children[8] = { nullptr };
        bool isLeaf = true;

        const glm::vec3 nodeMin;
        const glm::vec3 nodeMax;

        BoundingVolume bv = BoundingVolume();
    };

public:
    Octree(glm::vec3 sceneMin, glm::vec3 sceneMax);

    void Insert(std::shared_ptr<Chunk> chunk);

    void Build();

    bool Intersect(const Ray& ray, Intersection* isect);

private:
    void InsertChunkIntoNode(std::shared_ptr<Chunk> chunk, std::shared_ptr<OctreeNode> node, uint8_t depth = 0);

    void BuildBoundingVolumes(std::shared_ptr<OctreeNode> node);

    uint8_t maxDepth = 6;
    std::shared_ptr<OctreeNode> root;
    const glm::vec3 sceneMin;
    const glm::vec3 sceneMax;
};

class BVH
{
public:
    BVH(std::vector<std::shared_ptr<TriMesh>> meshes);

    bool Intersect(const Ray& ray, Intersection* isect) const;

    std::vector<std::shared_ptr<Chunk>> chunks;
    std::vector<std::shared_ptr<TriMesh>> meshes;
    std::shared_ptr<Octree> octree;
};


