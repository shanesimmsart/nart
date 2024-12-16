#pragma once

#include <queue>

#include "geometry.h"

class BVH
{
public:
    BVH(std::vector<TriMeshPtr>&& _meshes);

    bool Intersect(const Ray& ray, Intersection* isect) const;

private:
    struct BoundingVolume
    {
    public:
        // Extend bounding volume by another to ensure it is fully contained
        void ExtendBy(const BoundingVolume& bv);

        bool Intersect(const Ray& ray, Intersection* isect) const;

        glm::vec3 normals[3];
        float boundsMin[3];
        float boundsMax[3];

    private:
        BoundingVolume();

        friend class BVH;
    };

    // A chunk of triangles to be placed inside of a bounding volume
    class Chunk
    {
        // Stores a mesh pointer and an index into one of its triangles
        struct TriangleIndex
        {
        public:
            const TriMesh& mesh;
            const uint32_t index;

        private:
            TriangleIndex(const TriMesh& mesh, const uint32_t& index) : mesh(mesh), index(index) {}

            friend class Chunk;
        };

    public:
        void Insert(const TriMesh& mesh, const uint32_t& index);

        bool Intersect(const Ray& ray, Intersection* isect) const;

        void CalculateBounds();

        // Bounding box (used for octree position)
        glm::vec3 bboxMin = glm::vec3(Infinity);
        glm::vec3 bboxMax = glm::vec3(-Infinity);
        // Bounding volume (used for intersecting)
        BoundingVolume bv = BoundingVolume();

    private:
        Chunk() = default;

        std::vector<TriangleIndex> triangleIndices;

        friend class BVH;
    };

    using ChunkPtr = std::unique_ptr<Chunk>;
    
    // BVH uses an octree structure for containing bounding volumes
    class Octree
    {
        struct OctreeNode
        {
        public:
            std::vector<ChunkPtr> chunks;
            std::unique_ptr<OctreeNode> children[8] = { nullptr };
            bool isLeaf = true;

            const glm::vec3 nodeMin;
            const glm::vec3 nodeMax;

            BoundingVolume bv = BoundingVolume();

            void InsertChunk(ChunkPtr chunk, uint8_t depth = 0);

            void BuildBoundingVolumes();

        private:
            const uint8_t maxDepth = 5;

            OctreeNode(glm::vec3 nodeMin, glm::vec3 nodeMax);

            friend class Octree;
        };

        using OctreeNodePtr = std::unique_ptr<OctreeNode>;

    public:
        void Insert(ChunkPtr chunk);

        void Build();

        bool Intersect(const Ray& ray, Intersection* isect);

    private:
        Octree(glm::vec3 sceneMin, glm::vec3 sceneMax);

        friend class BVH;

        OctreeNodePtr root;
        const glm::vec3 sceneMin;
        const glm::vec3 sceneMax;
    };

    using OctreePtr = std::unique_ptr<Octree>;

    std::vector<TriMeshPtr> meshes;
    OctreePtr octree;
};

using BVHPtr = std::unique_ptr<BVH>;


