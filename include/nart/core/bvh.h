#pragma once

#include <queue>

#include "geometry.h"

class BVH {
 public:
  BVH(std::vector<TriMeshPtr>&& _meshes);

  bool Intersect(const Ray& ray, Intersection& isect) const;

 private:
  struct BoundingVolume {
    BoundingVolume();

    // Extend bounding volume by another to ensure it is fully contained
    void ExtendBy(const BoundingVolume& bv);

    bool Intersect(const Ray& ray, Intersection& isect) const;

    glm::vec3 normals[3];
    float boundsMin[3];
    float boundsMax[3];
  };

  // A chunk of triangles to be placed inside of a bounding volume
  struct Chunk {
   public:
    void Insert(const TriMesh& mesh, const uint32_t& index);

    bool Intersect(const Ray& ray, Intersection& isect) const;

    void CalculateBounds();

    // Bounding box (used for octree position)
    glm::vec3 bboxMin = glm::vec3(Infinity);
    glm::vec3 bboxMax = glm::vec3(-Infinity);
    // Bounding volume (used for intersecting)
    BoundingVolume bv = BoundingVolume();

   private:
    // Stores a mesh pointer and an index into one of its triangles
    struct TriangleIndex {
      TriangleIndex(const TriMesh& mesh, const uint32_t& index)
          : mesh(mesh), index(index) {}
      const TriMesh& mesh;
      const uint32_t index;
    };

    std::vector<TriangleIndex> triangleIndices;
  };

  using ChunkPtr = std::unique_ptr<Chunk>;

  // BVH uses an octree structure for containing bounding volumes
  struct Octree {
   private:
    // Forward-declaration
    struct OctreeNode;

    using OctreeNodePtr = std::unique_ptr<OctreeNode>;

   public:
    Octree(const glm::vec3& sceneMin, const glm::vec3& sceneMax);

    void Insert(ChunkPtr chunk);

    void Build();

    bool Intersect(const Ray& ray, Intersection& isect);

    OctreeNodePtr root;
    const glm::vec3 sceneMin;
    const glm::vec3 sceneMax;

   private:
    struct OctreeNode {
     public:
      OctreeNode(const glm::vec3& nodeMin, const glm::vec3& nodeMax);

      std::vector<ChunkPtr> chunks;
      std::unique_ptr<OctreeNode> children[8] = {nullptr};
      bool isLeaf = true;

      const glm::vec3 nodeMin;
      const glm::vec3 nodeMax;

      BoundingVolume bv = BoundingVolume();

      void InsertChunk(ChunkPtr chunk, uint8_t depth = 0);

      void BuildBoundingVolumes();

     private:
      const uint8_t maxDepth = 5;
    };

    using OctreeNodePtr = std::unique_ptr<OctreeNode>;
  };

  using OctreePtr = std::unique_ptr<Octree>;

  std::vector<TriMeshPtr> meshes;
  OctreePtr octree;
};

using BVHPtr = std::unique_ptr<BVH>;
