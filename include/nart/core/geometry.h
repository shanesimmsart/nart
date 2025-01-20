#pragma once

#include <glm/glm.hpp>
#include <glm/vec3.hpp>
#include <memory>
#include <vector>

#include "material.h"

#define Infinity std::numeric_limits<float>::infinity()

struct Ray {
    Ray(const glm::vec3& o, const glm::vec3& d);

    glm::vec3 o = glm::vec3(0.f, 0.f, 0.f);
    glm::vec3 d = glm::vec3(0.f, 0.f, 1.f);

    uint8_t majorAxis = 2;
    float Sx = 1.f;
    float Sy = 1.f;
    float Sz = 1.f;
};

struct Intersection {
    Material* material;
    // Barycentric coords
    float u, v;
    // Point, geometric normal, shading normal
    glm::vec3 p, gn, sn;
    // UVs
    glm::vec2 st;

    // Upper and lower bounds of ray to be considered
    float tMin = 0.f;
    float tMax = Infinity;

    // Flag for area light intersection
    bool isLight = false;
};

struct Triangle {
    Triangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2,
             const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2,
             const glm::vec2& uv0 = glm::vec2(0.f),
             const glm::vec2& uv1 = glm::vec2(0.f),
             const glm::vec2& uv2 = glm::vec2(0.f));

    bool Intersect(const Ray& ray, Intersection& isect) const;

    const glm::vec3  v0,  v1,  v2;
    const glm::vec3  n0,  n1,  n2;
    const glm::vec2 uv0, uv1, uv2;
};

class TriMesh {
public:
    TriMesh(MaterialPtr&& material);

    // Passing as pointer instead of const ref as Intersect needs to hold the
    // pointer and Material will call CreateBSDF() (non-const)
    Material* GetMaterial() const;

    // Triangles are passed directly, outside of the ctor, to reduce memory
    // allocations
    std::vector<Triangle> triangles;
    MaterialPtr material;
};

using TriMeshPtr = std::unique_ptr<TriMesh>;
