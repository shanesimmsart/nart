#pragma once

#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <vector>
#include <memory>

#include "nart.h"
#include "materials.h"

#define Infinity std::numeric_limits<float>::infinity()

class Ray
{
public:
    Ray(glm::vec3 o, glm::vec3 d) : o(o), d(d) {}

    glm::vec3 o = glm::vec3(0.f, 0.f, 0.f);
    glm::vec3 d = glm::vec3(0.f, 1.f, 0.f);
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

    // Flag for area light intersection
    bool isLight = false;
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

        // Back-face culling
        if (glm::dot(n, ray.d) >= 0.f) return false;

        // Find t of intersection with plane
        // Note: n doesn't need to be normalised yet as denominator divides by n's length
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


