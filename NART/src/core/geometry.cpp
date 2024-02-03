#include "geometry.h"

Ray::Ray(glm::vec3 o, glm::vec3 d) : o(o), d(d) {}

Triangle::Triangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2) :
    v0(v0), v1(v1), v2(v2), n0(n0), n1(n1), n2(n2) {}

bool Triangle::Intersect(const Ray& ray, Intersection* isect)
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

TriMesh::TriMesh(std::shared_ptr<Material> material) : material(material) {}


