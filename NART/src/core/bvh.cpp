#include "bvh.h"

BVH::BoundingVolume::BoundingVolume()
{
    // Bounding volumes are defined by slabs along the cardinal direction vectors
    normals[0] = glm::vec3(1, 0, 0);
    normals[1] = glm::vec3(0, 1, 0);
    normals[2] = glm::vec3(0, 0, 1);

    for (uint8_t i = 0; i < 3; ++i)
    {
        boundsMin[i] = Infinity;
        boundsMax[i] = -Infinity;
    }
}

void BVH::BoundingVolume::ExtendBy(const BoundingVolume& bv)
{
    for (uint8_t i = 0; i < 3; ++i)
    {
        boundsMin[i] = glm::min(boundsMin[i], bv.boundsMin[i]);
        boundsMax[i] = glm::max(boundsMax[i], bv.boundsMax[i]);
    }
}

bool BVH::BoundingVolume::Intersect(const Ray& ray, Intersection* isect)
{
    float tMin = -Infinity;
    float tMax = Infinity;
    float dotLow = Infinity;

    for (uint8_t i = 0; i < 3; ++i)
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



void BVH::Chunk::Insert(const TriMesh& mesh, const uint32_t& index)
{
    triangleIndices.push_back(TriangleIndex(mesh, index));
}

bool BVH::Chunk::Intersect(const Ray& ray, Intersection* isect) const
{
    bool hit = false;
    for (uint32_t i = 0; i < triangleIndices.size(); ++i)
    {
        if (triangleIndices[i].mesh.triangles[triangleIndices[i].index].Intersect(ray, isect))
        {
            hit = true;
            isect->material = triangleIndices[i].mesh.material;
        }
    }
    return hit;
}

void BVH::Chunk::CalculateBounds()
{
    glm::vec3 normals[3];
    normals[0] = glm::vec3(1, 0, 0);
    normals[1] = glm::vec3(0, 1, 0);
    normals[2] = glm::vec3(0, 0, 1);

    for (TriangleIndex triIndex : triangleIndices)
    {
        Triangle triangle = triIndex.mesh.triangles[triIndex.index];

        bboxMin = glm::min(bboxMin, triangle.v0);
        bboxMax = glm::max(bboxMax, triangle.v0);

        bboxMin = glm::min(bboxMin, triangle.v1);
        bboxMax = glm::max(bboxMax, triangle.v1);

        bboxMin = glm::min(bboxMin, triangle.v2);
        bboxMax = glm::max(bboxMax, triangle.v2);

        for (uint8_t i = 0; i < 3; ++i)
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



BVH::Octree::OctreeNode::OctreeNode(glm::vec3 nodeMin, glm::vec3 nodeMax) : nodeMin(nodeMin), nodeMax(nodeMax) {}

BVH::Octree::Octree(glm::vec3 sceneMin, glm::vec3 sceneMax) : sceneMin(sceneMin), sceneMax(sceneMax)
{
    root = std::shared_ptr<OctreeNode>(new OctreeNode(sceneMin, sceneMax));
}

void BVH::Octree::Insert(ChunkPtr chunk)
{
    chunk->CalculateBounds();
    InsertChunkIntoNode(std::move(chunk), root);
}

void BVH::Octree::Build()
{
    BuildBoundingVolumes(root);
}

bool BVH::Octree::Intersect(const Ray& ray, Intersection* isect)
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
                                if (triIsect.tMax < isect->tMax) hit = true;
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

void BVH::Octree::InsertChunkIntoNode(ChunkPtr chunk, std::shared_ptr<OctreeNode> node, uint8_t depth)
{
    if (node->isLeaf)
    {
        if (node->chunks.empty() || depth >= maxDepth)
        {
            node->chunks.push_back(std::move(chunk));
        }

        else
        {
            node->isLeaf = false;

            for (uint32_t i = 0; i < node->chunks.size(); ++i)
            {
                InsertChunkIntoNode(std::move(node->chunks.back()), node, depth++);
                node->chunks.pop_back();
            }

            InsertChunkIntoNode(std::move(chunk), node, depth++);
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

            node->children[nodeIndex] = std::shared_ptr<OctreeNode>(new OctreeNode(childNodeMin, childNodeMax));
        }

        InsertChunkIntoNode(std::move(chunk), node->children[nodeIndex], depth++);
    }
}

void BVH::Octree::BuildBoundingVolumes(std::shared_ptr<OctreeNode> node)
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



BVH::BVH(std::vector<TriMeshPtr>&& _meshes) : meshes(std::move(_meshes))
{
    // Calculate scene extents and number of triangles
    uint32_t numTriangles = 0;
    glm::vec3 sceneMax(-Infinity);
    glm::vec3 sceneMin(Infinity);
    for (uint32_t i = 0; i < meshes.size(); ++i)
    {
        const TriMesh& mesh = *(meshes[i]);
        numTriangles += mesh.triangles.size();
        for (Triangle triangle : mesh.triangles)
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
    std::vector<ChunkPtr> chunks(numChunks + 1);

    // Add triangles to chunks
    for (uint32_t j = 0; j < meshes.size(); ++j)
    {
        for (uint32_t i = 0; i < meshes[j]->triangles.size(); ++i)
        {
            glm::vec3 triangleMin = glm::min(glm::min(meshes[j]->triangles[i].v0, meshes[j]->triangles[i].v1), meshes[j]->triangles[i].v2) - sceneMin;
            glm::vec3 chunkCoords = glm::floor((triangleMin / sceneSize) * resolution);
            uint32_t chunkIndex = static_cast<uint32_t>(glm::floor(chunkCoords.x * resolution.y * resolution.z + chunkCoords.y * resolution.z + chunkCoords.x));
            // TODO: Find source of overshoots
            chunkIndex = glm::min(chunkIndex, numChunks);
            if (!chunks[chunkIndex])
            {
                // std::shared_ptr<Chunk> chunk = std::shared_ptr<Chunk>(new Chunk);
                chunks[chunkIndex] = std::unique_ptr<Chunk>(new Chunk);
            }
            chunks[chunkIndex]->Insert(*(meshes[j]), i);
        }
    }

    // Create octree structure from chunks
    octree = std::unique_ptr<Octree>(new Octree(sceneMin, sceneMax));

    for (uint32_t i = 0; i < chunks.size(); ++i)
    {
        if (chunks[i])
        {
            octree->Insert(std::move(chunks[i]));
        }
    }

    // Create bounding volume for each octree node, from leaves to root
    octree->Build();
}

bool BVH::Intersect(const Ray& ray, Intersection* isect) const
{
    return octree->Intersect(ray, isect);
}


