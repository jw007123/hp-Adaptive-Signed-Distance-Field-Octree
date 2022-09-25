#include "Meshing/BVH.h"

namespace Meshing
{
    BVH::Node::Node()
    {
        isLeaf   = true;
        childIdx = -1;
    }


    BVH::Node::~Node()
    {

    }


    BVH::BVH()
    {
        mesh = nullptr;
    }


    BVH::~BVH()
    {

    }

    void BVH::Clear()
    {
        mesh = nullptr;

        nodes.clear();
    }


    void BVH::CreateLeaves(std::vector<ParentInfo>& parents_)
    {
        struct LeafInfo
        {
            Eigen::AlignedBox3f aabb;
            Eigen::Vector3f     centroid;
            bool                available;
        };
        std::vector<LeafInfo> leavesInfo;

        // Create octree
        NNOctree octree;
        octree.Create(nodes[0].aabb);

        // Create initial leaves info and and setup octree
        for (u32 i = 0; i < (mesh->triIndices.size() / 3); ++i)
        {
            // Form triangle
            const Eigen::Vector3f a = mesh->vertices[mesh->triIndices[3 * i + 0]];
            const Eigen::Vector3f b = mesh->vertices[mesh->triIndices[3 * i + 1]];
            const Eigen::Vector3f c = mesh->vertices[mesh->triIndices[3 * i + 2]];

            // Set information and push
            LeafInfo leafInfo;
            leafInfo.aabb      = CalculateAABB(a, b, c);
            leafInfo.centroid  = (a + b + c) * (1.0f / 3.0f);
            leafInfo.available = true;

            leavesInfo.push_back(leafInfo);
            octree.InsertPoint(leafInfo.centroid, i);
        }

        // Create leaves
        u32 lastA = 0;
        while (1)
        {
            // Determine first available aabb
            u32 aIdx = -1;
            for (u32 i = lastA; i < leavesInfo.size(); ++i)
            {
                if (leavesInfo[i].available)
                {
                    aIdx = i;
                    break;
                }
            }
            if (aIdx == -1)
            {
                // Finished
                break;
            }
            lastA = aIdx;

            // Use octree to determine best pairing between leaves
            const LeafInfo& a = leavesInfo[aIdx];
            octree.RemovePoint(a.centroid);

            u32 bIdx = octree.NearestPoint(a.centroid, 10.0f).second;
            if (bIdx == -1)
            {
                // Octree failed due to empty node, need to find manually
                f32 bestDist = std::numeric_limits<f32>::max();
                for (u32 i = lastA; i < leavesInfo.size(); ++i)
                {
                    const f32 dist = leavesInfo[i].available ? (a.centroid - leavesInfo[i].centroid).squaredNorm() : bestDist;
                    if (dist < bestDist && i != aIdx)
                    {
                        bestDist = dist;
                        bIdx     = i;
                    }
                }
            }
            if (bIdx == -1)
            {
                // Special case if we have an odd number of leaves 
                bIdx = aIdx;
            }

            // Obtain pairing!
            const LeafInfo& b = leavesInfo[bIdx];
            octree.RemovePoint(b.centroid);

            // Create leaf nodes out of pair
            Node aLeafNode;
            aLeafNode.aabb   = a.aabb;
            aLeafNode.triIdx = aIdx;
            nodes.push_back(aLeafNode);
            Node bLeafNode;
            bLeafNode.aabb   = b.aabb;
            bLeafNode.triIdx = bIdx;
            nodes.push_back(bLeafNode);

            // Mark nodes as in BVH
            leavesInfo[aIdx].available = false;
            leavesInfo[bIdx].available = false;

            // Add expected parent to array
            ParentInfo parentInfo;
            parentInfo.childIdx   = (u32)nodes.size() - 2;
            parentInfo.parentAABB = (a.aabb).merged(b.aabb);
            parentInfo.available  = true;
            parents_.push_back(parentInfo);
        }
    }


    void BVH::CreateParents(std::vector<ParentInfo>& parents_, std::vector<ParentInfo>& newParents_)
    {
        // Setup octree
        NNOctree octree;
        octree.Create(nodes[0].aabb);
        for (u32 i = 0; i < parents_.size(); ++i)
        {
            octree.InsertPoint(parents_[i].parentAABB.center(), i);
        }

        u32 lastA = 0;
        while (1)
        {
            // Determine first available aabb
            u32 aIdx = -1;
            for (u32 i = lastA; i < parents_.size(); ++i)
            {
                if (parents_[i].available)
                {
                    aIdx = i;
                    break;
                }
            }
            if (aIdx == -1)
            {
                // Finished
                break;
            }
            lastA = aIdx;

            // Use octree to determine best pairing
            const Eigen::AlignedBox3f& a = parents_[aIdx].parentAABB;
            octree.RemovePoint(a.center());

            u32 bIdx = octree.NearestPoint(a.center(), 10.0f).second;
            if (bIdx == -1)
            {
                // Octree failed due to empty node, need to find manually
                f32 bestDist = std::numeric_limits<f32>::max();
                for (u32 i = lastA; i < parents_.size(); ++i)
                {
                    const f32 dist = parents_[i].available ? (a.center() - parents_[i].parentAABB.center()).squaredNorm() : bestDist;
                    if (dist < bestDist && i != aIdx)
                    {
                        bestDist = dist;
                        bIdx     = i;
                    }
                }
            }
            if (bIdx == -1)
            {
                // Special case if we have an odd number of leaves 
                bIdx = aIdx;
            }

            // Obtain pairing!
            const Eigen::AlignedBox3f& b = parents_[bIdx].parentAABB;
            octree.RemovePoint(b.center());

            // Create parents
            Node aParentNode;
            aParentNode.aabb     = a;
            aParentNode.childIdx = parents_[aIdx].childIdx;
            aParentNode.isLeaf   = false;
            nodes.push_back(aParentNode);
            Node bParentNode;
            bParentNode.aabb     = b;
            bParentNode.childIdx = parents_[bIdx].childIdx;
            bParentNode.isLeaf   = false;
            nodes.push_back(bParentNode);

            // Mark nodes as in BVH
            parents_[aIdx].available = false;
            parents_[bIdx].available = false;

            // Add expected parent to array
            ParentInfo parentInfo;
            parentInfo.childIdx   = (u32)nodes.size() - 2;
            parentInfo.parentAABB = a.merged(b);
            parentInfo.available  = true;
            newParents_.push_back(parentInfo);
        }
    }


    bool BVH::Create(const Mesh& mesh_)
    {
        Clear();

        // Set ptrs
        mesh = &mesh_;

        // Create root
        Node rootNode;
        rootNode.aabb = mesh->CalculateMeshAABB();
        nodes.push_back(rootNode);

        // Create leaves and obtain initial parent set
        std::vector<ParentInfo> parentIdxs;
        CreateLeaves(parentIdxs);

        // Create parent levels
        std::vector<ParentInfo> newParentIdxs;
        while (1)
        {
            if (parentIdxs.size() == 1)
            {
                // Done
                nodes[0].isLeaf   = false;
                nodes[0].childIdx = (u32)nodes.size() - 2;
                break;
            }

            // Add parents at this level
            CreateParents(parentIdxs, newParentIdxs);

            // Swap arrays to perform next iteration
            parentIdxs.swap(newParentIdxs);
            newParentIdxs.clear();
        }

        // Setup traversalQueues for querying based on max threads
        for (u32 i = 0; i < std::thread::hardware_concurrency(); ++i)
        {
            traversalQueues.push_back({});
        }
        
        return true;
    }


    ClosestSimplexInfo BVH::ClosestTriangleToPt(const Eigen::Vector3f& pt_, u32& closestTriIdx_, const u32 threadIdx_) const
    {
        if (!nodes.size())
        {
            // Your time is just as important as mine
            assert(0);
        }

        // Queue in .h so that its alloc calls on resizing are only done once
        while (!traversalQueues[threadIdx_].empty())
        {
            traversalQueues[threadIdx_].pop();
        }

        // Set as maximum distance
        ClosestSimplexInfo bestSimplex = {};
        u32                bestTri     = -1;
        f32 bestDistance               = std::numeric_limits<f32>::max();

        // Start at root
        traversalQueues[threadIdx_].push(std::pair<u32, f32>(0, (ClosestPtOnAABB(pt_, nodes[0].aabb) - pt_).squaredNorm()));
        while (!traversalQueues[threadIdx_].empty())
        {
            // Get best untested candidate
            const std::pair<u32, f32> curPair = traversalQueues[threadIdx_].top();
            traversalQueues[threadIdx_].pop();
            
            // Avoid doing any extra work
            if (curPair.second > bestDistance)
            {
                continue;
            }
            
            // Test both children
            for (u8 i = 0; i < 2; ++i)
            {
                // Determine whether any triangles inside AABB can be better than current best
                const u32 childIdx              = nodes[curPair.first].childIdx + i;
                const Node& childNode           = nodes[childIdx];
                const Eigen::Vector3f closestPt = ClosestPtOnAABB(pt_, childNode.aabb);
                const f32 dist                  = (closestPt - pt_).squaredNorm();

                // Nothing inside this branch can be better than our current triangle
                if (dist > bestDistance)
                {
                    continue;
                }

                if (childNode.isLeaf)
                {
                    // Test against triangle in leaf and potentially update bestDistance
                    const u32 t = 3 * childNode.triIdx;
                    const u32 a = mesh->triIndices[t + 0];
                    const u32 b = mesh->triIndices[t + 1];
                    const u32 c = mesh->triIndices[t + 2];

                    const ClosestSimplexInfo simplexInfo = ClosestSimplexToPt(pt_, mesh->vertices[a],
                                                                                    mesh->vertices[b],
                                                                                    mesh->vertices[c]);
                    const f32 dist = (simplexInfo.closestPt - pt_).squaredNorm();

                    // Update best distance w/ new candidate
                    if (dist < bestDistance)
                    {
                        bestDistance = dist;
                        bestSimplex  = simplexInfo;
                        bestTri      = t / 3;
                    }
                }
                else
                {
                    // There's hope yet
                    traversalQueues[threadIdx_].push(std::pair<u32, f32>(childIdx, dist));
                }
            }
        }

        closestTriIdx_ = bestTri;
        return bestSimplex;
    }
}