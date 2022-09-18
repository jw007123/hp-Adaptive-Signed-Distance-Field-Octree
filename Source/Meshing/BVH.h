#pragma once

#include <vector>
#include <queue>
#include <thread>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "Literals.h"

#include "Mesh.h"
#include "Utility.h"
#include "NNOctree.h"

namespace Meshing
{
    class BVH
    {
    public:
        BVH();
        ~BVH();

        /// Clears the internal state
        void Clear();

        /// Creates a BVH from a set of triangle indices and vertices
        bool Create(const Mesh& mesh_);

        /// Returns the closest triangle idx (from triIndices_ in Create) to pt_
        ClosestSimplexInfo ClosestTriangleToPt(const Eigen::Vector3f& pt_, u32& closestTriIdx_, const u32 threadIdx_ = 0) const;

    private:
        struct ParentInfo
        {
            u32                 childIdx;
            Eigen::AlignedBox3f parentAABB;
            bool                available;
        };

        struct Node
        {
            Eigen::AlignedBox3f aabb;

            // Left child = childIdx, right child = childIdx + 1
            union
            {
                u32 childIdx;
                u32 triIdx;
            };

            // Odd num of bytes
            u8 isLeaf;

            Node();
            ~Node();
        };
        std::vector<Node> nodes;

        /*
            Priority queue for traversal. The assumption is that if we follow paths with the closest/best 
            distances first, we can skip more pointless traversal through the tree later on.
        */
        struct PriorityQueuePredicate
        {
            bool operator()(const std::pair<u32, f32>& v1_, const std::pair<u32, f32>& v2_)
            {
                return v1_.second > v2_.second;
            }
        };
        mutable std::vector<std::priority_queue<std::pair<u32, f32>, std::vector<std::pair<u32, f32>>, PriorityQueuePredicate>> traversalQueues;

        const Mesh* mesh;

        /// Creates the lowest level leaves of the BVH
        void CreateLeaves(std::vector<ParentInfo>& parents_);

        /// Creates a parent level 
        void CreateParents(std::vector<ParentInfo>& parents_, std::vector<ParentInfo>& newParents_);
    };
}