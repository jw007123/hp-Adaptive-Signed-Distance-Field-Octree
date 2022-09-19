#pragma once

#include <vector>
#include <queue>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "Utility/Literals.h"

#include "Meshing/Utility.h"

namespace Meshing
{
    class NNOctree
    {
    public:
        NNOctree();
        ~NNOctree();

        /// Clears the internal state
        void Clear();

        /// Creates the tree w/ the given domain
        void Create(const Eigen::AlignedBox3f& root_);

        /// Insert a point and stores the u32 alongside the point
        void InsertPoint(const Eigen::Vector3f& pt_, const u32& data_);

        /// Removes a point and its data from the tree
        bool RemovePoint(const Eigen::Vector3f& pt_);

        /// Finds and returns the (approximate) nearest point and its data 
        std::pair<Eigen::Vector3f, u32> NearestPoint(const Eigen::Vector3f& pt_, const f32 maxDistance_) const;

    private:
        typedef std::pair<Eigen::Vector3f, u32> TYPE;

        static constexpr u8  LEAF_NODE_MAX_POINTS = 25;
        static constexpr u32 TYPE_SZ              = sizeof(TYPE);

        struct Node
        {
            Eigen::AlignedBox3f aabb;

            union
            {
                u32   childIdx;
                TYPE* data;
            };

            u8 isLeaf;
            u8 nIdxs;

            Node();
            ~Node();
        };
        std::vector<Node> nodes;

        struct PriorityQueuePredicate
        {
            bool operator()(const std::pair<u32, f32>& v1_, const std::pair<u32, f32>& v2_)
            {
                return v1_.second > v2_.second;
            }
        };
        mutable std::priority_queue<std::pair<u32, f32>, std::vector<std::pair<u32, f32>>, PriorityQueuePredicate> traversalQueue;

        /// Finds the idx and parent idxs of the node that contains pt_ 
        std::pair<u32, u32> FindLeaf(const Eigen::Vector3f& pt_) const;

        /// Subdivides the current node 
        void Subdivide(const u32 nodeIdx_);

        /// Returns the ith child AABB of a node
        Eigen::AlignedBox3f CornerAABB(const Eigen::AlignedBox3f& aabb_, const u8 i_) const;

        /// Performs a NN on a leaf's points, returning the idx of the closest (and second closest if notThis != -1)
        u8 NearestPointInLeaf(const u32 leafIdx_, const Eigen::Vector3f& pt_) const;

        /// Adds the given entry to a leaf's data array
        void PushToLeaf(const u32 leafIdx_, const TYPE& data_);
    };
}