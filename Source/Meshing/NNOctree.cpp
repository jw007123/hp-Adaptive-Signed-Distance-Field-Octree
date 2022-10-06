#include "Meshing/NNOctree.h"

namespace Meshing
{
    NNOctree::Node::Node()
    {
        isLeaf = true;
        data   = nullptr;
        nIdxs  = 0;
    }


    NNOctree::Node::~Node()
    {
        data  = nullptr;
        nIdxs = 0;
    }


    NNOctree::~NNOctree()
    {
        for (u32 i = 0; i < nodes.size(); ++i)
        {
            if (nodes[i].isLeaf && nodes[i].data != nullptr)
            {
                free(nodes[i].data);
            }
        }
    }

    
    void NNOctree::Clear()
    {
        for (u32 i = 0; i < nodes.size(); ++i)
        {
            if (nodes[i].isLeaf && nodes[i].data != nullptr)
            {
                free(nodes[i].data);
            }
        }
        nodes.clear();
    }


    void NNOctree::Create(const Eigen::AlignedBox3f& root_)
    {
        Clear();

        // Create root and subdivide to avoid special casing elsewhere
        Node root;
        root.aabb = root_;
        nodes.push_back(root);
        Subdivide(0);
    }

    
    void NNOctree::InsertPoint(const Eigen::Vector3f& pt_, const u32& data_)
    {
        if (!nodes.size())
        {
            assert(0);
        }

        // Find leaf that can (potentially) hold pt_
        u32 leafIdx = FindLeaf(pt_).first;
        if (leafIdx == -1)
        {
            return;
        }

        // If data array full, subdivide
        if (nodes[leafIdx].nIdxs >= LEAF_NODE_MAX_POINTS)
        {
            Subdivide(leafIdx);
            leafIdx  = FindLeaf(pt_).first;
        }
        PushToLeaf(leafIdx, TYPE(pt_, data_));
    }


    bool NNOctree::RemovePoint(const Eigen::Vector3f& pt_)
    {
        if (!nodes.size())
        {
            assert(0);
        }

        // Obtain leaf and check pt_ is within domain
        const std::pair<u32, u32> idxs = FindLeaf(pt_);
        if (idxs.first == -1)
        {
            return false;
        }

        // Shuffle points in leaf and resize if leaf has data
        Node& leafNode     = nodes[idxs.first];
        const u8 removeIdx = NearestPointInLeaf(idxs.first, pt_);
        if (removeIdx == (u8)-1)
        {
            return false;
        }

        for (u8 i = removeIdx; i < (leafNode.nIdxs - 1); ++i)
        {
            leafNode.data[i] = leafNode.data[i + 1];
        }

        // Realloc
        leafNode.nIdxs -= 1;
        if (leafNode.nIdxs)
        {
            leafNode.data = (TYPE*)realloc(leafNode.data, TYPE_SZ * leafNode.nIdxs);
        }
        else
        {
            free(leafNode.data);
            leafNode.data = nullptr;
        }

        return true;
    }

 
    std::pair<Eigen::Vector3f, u32> NNOctree::NearestPoint(const Eigen::Vector3f& pt_, const f32 maxDistance_) const
    {
        if (!nodes.size())
        {
            assert(0);
        }

        // Empty queue
        while (!traversalQueue.empty())
        {
            traversalQueue.pop();
        }

        // Start at root
        traversalQueue.push(std::pair<u32, f32>(0, (ClosestPtOnAABB(pt_, nodes[0].aabb) - pt_).squaredNorm()));

        // Search tree
        TYPE bestItem = TYPE(pt_, -1);
        f32  bestDist = maxDistance_;
        while (!traversalQueue.empty())
        {
            const std::pair<u32, f32> curItem = traversalQueue.top();
            traversalQueue.pop();

            // Pointless work
            if (curItem.second > bestDist)
            {
                continue;
            }

            if (nodes[curItem.first].isLeaf)
            {
                // Get best point in node and update bests
                const u8 bestPtIdx = NearestPointInLeaf(curItem.first, pt_);
                if (bestPtIdx != (u8)-1)
                {
                    const f32 dist = (nodes[curItem.first].data[bestPtIdx].first - pt_).squaredNorm();
                    if (dist < bestDist)
                    {
                        bestDist = dist;
                        bestItem = nodes[curItem.first].data[bestPtIdx];
                    }
                }
            }
            else
            {
                // Add compatible child nodes to queue
                for (u8 i = 0; i < 8; ++i)
                {
                    const u32 childIdx    = nodes[curItem.first].childIdx + i;
                    const Node& childNode = nodes[childIdx];
                    const f32 dist        = (ClosestPtOnAABB(pt_, childNode.aabb) - pt_).squaredNorm();

                    if (dist < bestDist)
                    {
                        traversalQueue.push(std::pair<u32, f32>(childIdx, dist));
                    }
                }
            }
        }

        return bestItem;
    }


    std::pair<u32, u32> NNOctree::FindLeaf(const Eigen::Vector3f& pt_) const
    {
        // Not in volume
        if (!nodes[0].aabb.contains(pt_))
        {
            return std::pair<u32, u32>(-1, -1);
        }

        // Start at root
        u32 curNodeIdx = 0;
        while (1)
        {
            const Eigen::Vector3f& aabbMin = nodes[curNodeIdx].aabb.min();
            const Eigen::Vector3f& aabbMax = nodes[curNodeIdx].aabb.max();
            const f32 curNodeAABBHalfX     = (aabbMax.x() - aabbMin.x()) * 0.5f;
            const f32 curNodeAABBHalfY     = (aabbMax.y() - aabbMin.y()) * 0.5f;
            const f32 curNodeAABBHalfZ     = (aabbMax.z() - aabbMin.z()) * 0.5f;

            const u32 xIdx = (pt_.x() >= (aabbMin.x() + curNodeAABBHalfX));
            const u32 yIdx = (pt_.y() >= (aabbMin.y() + curNodeAABBHalfY)) << 1;
            const u32 zIdx = (pt_.z() >= (aabbMin.z() + curNodeAABBHalfZ)) << 2;

            const u32 childIdx   = (u32)(nodes[curNodeIdx].childIdx + xIdx + yIdx + zIdx);
            const Node& curChild = nodes[childIdx];
            if (curChild.isLeaf)
            {
                return std::pair<u32, u32>(childIdx, curNodeIdx);
            }
            else
            {
                // Not leaf
                curNodeIdx = childIdx;
            }
        }
    }


    void NNOctree::Subdivide(const u32 nodeIdx_)
    {
        TYPE* dataPtr = nodes[nodeIdx_].data;

        // Create children
        for (u8 i = 0; i < 8; ++i)
        {
            nodes.push_back({});

            Node& newNode = nodes.back();
            newNode.aabb  = CornerAABB(nodes[nodeIdx_].aabb, i);
        }

        // Make non-leaf. Defer data free until points re-added
        nodes[nodeIdx_].childIdx = (u32)nodes.size() - 8;
        nodes[nodeIdx_].isLeaf   = false;

        // Insert parent points in new leaves
        for (u8 i = 0; i < nodes[nodeIdx_].nIdxs; ++i)
        {
            const u32 leafIdx = FindLeaf(dataPtr[i].first).first;
            PushToLeaf(leafIdx, dataPtr[i]);
        }
        
        // Free mem
        if (dataPtr)
        {
            free(dataPtr);
        }
        nodes[nodeIdx_].nIdxs = 0;
    }


    Eigen::AlignedBox3f NNOctree::CornerAABB(const Eigen::AlignedBox3f& aabb_, const u8 i_) const
    {
        Eigen::AlignedBox3f cornerAABB = aabb_;
        for (u8 d = 0; d < 3; ++d)
        {
            if (i_ & (1 << d))
            {
                cornerAABB.min().coeffRef(d) = (aabb_.max().coeff(d) + aabb_.min().coeff(d)) * 0.5f;
            }
            else
            {
                cornerAABB.max().coeffRef(d) = (aabb_.max().coeff(d) + aabb_.min().coeff(d)) * 0.5f;
            }
        }

        return cornerAABB;
    }


    u8 NNOctree::NearestPointInLeaf(const u32 leafIdx_, const Eigen::Vector3f& pt_) const
    {
        const Node& leafNode = nodes[leafIdx_];
        u8 bestIdx           = -1;
        f32 bestDist         = std::numeric_limits<f32>::max();

        for (u8 i = 0; i < leafNode.nIdxs; ++i)
        {
            const f32 dist = (leafNode.data[i].first - pt_).squaredNorm();
            if (dist < bestDist)
            {
                bestIdx  = i;
                bestDist = dist;
            }
        }

        return bestIdx;
    }


    void NNOctree::PushToLeaf(const u32 leafIdx_, const TYPE& data_)
    {
        Node& leafNode = nodes[leafIdx_];

        if (leafNode.data != nullptr)
        {
            leafNode.data = (TYPE*)realloc(leafNode.data, TYPE_SZ * (leafNode.nIdxs + 1));
        }
        else
        {
            leafNode.data = (TYPE*)malloc(TYPE_SZ);
        }

        leafNode.data[leafNode.nIdxs] = data_;
        leafNode.nIdxs++;
    }
}
