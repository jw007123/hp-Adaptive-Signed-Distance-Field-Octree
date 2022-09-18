#pragma once

#include "Literals.h"

#include "Eigen/Core"
#include "Eigen/Geometry"

namespace Meshing
{
    enum Simplex : u8
    {
        Vertex = 0,
        Edge   = 1,
        Face   = 2
    };

    enum SimplexID : u8
    {
        VertexA = 0,
        VertexB = 1,
        VertexC = 2,
        EdgeAB  = 0,
        EdgeBC  = 1,
        EdgeCA  = 2,
        FaceABC = 0
    };

    /// http://www.r-5.org/files/books/computers/algo-list/realtime-3f/Christer_Ericson-Real-Time_Collision_Detection-EN.pdf
    struct ClosestSimplexInfo
    {
        Simplex         simplex;
        SimplexID       simplexIdx;
        Eigen::Vector3f closestPt;
    };
    ClosestSimplexInfo ClosestSimplexToPt(const Eigen::Vector3f& pt_, 
                                          const Eigen::Vector3f& a_, const Eigen::Vector3f& b_, const Eigen::Vector3f& c_);

    /// Calculates and returns the bounding AABB of the triangle
    Eigen::AlignedBox3f CalculateAABB(const Eigen::Vector3f& a_, const Eigen::Vector3f& b_, const Eigen::Vector3f& c_);

    /// Returns the closest point in the aabb to p_
    Eigen::Vector3f ClosestPtOnAABB(const Eigen::Vector3f& p_, const Eigen::AlignedBox3f& aabb_);
}