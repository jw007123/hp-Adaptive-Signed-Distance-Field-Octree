#pragma once

#include "Eigen/Core"
#include "Eigen/Geometry"

namespace SDF
{
    /// Defines a line R(t) = O + t * D, where ||D||_2 = 1
    struct Ray
    {
        Ray(const Eigen::Vector3d& origin_, const Eigen::Vector3d& direction_);

        Eigen::Vector3d origin;
        Eigen::Vector3d direction;
        
        Eigen::Vector3d invDirection;
        Eigen::Vector3i sign;

        /// Returns whether the ray intersects the box and, if so, stores the intersection points in a_ and b_
        bool IntersectAABB(const Eigen::AlignedBox3d& aabb_, Eigen::Vector3d& a_, Eigen::Vector3d& b_) const;
    };
}