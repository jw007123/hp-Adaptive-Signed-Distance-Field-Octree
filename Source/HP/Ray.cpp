#include "Ray.h"

namespace SDF
{
    Ray::Ray(const Eigen::Vector3d& origin_, const Eigen::Vector3d& direction_)
    {
        origin = origin_;
        direction = direction_;

        invDirection = direction_.cwiseInverse();
        sign.x()     = invDirection.x() < 0.0;
        sign.y()     = invDirection.y() < 0.0;
        sign.z()     = invDirection.z() < 0.0;
    }


    bool Ray::IntersectAABB(const Eigen::AlignedBox3d& aabb_, Eigen::Vector3d& a_, Eigen::Vector3d& b_) const
    {
        // c.f. https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
        const Eigen::Vector3d bounds[2] = 
        {
            aabb_.min(),
            aabb_.max()
        };

        a_.x() = (bounds[sign.x()].x()     - origin.x()) * invDirection.x();
        b_.x() = (bounds[1 - sign.x()].x() - origin.x()) * invDirection.x();
        a_.y() = (bounds[sign.y()].y()     - origin.y()) * invDirection.y();
        b_.y() = (bounds[1 - sign.y()].y() - origin.y()) * invDirection.y();

        if ((a_.x() > b_.y()) || (a_.y() > b_.x()))
        {
            return false;
        }

        if (a_.y() > a_.x())
        {
            a_.x() = a_.y();
        }

        if (b_.y() < b_.x())
        {
            b_.x() = b_.y();
        }

        a_.z() = (bounds[sign.z()].z()     - origin.z()) * invDirection.z();
        b_.z() = (bounds[1 - sign.z()].z() - origin.z()) * invDirection.z();

        if ((a_.x() > b_.z()) || (a_.z() > b_.x()))
        {
            return false;
        }

        if (a_.z() > a_.x())
        {
            a_.x() = a_.z();
        }

        if (b_.z() < b_.x())
        {
            b_.x() = b_.z();
        }

        return true;
    }
}