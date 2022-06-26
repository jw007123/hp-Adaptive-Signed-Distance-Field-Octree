#include "Simplex.h"

namespace Meshing
{
    ClosestSimplexInfo ClosestSimplexToPt(const Eigen::Vector3f& pt_,
                                          const Eigen::Vector3f& a_, const Eigen::Vector3f& b_, const Eigen::Vector3f& c_)
    {
        // c.f. page 139
        const Eigen::Vector3f ab = b_ - a_;
        const Eigen::Vector3f ac = c_ - a_;
        const Eigen::Vector3f bc = c_ - b_;

        const f32 snom   = (pt_ - a_).dot(ab);
        const f32 sdenom = (pt_ - b_).dot(a_ - b_);
        const f32 tnom   = (pt_ - a_).dot(ac);
        const f32 tdenom = (pt_ - c_).dot(a_ - c_);

        if (snom < EPSILON_F32 && tnom < EPSILON_F32)
        {
            ClosestSimplexInfo simpInfo;
            simpInfo.closestPt  = a_;
            simpInfo.simplex    = Simplex::Vertex;
            simpInfo.simplexIdx = SimplexID::VertexA;

            return simpInfo;
        }

        const f32 unom   = (pt_ - b_).dot(bc);
        const f32 udenom = (pt_ - c_).dot(b_ - c_);

        if (sdenom < EPSILON_F32 && unom < EPSILON_F32)
        {
            ClosestSimplexInfo simpInfo;
            simpInfo.closestPt  = b_;
            simpInfo.simplex    = Simplex::Vertex;
            simpInfo.simplexIdx = SimplexID::VertexB;

            return simpInfo;
        }

        if (tdenom < EPSILON_F32 && udenom < EPSILON_F32)
        {
            ClosestSimplexInfo simpInfo;
            simpInfo.closestPt  = c_;
            simpInfo.simplex    = Simplex::Vertex;
            simpInfo.simplexIdx = SimplexID::VertexC;

            return simpInfo;
        }

        const Eigen::Vector3f n = (b_ - a_).cross(c_ - a_);
        const f32 vc            = n.dot((a_ - pt_).cross(b_ - pt_));

        if (vc < EPSILON_F32 && snom > EPSILON_F32 && sdenom > EPSILON_F32)
        {
            ClosestSimplexInfo simpInfo;
            simpInfo.closestPt  = a_ + snom / (snom + sdenom) * ab;
            simpInfo.simplex    = Simplex::Edge;
            simpInfo.simplexIdx = SimplexID::EdgeAB;

            return simpInfo;
        }

        const f32 va = n.dot((b_ - pt_).cross(c_ - pt_));

        if (va < EPSILON_F32 && unom > EPSILON_F32 && udenom > EPSILON_F32)
        {
            ClosestSimplexInfo simpInfo;
            simpInfo.closestPt  = b_ + unom / (unom + udenom) * bc;
            simpInfo.simplex    = Simplex::Edge;
            simpInfo.simplexIdx = SimplexID::EdgeBC;

            return simpInfo;
        }

        const f32 vb = n.dot((c_ - pt_).cross(a_ - pt_));

        if (vb < EPSILON_F32 && tnom > EPSILON_F32 && tdenom > EPSILON_F32)
        {
            ClosestSimplexInfo simpInfo;
            simpInfo.closestPt  = a_ + tnom / (tnom + tdenom) * ac;
            simpInfo.simplex    = Simplex::Edge;
            simpInfo.simplexIdx = SimplexID::EdgeCA;

            return simpInfo;
        }

        const f32 u = va / (va + vb + vc);
        const f32 v = vb / (va + vb + vc);
        const f32 w = 1.0 - u - v;

        ClosestSimplexInfo simpInfo;
        simpInfo.closestPt  = u * a_ + v * b_ + w * c_;
        simpInfo.simplex    = Simplex::Face;
        simpInfo.simplexIdx = SimplexID::FaceABC;

        return simpInfo;
    }


    Eigen::AlignedBox3f CalculateAABB(const Eigen::Vector3f& a_, const Eigen::Vector3f& b_, const Eigen::Vector3f& c_)
    {
        Eigen::AlignedBox3f ret;
        ret.setEmpty();

        for (u8 i = 0; i < 3; ++i)
        {
            ret.min()(i) = std::min<f32>(a_(i), std::min<f32>(b_(i), c_(i)));
            ret.max()(i) = std::max<f32>(a_(i), std::max<f32>(b_(i), c_(i)));
        }

        return ret;
    }


    Eigen::Vector3f ClosestPtOnAABB(const Eigen::Vector3f& p_, const Eigen::AlignedBox3f& aabb_)
    {
        Eigen::Vector3f ret;

        for (u8 i = 0; i < 3; ++i)
        {
            f32 v = p_(i);
            const f32 min = aabb_.min()(i);
            const f32 max = aabb_.max()(i);

            if (v < min)
            {
                v = min;
            }

            if (v > max)
            {
                v = max;
            }

            ret(i) = v;
        }

        return ret;
    }
}