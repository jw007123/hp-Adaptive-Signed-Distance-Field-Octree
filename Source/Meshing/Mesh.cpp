#include "Mesh.h"

namespace Meshing
{
    Mesh::Mesh()
    {
        hasValidState = false;
    }


    Mesh::~Mesh()
    {

    }


    void Mesh::Clear()
    {
        hasValidState = false;

        halfEdges.clear();
        triIndices.clear();
        vertexNormals.clear();
        vertices.clear();
    }


    void Mesh::CreateFromObj(const char* objPath_)
    {
        Clear();

        ObjParser objParser;
        objParser.Load(objPath_);

        // Only need to fill in halfEdges array. Rest can be taken straight away
        const std::vector<Eigen::Vector3f>& objVerts = objParser.GetVertices();
        vertices.resize(objVerts.size());
        memcpy(vertices.data(), objVerts.data(), sizeof(Eigen::Vector3f) * objVerts.size());

        const std::vector<Eigen::Vector3f>& objVertNorms = objParser.GetVertexNormals();
        vertexNormals.resize(objVertNorms.size());
        memcpy(vertexNormals.data(), objVertNorms.data(), sizeof(Eigen::Vector3f) * objVertNorms.size());

        const std::vector<u32>& objTriIndices = objParser.GetTriIndices();
        triIndices.resize(objTriIndices.size());
        memcpy(triIndices.data(), objTriIndices.data(), sizeof(u32) * objTriIndices.size());

        hasValidState = CreateHalfEdges();
    }


    f32 Mesh::SignedDistanceAtPt(const Eigen::Vector3f& pt_)
    {
        if (!hasValidState)
        {
            return std::numeric_limits<f32>::max();
        }

        u32 closestTriIdx                    = -1;
        const ClosestSimplexInfo simplexInfo = ClosestTriangleToPt(pt_, closestTriIdx);
        const Eigen::Vector3f pseudoNormal   = PseudoNormal(closestTriIdx, simplexInfo.simplex, simplexInfo.simplexIdx);

        // Use info to calculate final signed distance
        const f32 sign = pseudoNormal.dot(pt_ - simplexInfo.closestPt) > 0.0 ? 1.0 : -1.0;
        return sign * (pt_ - simplexInfo.closestPt).norm();
    }


    f32 Mesh::SignedDistanceAtPt(const Eigen::Vector3f& pt_, const BVH& bvh_, const u32 threadIdx_)
    {
        if (!hasValidState)
        {
            return std::numeric_limits<f32>::max();
        }

        u32 closestTriIdx                    = -1;
        const ClosestSimplexInfo simplexInfo = bvh_.ClosestTriangleToPt(pt_, closestTriIdx, threadIdx_);
        const Eigen::Vector3f pseudoNormal   = PseudoNormal(closestTriIdx, simplexInfo.simplex, simplexInfo.simplexIdx);

        // Use info to calculate final signed distance
        const f32 sign = pseudoNormal.dot(pt_ - simplexInfo.closestPt) > 0.0 ? 1.0 : -1.0;
        return sign * (pt_ - simplexInfo.closestPt).norm();
    }


    Eigen::AlignedBox3f Mesh::CalculateMeshAABB() const
    {
        Eigen::AlignedBox3f ret;
        ret.setEmpty();

        for (u32 i = 0; (i < triIndices.size() / 3); ++i)
        {
            // Form triangle
            const Eigen::Vector3f& a = vertices[triIndices[3 * i + 0]];
            const Eigen::Vector3f& b = vertices[triIndices[3 * i + 1]];
            const Eigen::Vector3f& c = vertices[triIndices[3 * i + 2]];

            // Obtain AABB from triangle
            const Eigen::AlignedBox3f triAABB = CalculateAABB(a, b, c);
            ret.extend(triAABB);
        }

        return ret;
    }


    bool Mesh::CreateHalfEdges()
    {
        std::map<std::pair<u32, u32>, u32> edgeMap;
        halfEdges.assign(triIndices.size(), -1);

        for (u32 i = 0; i < triIndices.size(); ++i)
        {
            // Create edge pair
            std::pair<u32, u32> edge;
            std::pair<u32, u32> edgeRev;
            if (i % 3 == 2)
            {
                edge = std::pair<u32, u32>(triIndices[i], triIndices[i - 2]);
            }
            else
            {
                edge = std::pair<u32, u32>(triIndices[i], triIndices[i + 1]);
            }
            edgeRev = std::pair<u32, u32>(edge.second, edge.first);

            // Use map to look for edgeRev
            const auto& edgeFind = edgeMap.find(edgeRev);
            if (edgeFind != edgeMap.end())
            {
                halfEdges[(*edgeFind).second] = i;
                halfEdges[i]                  = (*edgeFind).second;
            }
            else
            {
                // Insert edge and its position in triIndices into map
                edgeMap.insert(std::pair<std::pair<u32, u32>, u32>(edge, i));
            }
        }

        // Check array has no -1 values. Otherwise unable to do SDF on mesh
        for (u32 i = 0; i < halfEdges.size(); ++i)
        {
            if (halfEdges[i] == -1)
            {
                return false;
            }
        }

        return true;
    }


    ClosestSimplexInfo Mesh::ClosestTriangleToPt(const Eigen::Vector3f& pt_, u32& closestTriIdx_) const
    {
        ClosestSimplexInfo bestSimplex = {};
        u32                bestTri     = -1;
        f32                bestDist    = std::numeric_limits<f32>::max();

        for (u32 i = 0; i < (triIndices.size() / 3); ++i)
        {
            // Determine dist to closest pt on triangle
            const ClosestSimplexInfo cInfo = ClosestSimplexToPt(pt_, 
                                             vertices[triIndices[i * 3 + 0]], 
                                             vertices[triIndices[i * 3 + 1]], 
                                             vertices[triIndices[i * 3 + 2]]);
            const f32 thisDist             = (pt_ - cInfo.closestPt).squaredNorm();

            if (thisDist < bestDist)
            {
                bestSimplex = cInfo;
                bestTri     = i;
                bestDist    = thisDist;
            }
        }

        closestTriIdx_ = bestTri;
        return bestSimplex;
    }


    Eigen::Vector3f Mesh::PseudoNormal(const u32 triIndex_, const Simplex simplex_, const SimplexID simplexIdx_)
    {
        switch (simplex_)
        {
            case Simplex::Vertex:
            {
                return PseudoNormalVertex(triIndex_, simplexIdx_);
            }

            case Simplex::Edge:
            {
                return PseudoNormalEdge(triIndex_, simplexIdx_);
            }

            case Simplex::Face:
            {
                return PseudoNormalFace(triIndex_);
            }

            default:
                assert(0);
        }
    }


    Eigen::Vector3f Mesh::PseudoNormalFace(const u32 triIndex_)
    {
        const Eigen::Vector3f ab = vertices[triIndices[triIndex_ * 3 + 1]] - vertices[triIndices[triIndex_ * 3 + 0]];
        const Eigen::Vector3f ac = vertices[triIndices[triIndex_ * 3 + 2]] - vertices[triIndices[triIndex_ * 3 + 0]];
        const Eigen::Vector3f n  = ab.cross(ac);

        // Regular face normal
        return n.normalized();
    }


    Eigen::Vector3f Mesh::PseudoNormalEdge(const u32 triIndex_, const SimplexID simplexIdx_)
    {
        // Determine adjacent triangle
        const u32 adjEdgeIdx = halfEdges[triIndex_ * 3 + simplexIdx_];
        const u32 adjTriIdx  = (adjEdgeIdx - (adjEdgeIdx % 3)) / 3;

        // Determine weighted average of both normals
        const Eigen::Vector3f nA = PseudoNormalFace(triIndex_);
        const Eigen::Vector3f nB = PseudoNormalFace(adjTriIdx);
        const Eigen::Vector3f n  = nA * PI + nB * PI;

        return n.normalized();
    }


    Eigen::Vector3f Mesh::PseudoNormalVertex(const u32 triIndex_, const SimplexID simplexIdx_)
    {
        Eigen::Vector3f n = Eigen::Vector3f::Zero();
        u32 curHEIdx      = triIndex_ * 3 + simplexIdx_;
        u32 curTriIdx     = triIndex_;

        do
        {
            Eigen::Vector3f curTri[3];
            curTri[0] = vertices[triIndices[curTriIdx * 3]];
            curTri[1] = vertices[triIndices[curTriIdx * 3 + 1]];
            curTri[2] = vertices[triIndices[curTriIdx * 3 + 2]];

            const u32 hemod3  = curHEIdx % 3;
            const u32 he1mod3 = (curHEIdx + 1) % 3;
            const u32 he2mod3 = (curHEIdx + 2) % 3;

            const Eigen::Vector3f ab = curTri[he1mod3] - curTri[hemod3];
            const Eigen::Vector3f ac = curTri[he2mod3] - curTri[hemod3];
            const f32 ang            = acosf(ab.normalized().dot(ac.normalized()));

            n += PseudoNormalFace(curTriIdx) * ang;

            curHEIdx = halfEdges[curHEIdx];
            curHEIdx = ((curHEIdx % 3) == 2) ? (curHEIdx - 2) : (curHEIdx + 1);

            curTriIdx = (curHEIdx - (curHEIdx % 3)) / 3;

        } while (curTriIdx != triIndex_);

        return n.normalized();
    }
}