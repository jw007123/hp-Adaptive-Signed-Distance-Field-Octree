#pragma once

#include <vector>
#include <tuple>
#include <map>
#include <cassert>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "Utility/Literals.h"
#include "Meshing/ObjParser.h"
#include "Meshing/Utility.h"

namespace Meshing
{
    // Forward declare so BVH can access mesh vertices etc
    class BVH;

    /*
        Mesh class using a half-edge structure. The halfEdges array can be interpreted as follows:
            
            triIndices[i] denotes an edge starting at points[triIndices[i]] to points[triIndices[i + 1]], ignoring 
            mod_3 special-casing. 

            halfEdges[i] can be used to find the edge of any triangle that shares the edge starting at triIndices[i]. If
            there is no adjacent edge, then halfEdges[i] = -1.

            0.......................1       For example, consider the following mesh on the left.
            . .                     .       Then a possible triIndices and halfEdges arrays structure 
            .   .                   .       is as follows:
            .     .                 .       
            .       .               .       halfEdges  := { -1, -1, 3, 2, -1, -1 }
            .         .             .       triIndices := {  0,  2, 3, 0,  3,  1 }
            .           .           .
            .             .         .       An additional constraint is that the ordering of indices in
            .               .       .       triIndices gives a CCW orientation.
            .                 .     .
            .                   .   .
            .                     . .
            2.......................3
    */
    class Mesh
    {
    public:
        /// Clears the internal state
        void Clear();

        /// Creates a mesh from an ObjParser object or a .obj filepath
        bool CreateFromObj(const char* objPath_);

        /// > 0 implies outside mesh
        f32 SignedDistanceAtPt(const Eigen::Vector3f& pt_);
        f32 SignedDistanceAtPt(const Eigen::Vector3f& pt_, const BVH& bvh_, const u32 threadIdx_ = 0);

        /// Returns a bounding volume for the mesh
        Eigen::AlignedBox3f CalculateMeshAABB() const;

    private:
        friend class BVH;

        /// Fills in the halfEdges array
        bool CreateHalfEdges();

        /// Finds the closest triangle to the point in O(n) time
        ClosestSimplexInfo ClosestTriangleToPt(const Eigen::Vector3f& pt_, u32& closestTriIdx_) const;

        /// https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.107.9173&rep=rep1&type=pdf
        Eigen::Vector3f PseudoNormal(const u32 triIndex_, const Simplex simplex_, const SimplexID simplexIdx_);
        Eigen::Vector3f PseudoNormalFace(const u32 triIndex_);
        Eigen::Vector3f PseudoNormalEdge(const u32 triIndex_, const SimplexID simplexIdx_);
        Eigen::Vector3f PseudoNormalVertex(const u32 triIndex_, const SimplexID simplexIdx_);

        std::vector<u32>             halfEdges;
        std::vector<u32>             triIndices;
        std::vector<Eigen::Vector3f> vertexNormals;
        std::vector<Eigen::Vector3f> vertices;
    };
}
