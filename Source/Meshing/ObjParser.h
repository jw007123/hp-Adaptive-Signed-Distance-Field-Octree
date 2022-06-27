#pragma once

#include <fstream>
#include <vector>

#include "Eigen/Core"

#include "../Literals.h"

namespace Meshing
{
    class ObjParser
    {
    public:
        ObjParser();
        ~ObjParser();

        /// Loads a set of vertices, normals and triangles into the object
        void Load(const char* objPath_);

        /// Clears the internal state
        void Clear();

        /// Is true if the internal arrays contain valid mesh data
        const bool HasValidData() const { return hasValidData; }

        /// Allows access to the internal arrays once Load has successfully returned
        const std::vector<u32>&             GetTriIndices()    const { return triIndices; }
        const std::vector<Eigen::Vector3f>& GetVertices()      const { return vertices; }
        const std::vector<Eigen::Vector3f>& GetVertexNormals() const { return vertexNormals; }

    private:
        enum VertFlags : u8
        {
            V  = 0,
            VN = 1,
            VT = 2,
            Num
        };
        u8 vertInfoFlags[3];

        bool hasValidData;

        std::vector<u32>             triIndices;
        std::vector<Eigen::Vector3f> vertices;
        std::vector<Eigen::Vector3f> vertexNormals;

        /// Parses a single line of the file
        void ParseLine(const char* lineBuff_);

        /// Calculates vertex normals from just vertices and indices alone
        void CalculateVertexNormals();
    };
}