#include "ObjParser.h"

namespace Meshing
{
    ObjParser::ObjParser()
    {
        memset(vertInfoFlags, 0, sizeof(vertInfoFlags));
    }


    ObjParser::~ObjParser()
    {

    }


    bool ObjParser::Load(const char* objPath_)
    {
        Clear();

        std::ifstream fStream;
        constexpr usize buffSize = 1024;
        char lineBuff[buffSize];

        fStream.open(objPath_, std::ios::in);
        while (fStream)
        {
            fStream.getline(lineBuff, buffSize);
            ParseLine(lineBuff);
        }
        fStream.close();

        // Fill in vertex normals
        if (vertices.size() && triIndices.size())
        {
            CalculateVertexNormals();
        }

        // Everything worked?
        return (vertices.size() && triIndices.size() && vertexNormals.size());
    }


    void ObjParser::Clear()
    {
        memset(vertInfoFlags, 0, sizeof(vertInfoFlags));

        vertices.clear();
        vertexNormals.clear();
        triIndices.clear();
    }


    void ObjParser::ParseLine(const char* lineBuff_)
    {
        char type[8] = { '\0' };
        u8 pos = 0;

        // Determine identifier
        while (lineBuff_[pos] != ' ')
        {
            type[pos] = lineBuff_[pos];
            pos++;
        }
        
        // Only interested in 'v' and 'f'
        if (!strcmp(type, "v"))
        {
            Eigen::Vector3f newVertex;

            std::sscanf(lineBuff_ + 2, "%f %f %f", &newVertex.x(), &newVertex.y(), &newVertex.z());
            vertices.push_back(newVertex);

            vertInfoFlags[VertFlags::V] = 1;
        }
        else if (!strcmp(type, "vn"))
        {
            vertInfoFlags[VertFlags::VN] = 1;
        }
        else if (!strcmp(type, "vt"))
        {
            vertInfoFlags[VertFlags::VT] = 1;
        }
        else if (!strcmp(type, "f"))
        {
            u8 formatType = 0;
            for (u8 i = 0; i < VertFlags::Num; ++i)
            {
                formatType += vertInfoFlags[i];
            }

            switch (formatType)
            {
                case 1:
                {

                    Eigen::Vector<u32, 3> newIndices;
                    std::sscanf(lineBuff_ + 2, "%lu %lu %lu", &newIndices.coeffRef(0),
                                                              &newIndices.coeffRef(1),
                                                              &newIndices.coeffRef(2));

                    triIndices.push_back(newIndices.coeff(0) - 1);
                    triIndices.push_back(newIndices.coeff(1) - 1);
                    triIndices.push_back(newIndices.coeff(2) - 1);


                    break;
                }

                case 2:
                {
                    Eigen::Vector<u32, 6> newIndices;
                    std::sscanf(lineBuff_ + 2, "%lu//%lu %lu//%lu %lu//%lu",
                                                &newIndices.coeffRef(0), &newIndices.coeffRef(1), &newIndices.coeffRef(2),
                                                &newIndices.coeffRef(3), &newIndices.coeffRef(4), &newIndices.coeffRef(5));

                    triIndices.push_back(newIndices.coeff(0) - 1);
                    triIndices.push_back(newIndices.coeff(2) - 1);
                    triIndices.push_back(newIndices.coeff(4) - 1);

                    break;
                }

                case 3:
                {
                    Eigen::Vector<u32, 9> newIndices;
                    std::sscanf(lineBuff_ + 2, "%lu/%lu/%lu %lu/%lu/%lu %lu/%lu/%lu",
                                                &newIndices.coeffRef(0), &newIndices.coeffRef(1), &newIndices.coeffRef(2),
                                                &newIndices.coeffRef(3), &newIndices.coeffRef(4), &newIndices.coeffRef(5),
                                                &newIndices.coeffRef(6), &newIndices.coeffRef(7), &newIndices.coeffRef(8));

                    triIndices.push_back(newIndices.coeff(0) - 1);
                    triIndices.push_back(newIndices.coeff(3) - 1);
                    triIndices.push_back(newIndices.coeff(6) - 1);

                    break;
                }

                default:
                    assert(0);
            }
        }
    }

  
    void ObjParser::CalculateVertexNormals()
    {
        const u32 nTris = (u32)triIndices.size() / 3;
        
        // Alloc and set all to 0
        vertexNormals.resize(vertices.size());
        memset(vertexNormals.data(), 0, sizeof(Eigen::Vector3f) * vertices.size());

        for (u32 i = 0; i < nTris; ++i)
        {
            const Eigen::Vector3f ab = vertices[triIndices[i * 3 + 1]] - vertices[triIndices[i * 3 + 0]];
            const Eigen::Vector3f ac = vertices[triIndices[i * 3 + 2]] - vertices[triIndices[i * 3 + 0]];
            const Eigen::Vector3f n  = (ab.cross(ac)).normalized();

            vertexNormals[triIndices[i * 3 + 0]] += n;
            vertexNormals[triIndices[i * 3 + 1]] += n;
            vertexNormals[triIndices[i * 3 + 2]] += n;
        }

        for (u32 i = 0; i < vertices.size(); ++i)
        {
            vertexNormals[i].normalize();
        }
    }
}