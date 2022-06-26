#pragma once
#include "../Literals.h"

#include "ObjParser.h"
#include "BVH.h"
#include "Mesh.h"

namespace Meshing
{
    class Benchmarks
    {
    public:
        Benchmarks();
        ~Benchmarks();

        void Run();

    private:
        enum Benchmark : u8
        {
            ObjParsing      = 0,
            MeshFromObj     = 1,
            BVHConstruction = 2,
            BVHQuerying     = 3,
            MeshQuerying    = 4,
            Num
        };

        const char* BenchmarkStrings[Benchmark::Num] =
        {
            "1.6M tri .obj File Parsing",
            "1.6M tri .obj Mesh Creation",
            "1.6M tri Mesh -> BVH",
            "1.6M tri BVH 10K Queries",
            "1.6M tri Mesh O(n) 100 Queries"
        };

        f64 BenchmarkObjFileParsing();

        f64 BenchmarkMeshFromObj();

        f64 BenchmarkBVHFromMesh();

        f64 BenchmarkBVHQuerying();

        f64 BenchmarkNaiveMeshQuerying();
    };
}