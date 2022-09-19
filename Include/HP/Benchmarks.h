#pragma once

#include "Utility/Literals.h"
#include "HP/Octree.h"

namespace SDF
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
            Creation       = 0,
            Continuity     = 1,
            QueryDistanceR = 2,
            QueryDistanceU = 3,
            QueryNormal    = 4,
            SDFOperations  = 5,
            Num
        };

        const char* BenchmarkStrings[Benchmark::Num] =
        {
            "Octree Creation",
            "Octree Continuity",
            "Octree 8M Random Distance Querying",
            "Octree 8M Uniform Distance Querying",
            "Octree 8M Normal Querying",
            "Octree SDF Operations"
        };

        f64 BenchmarkOctreeCreation();
        f64 BenchmarkOctreeContinuity();
        f64 BenchmarkOctreeRandomDistanceQuerying();
        f64 BenchmarkOctreeUniformDistanceQuerying();
        f64 BenchmarkOctreeNormalQuerying();
        f64 BenchmarkOctreeSDFOperations();
    };
}