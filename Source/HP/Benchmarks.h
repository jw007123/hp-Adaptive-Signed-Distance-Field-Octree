#pragma once

#include "Octree.h"

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
            Creation      = 0,
            Continuity    = 1,
            QueryDistance = 2,
            QueryNormal   = 3,
            SDFOperations = 4,
            Num
        };

        const char* BenchmarkStrings[Benchmark::Num] =
        {
            "Octree Creation",
            "Octree Continuity",
            "Octree Distance Querying",
            "Octree Normal Querying",
            "Octree SDF Operations"
        };

        f64 BenchmarkOctreeCreation();

        f64 BenchmarkOctreeContinuity();

        f64 BenchmarkOctreeDistanceQuerying();

        f64 BenchmarkOctreeNormalQuerying();

        f64 BenchmarkOctreeSDFOperations();
    };
}