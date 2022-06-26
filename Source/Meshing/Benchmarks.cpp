#include "Benchmarks.h"

namespace Meshing
{
    Benchmarks::Benchmarks()
    {

    }


    Benchmarks::~Benchmarks()
    {

    }


    void Benchmarks::Run()
    {
        usize benchmarksRan = 0;
        for (usize i = 0; i < Benchmarks::Num; ++i)
        {
            f64 timeDiff = 0.0;

            switch (i)
            {
                case Benchmark::ObjParsing:
                {
                    timeDiff = BenchmarkObjFileParsing();
                    break;
                }

                case Benchmark::MeshFromObj:
                {
                    timeDiff = BenchmarkMeshFromObj();
                    break;
                }

                case Benchmark::BVHConstruction:
                {
                    timeDiff = BenchmarkBVHFromMesh();
                    break;
                }

                case Benchmark::BVHQuerying:
                {
                    timeDiff = BenchmarkBVHQuerying();
                    break;
                }

                case Benchmark::MeshQuerying:
                {
                    timeDiff = BenchmarkNaiveMeshQuerying();
                    break;
                }

                default:
                    assert(0);
            }

            printf("\n\nCompleted %s in %lfs", BenchmarkStrings[i], timeDiff);
        }
    }


    f64 Benchmarks::BenchmarkObjFileParsing()
    {
        return 0.0;
    }


    f64 Benchmarks::BenchmarkMeshFromObj()
    {
        return 0.0;
    }


    f64 Benchmarks::BenchmarkBVHFromMesh()
    {
        return 0.0;
    }


    f64 Benchmarks::BenchmarkBVHQuerying()
    {
        return 0.0;
    }


    f64 Benchmarks::BenchmarkNaiveMeshQuerying()
    {
        return 0.0;
    }
}