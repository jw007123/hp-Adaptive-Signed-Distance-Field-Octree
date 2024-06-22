#include "MeshingBenchmarks.h"

namespace Meshing
{
    void Benchmarks::Run()
    {
        std::function<f64()> BenchFuncs[Benchmarks::Num] =
        {
            std::bind(&Benchmarks::BenchmarkObjFileParsing,    this),
            std::bind(&Benchmarks::BenchmarkMeshFromObj,       this),
            std::bind(&Benchmarks::BenchmarkBVHFromMesh,       this),
            std::bind(&Benchmarks::BenchmarkBVHQuerying,       this),
            std::bind(&Benchmarks::BenchmarkNaiveMeshQuerying, this)
        };

        for (usize i = 0; i < Benchmarks::Num; ++i)
        {
            const f64 timeDiff = BenchFuncs[i]();
            printf("Completed %s in %lfs\n\n", BenchmarkStrings[i], timeDiff);
        }
    }


    f64 Benchmarks::BenchmarkObjFileParsing()
    {
        const std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();
        {
            ObjParser objParser;
            objParser.Load(".\\Resources\\Ramesses.obj");
        }
        const std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<f64> timeDiff                                 = endTime - startTime;

        return timeDiff.count();
    }


    f64 Benchmarks::BenchmarkMeshFromObj()
    {
        const std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();
        {
            Mesh objMesh;
            objMesh.CreateFromObj(".\\Ramesses.obj");
        }
        const std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<f64> timeDiff                                 = endTime - startTime;

        return timeDiff.count();
    }


    f64 Benchmarks::BenchmarkBVHFromMesh()
    {
        Mesh objMesh;
        objMesh.CreateFromObj(".\\Resources\\Ramesses.obj");

        const std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();
        {
            BVH objBVH;
            objBVH.Create(objMesh);
        }
        const std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<f64> timeDiff                                 = endTime - startTime;

        return timeDiff.count();
    }


    f64 Benchmarks::BenchmarkBVHQuerying()
    {
        Mesh objMesh;
        objMesh.CreateFromObj(".\\Resources\\Ramesses.obj");
        const Eigen::AlignedBox3f meshRoot = objMesh.CalculateMeshAABB();

        BVH objBVH;
        objBVH.Create(objMesh);

        const std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();
        {
            for (u32 i = 0; i < 10000; ++i)
            {
                const Eigen::Vector3f sample = meshRoot.sample();
                const f32 d                  = objMesh.SignedDistanceAtPt(sample, objBVH);
            }
        }
        const std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<f64> timeDiff                                 = endTime - startTime;

        return timeDiff.count();
    }


    f64 Benchmarks::BenchmarkNaiveMeshQuerying()
    {
        Mesh objMesh;
        objMesh.CreateFromObj(".\\Resources\\Ramesses.obj");
        const Eigen::AlignedBox3f meshRoot = objMesh.CalculateMeshAABB();

        const std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();
        {
            for (u32 i = 0; i < 100; ++i)
            {
                const Eigen::Vector3f sample = meshRoot.sample();
                const f32 d = objMesh.SignedDistanceAtPt(sample);
            }
        }
        const std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<f64> timeDiff                                 = endTime - startTime;

        return timeDiff.count();
    }
}
