#include "Benchmarks.h"

namespace SDF
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
                case Benchmark::Creation:
                {
                    timeDiff = BenchmarkOctreeCreation();
                    break;
                }

                case Benchmark::Continuity:
                {
                    timeDiff = BenchmarkOctreeContinuity();
                    break;
                }

                case Benchmark::QueryDistance:
                {
                    timeDiff = BenchmarkOctreeDistanceQuerying();
                    break;
                }

                case Benchmark::QueryNormal:
                {
                    timeDiff = BenchmarkOctreeNormalQuerying();
                    break;
                }

                case Benchmark::SDFOperations:
                {
                    timeDiff = BenchmarkOctreeSDFOperations();
                    break;
                }

                default:
                    assert(0);
            }

            printf("\n\nCompleted %s in %lfs", BenchmarkStrings[i], timeDiff);
        }
    }


    f64 Benchmarks::BenchmarkOctreeCreation()
    {
        const std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();

        auto SphereFunc = [](const Eigen::Vector3d& pt_, const u32 threadIdx_) -> f64
        {
            return (pt_ - Eigen::Vector3d(0.25, 0, 0)).norm() - 0.5;
        };

        Config hpConfig;
        hpConfig.targetErrorThreshold       = pow(10, -10);
        hpConfig.nearnessWeighting.type     = Config::NearnessWeighting::Type::Exponential;
        hpConfig.nearnessWeighting.strength = 3.0;
        hpConfig.continuity.enforce         = false;
        hpConfig.threadCount                = std::thread::hardware_concurrency() != 0 ? std::thread::hardware_concurrency() : 1;

        Octree hpOctree;
        hpOctree.Create(hpConfig, SphereFunc);

        const std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<f64> timeDiff                                 = endTime - startTime;

        return timeDiff.count();
    }


    f64 Benchmarks::BenchmarkOctreeContinuity()
    {
        const std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();

        auto SphereFunc = [](const Eigen::Vector3d& pt_, const u32 threadIdx_) -> f64
        {
            return (pt_ - Eigen::Vector3d(0.25, 0, 0)).norm() - 0.5;
        };

        Config hpConfig;
        hpConfig.targetErrorThreshold       = pow(10, -10);
        hpConfig.nearnessWeighting.type     = Config::NearnessWeighting::Type::Exponential;
        hpConfig.nearnessWeighting.strength = 3.0;
        hpConfig.continuity.enforce         = true;
        hpConfig.continuity.strength        = 8.0;
        hpConfig.threadCount                = std::thread::hardware_concurrency() != 0 ? std::thread::hardware_concurrency() : 1;

        Octree hpOctree;
        hpOctree.Create(hpConfig, SphereFunc);

        const std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<f64> timeDiff                                 = endTime - startTime;

        return timeDiff.count();
    }


    f64 Benchmarks::BenchmarkOctreeDistanceQuerying()
    {
        auto SphereFunc = [](const Eigen::Vector3d& pt_, const u32 threadIdx_) -> f64
        {
            return (pt_ - Eigen::Vector3d(0.25, 0, 0)).norm() - 0.5;
        };

        Config hpConfig;
        hpConfig.targetErrorThreshold       = pow(10, -10);
        hpConfig.nearnessWeighting.type     = Config::NearnessWeighting::Type::Exponential;
        hpConfig.nearnessWeighting.strength = 3.0;
        hpConfig.continuity.enforce         = false;
        hpConfig.threadCount                = std::thread::hardware_concurrency() != 0 ? std::thread::hardware_concurrency() : 1;

        Octree hpOctree;
        hpOctree.Create(hpConfig, SphereFunc);

        const std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();

        const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.5, -0.5, -0.5), Eigen::Vector3d(0.5, 0.5, 0.5));
        for (usize i = 0; i < 10000000; ++i)
        {
            const Eigen::Vector3d sample(box.sample());
            const f64 octS = hpOctree.Query(sample);
        }

        const std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<f64> timeDiff                                 = endTime - startTime;

        return timeDiff.count();
    }


    f64 Benchmarks::BenchmarkOctreeNormalQuerying()
    {
        auto SphereFunc = [](const Eigen::Vector3d& pt_, const u32 threadIdx_) -> f64
        {
            return (pt_ - Eigen::Vector3d(0.25, 0, 0)).norm() - 0.5;
        };

        Config hpConfig;
        hpConfig.targetErrorThreshold       = pow(10, -10);
        hpConfig.nearnessWeighting.type     = Config::NearnessWeighting::Type::Exponential;
        hpConfig.nearnessWeighting.strength = 3.0;
        hpConfig.continuity.enforce         = false;
        hpConfig.threadCount                = std::thread::hardware_concurrency() != 0 ? std::thread::hardware_concurrency() : 1;

        Octree hpOctree;
        hpOctree.Create(hpConfig, SphereFunc);

        const std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();

        const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.5, -0.5, -0.5), Eigen::Vector3d(0.5, 0.5, 0.5));
        for (usize i = 0; i < 10000000; ++i)
        {
            const Eigen::Vector3d sample(box.sample());

            Eigen::Vector3d octGrad;
            const f64 octS = hpOctree.QueryWithGradient(sample, octGrad);
        }

        const std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<f64> timeDiff                                 = endTime - startTime;

        return timeDiff.count();
    }


    f64 Benchmarks::BenchmarkOctreeSDFOperations()
    {
        auto SphereFunc = [](const Eigen::Vector3d& pt_, const u32 threadIdx_) -> f64
        {
            return (pt_ - Eigen::Vector3d(0.25, 0, 0)).norm() - 0.5;
        };

        auto OtherSphereFunc = [](const Eigen::Vector3d& pt_, const u32 threadIdx_) -> f64
        {
            return (pt_ + Eigen::Vector3d(0.25, 0, 0)).norm() - 0.5;
        };

        Config hpConfig;
        hpConfig.targetErrorThreshold       = pow(10, -9);
        hpConfig.nearnessWeighting.type     = Config::NearnessWeighting::Type::Exponential;
        hpConfig.nearnessWeighting.strength = 3.0;
        hpConfig.continuity.enforce         = false;
        hpConfig.threadCount                = std::thread::hardware_concurrency() != 0 ? std::thread::hardware_concurrency() : 1;

        Octree hpOctree;
        hpOctree.Create(hpConfig, SphereFunc);

        const std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();

        hpOctree.UnionSDF(OtherSphereFunc);

        const std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<f64> timeDiff                                 = endTime - startTime;

        return timeDiff.count();
    }
}