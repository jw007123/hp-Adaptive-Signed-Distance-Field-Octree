#include "HPBenchmarks.h"

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
        std::function<f64()> BenchFuncs[Benchmarks::Num] =
        {
            std::bind(&Benchmarks::BenchmarkOctreeCreation,                this),
            std::bind(&Benchmarks::BenchmarkOctreeContinuity,              this),
            std::bind(&Benchmarks::BenchmarkOctreeRandomDistanceQuerying,  this),
            std::bind(&Benchmarks::BenchmarkOctreeUniformDistanceQuerying, this),
            std::bind(&Benchmarks::BenchmarkOctreeNormalQuerying,          this),
            std::bind(&Benchmarks::BenchmarkOctreeSDFOperations,           this)
        };

        for (usize i = 0; i < Benchmarks::Num; ++i)
        {
            const f64 timeDiff = BenchFuncs[i]();
            printf("Completed %s in %lfs\n\n", BenchmarkStrings[i], timeDiff);
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


    f64 Benchmarks::BenchmarkOctreeRandomDistanceQuerying()
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

        const usize nQueries = 8000000;
        const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.5, -0.5, -0.5), Eigen::Vector3d(0.5, 0.5, 0.5));

        std::vector<Eigen::Vector3d> samples;
        samples.resize(nQueries);
        for (usize i = 0; i < nQueries; ++i)
        {
            samples.push_back(box.sample());
        }

        const std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();
        for (usize i = 0; i < samples.size(); ++i)
        {
            const f64 octS = hpOctree.Query(samples[i]);
        }

        const std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<f64> timeDiff                                 = endTime - startTime;

        return timeDiff.count();
    }


    f64 Benchmarks::BenchmarkOctreeUniformDistanceQuerying()
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

        const usize nQueries = 200;
        const f64 dt         = 1.0 / (nQueries + 1);
        const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.5, -0.5, -0.5), Eigen::Vector3d(0.5, 0.5, 0.5));

        std::vector<Eigen::Vector3d> samples;
        samples.resize(nQueries * nQueries * nQueries);
        for (usize i = 1; i <= nQueries; ++i)
        {
            for (usize j = 1; j <= nQueries; ++j)
            {
                for (usize k = 1; k <= nQueries; ++k)
                {
                    Eigen::Vector3d sample = box.min();
                    sample.x() += i * dt;
                    sample.y() += j * dt;
                    sample.z() += k * dt;

                    samples.push_back(sample);
                }
            }
        }

        const std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();
        for (usize i = 0; i < samples.size(); ++i)
        {
            const f64 octS = hpOctree.Query(samples[i]);
        }
        const std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<f64> timeDiff = endTime - startTime;

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
        for (usize i = 0; i < 8000000; ++i)
        {
            const Eigen::Vector3d sample(box.sample());

            // NOTE: Can't really do a unit test for normals as the normal is undefined at the SDF boundary...
            // i.e. we can always pick a point on the SDF where the ||approximated normal - 'true' normal|| > tol
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
        hpConfig.targetErrorThreshold       = pow(10, -8);
        hpConfig.nearnessWeighting.type     = Config::NearnessWeighting::Type::Exponential;
        hpConfig.nearnessWeighting.strength = 3.0;
        hpConfig.continuity.enforce         = false;
        hpConfig.threadCount                = std::thread::hardware_concurrency() != 0 ? std::thread::hardware_concurrency() : 1;

        Octree hpOctree;
        hpOctree.Create(hpConfig, SphereFunc);

        const std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();
        {
            hpOctree.UnionSDF(OtherSphereFunc);
        }
        const std::chrono::time_point<std::chrono::high_resolution_clock> endTime = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<f64> timeDiff                                 = endTime - startTime;

        return timeDiff.count();
    }
}