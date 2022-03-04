#include "UnitTests.h"

namespace SDF
{
	UnitTests::UnitTests()
	{

	}


	UnitTests::~UnitTests()
	{

	}


	bool UnitTests::Run()
	{
        usize testsPassed = 0;
		for (usize i = 0; i < Test::Num; ++i)
		{
			bool testPassed = false;

			switch (i)
			{
				case Test::Creation:
				{
					testPassed = TestOctreeCreation();
					break;
				}

				case Test::Continuity:
				{
					testPassed = TestOctreeContinuity();
					break;
				}

				case Test::Serialisation:
				{
					testPassed = TestOctreeSerialisation();
					break;
				}

                case Test::Copying:
                {
                    testPassed = TestOctreeCopying();
                    break;
                }

                case Test::SDFOperations:
                {
                    testPassed = TestOctreeSDFOperations();
                    break;
                }

                case Test::CustomDomains:
                {
                    testPassed = TestOctreeCustomDomains();
                    break;
                }

				default:
					assert(0);
			}

			if (testPassed)
			{
				printf("%sPassed\n\n", TestStrings[i]);
			}
			else
			{
				printf("%sFailed\n\n", TestStrings[i]);
			}

            testsPassed += testPassed;
		}

        if (testsPassed == Test::Num)
        {
            printf("All tests passed!");
        }
        else
        {
            printf("Some tests failed!");
        }

        return (testsPassed == Test::Num);
	}


	bool UnitTests::TestOctreeCreation()
	{
		auto SphereFunc = [](const Eigen::Vector3d& pt_) -> f64
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

		const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.5, -0.5, -0.5), Eigen::Vector3d(0.5, 0.5, 0.5));
		for (usize i = 0; i < 1000000; ++i)
		{
			const Eigen::Vector3d sample(box.sample());
			const f64 octS  = hpOctree.Query(sample);
			const f64 trueS = SphereFunc(sample);

			if (abs(octS - trueS) > 0.01)
			{
				return false;
			}
		}

		return true;
	}


	bool UnitTests::TestOctreeContinuity()
	{
		auto SphereFunc = [](const Eigen::Vector3d& pt_) -> f64
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

		const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.5, -0.5, -0.5), Eigen::Vector3d(0.5, 0.5, 0.5));
		for (usize i = 0; i < 1000000; ++i)
		{
			const Eigen::Vector3d sample(box.sample());
			const f64 octS  = hpOctree.Query(sample);
			const f64 trueS = SphereFunc(sample);

			if (abs(octS - trueS) > 0.01)
			{
				return false;
			}
		}

		return true;
	}


	bool UnitTests::TestOctreeSerialisation()
	{
		auto SphereFunc = [](const Eigen::Vector3d& pt_) -> f64
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

		MemoryBlock hpBlock;
		{
			Octree hpOctree;
			hpOctree.Create(hpConfig, SphereFunc);
			hpBlock = hpOctree.ToMemoryBlock();
		}
		Octree hpOctree;
		hpOctree.FromMemoryBlock(hpBlock);
		free(hpBlock.ptr);

		const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.5, -0.5, -0.5), Eigen::Vector3d(0.5, 0.5, 0.5));
		for (usize i = 0; i < 1000000; ++i)
		{
			const Eigen::Vector3d sample(box.sample());
			const f64 octS  = hpOctree.Query(sample);
			const f64 trueS = SphereFunc(sample);

			if (abs(octS - trueS) > 0.01)
			{
				return false;
			}
		}

		return true;
	}


    bool UnitTests::TestOctreeCopying()
    {
        auto SphereFunc = [](const Eigen::Vector3d& pt_) -> f64
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

        Octree otherOctree = hpOctree;

        const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.5, -0.5, -0.5), Eigen::Vector3d(0.5, 0.5, 0.5));
        for (usize i = 0; i < 1000000; ++i)
        {
            const Eigen::Vector3d sample(box.sample());
            const f64 octS  = otherOctree.Query(sample);
            const f64 trueS = SphereFunc(sample);

            if (abs(octS - trueS) > 0.01)
            {
                return false;
            }
        }

        otherOctree = std::move(hpOctree);

        for (usize i = 0; i < 1000000; ++i)
        {
            const Eigen::Vector3d sample(box.sample());
            const f64 octS  = otherOctree.Query(sample);
            const f64 trueS = SphereFunc(sample);

            if (abs(octS - trueS) > 0.01)
            {
                return false;
            }
        }

        return true;
    }


    bool UnitTests::TestOctreeSDFOperations()
    {
        auto SphereFunc = [](const Eigen::Vector3d& pt_) -> f64
        {
            return (pt_ - Eigen::Vector3d(0.25, 0, 0)).norm() - 0.5;
        };

        auto OtherSphereFunc = [](const Eigen::Vector3d& pt_) -> f64
        {
            return (pt_ + Eigen::Vector3d(0.25, 0, 0)).norm() - 0.5;
        };

        Config hpConfig;
        hpConfig.targetErrorThreshold       = pow(10, -9);
        hpConfig.nearnessWeighting.type     = Config::NearnessWeighting::Type::Exponential;
        hpConfig.nearnessWeighting.strength = 3.0;
        hpConfig.continuity.enforce         = false;
        hpConfig.threadCount                = std::thread::hardware_concurrency() != 0 ? std::thread::hardware_concurrency() : 1;

        {
            Octree hpOctree;
            hpOctree.Create(hpConfig, SphereFunc);
            hpOctree.UnionSDF(OtherSphereFunc);

            const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.5, -0.5, -0.5), Eigen::Vector3d(0.5, 0.5, 0.5));
            for (usize i = 0; i < 1000000; ++i)
            {
                const Eigen::Vector3d sample(box.sample());
                const f64 octS  = hpOctree.Query(sample);
                const f64 trueS = std::min(SphereFunc(sample), OtherSphereFunc(sample));

                if (abs(octS - trueS) > 0.01)
                {
                    return false;
                }
            }
        }
       
        {
            Octree hpOctree;
            hpOctree.Create(hpConfig, SphereFunc);
            hpOctree.IntersectSDF(OtherSphereFunc);

            const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.5, -0.5, -0.5), Eigen::Vector3d(0.5, 0.5, 0.5));
            for (usize i = 0; i < 1000000; ++i)
            {
                const Eigen::Vector3d sample(box.sample());
                const f64 octS  = hpOctree.Query(sample);
                const f64 trueS = std::max(SphereFunc(sample), OtherSphereFunc(sample));

                if (abs(octS - trueS) > 0.01)
                {
                    return false;
                }
            }
        }

        {
            Octree hpOctree;
            hpOctree.Create(hpConfig, SphereFunc);
            hpOctree.SubtractSDF(OtherSphereFunc);

            const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.5, -0.5, -0.5), Eigen::Vector3d(0.5, 0.5, 0.5));
            for (usize i = 0; i < 1000000; ++i)
            {
                const Eigen::Vector3d sample(box.sample());
                const f64 octS  = hpOctree.Query(sample);
                const f64 trueS = std::max(SphereFunc(sample) * -1.0, OtherSphereFunc(sample));

                if (abs(octS - trueS) > 0.01)
                {
                    return false;
                }
            }
        }

        return true;
    }


    bool UnitTests::TestOctreeCustomDomains()
    {
        auto SphereFunc = [](const Eigen::Vector3d& pt_) -> f64
        {
            return (pt_ - Eigen::Vector3d(0.25, 0, 0)).norm() - 0.75;
        };

        Config hpConfig;
        hpConfig.targetErrorThreshold       = pow(10, -10);
        hpConfig.nearnessWeighting.type     = Config::NearnessWeighting::Type::Exponential;
        hpConfig.nearnessWeighting.strength = 3.0;
        hpConfig.continuity.enforce         = true;
        hpConfig.continuity.strength        = 8.0;
        hpConfig.threadCount                = std::thread::hardware_concurrency() != 0 ? std::thread::hardware_concurrency() : 1;
        hpConfig.root                       = Eigen::AlignedBox3f(Eigen::Vector3f(-0.25, -0.25, -0.25), Eigen::Vector3f(0.5, 0.5, 0.5));

        Octree hpOctree;
        hpOctree.Create(hpConfig, SphereFunc);

        const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.25, -0.25, -0.25), Eigen::Vector3d(0.5, 0.5, 0.5));
        for (usize i = 0; i < 1000000; ++i)
        {
            const Eigen::Vector3d sample(box.sample());
            const f64 octS = hpOctree.Query(sample);
            const f64 trueS = SphereFunc(sample);

            if (abs(octS - trueS) > 0.01)
            {
                return false;
            }
        }

        return true;
    }
}