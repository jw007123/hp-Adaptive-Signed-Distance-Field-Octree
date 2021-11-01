#include "UnitTests.h"

namespace SDF
{
	UnitTests::UnitTests()
	{

	}


	UnitTests::~UnitTests()
	{

	}


	void UnitTests::Run()
	{
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

				default:
					assert(0);
			}

			if (testPassed)
			{
				printf("\n%sPassed", TestStrings[i]);
			}
			else
			{
				printf("\n%sFailed", TestStrings[i]);
			}
		}

        std::this_thread::sleep_for(std::chrono::seconds(5));
	}


	bool UnitTests::TestOctreeCreation()
	{
		auto SphereFunc = [](const Eigen::Vector3d& pt_) -> double
		{
			return (pt_ - Eigen::Vector3d(0.25, 0, 0)).norm() - 0.5;
		};

		Config hpConfig;
		hpConfig.targetErrorThreshold = pow(10, -10);
		hpConfig.nearnessWeighting.type = Config::NearnessWeighting::Type::Exponential;
		hpConfig.nearnessWeighting.strength = 3.0;
		hpConfig.continuity.enforce = false;
        hpConfig.threadCount = std::thread::hardware_concurrency() != 0 ? std::thread::hardware_concurrency() : 1;

		Octree hpOctree;
		hpOctree.Create(hpConfig, SphereFunc);

		const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.5, -0.5, -0.5), Eigen::Vector3d(0.5, 0.5, 0.5));
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


	bool UnitTests::TestOctreeContinuity()
	{
		auto SphereFunc = [](const Eigen::Vector3d& pt_) -> double
		{
			return (pt_ - Eigen::Vector3d(0.25, 0, 0)).norm() - 0.5;
		};

		Config hpConfig;
		hpConfig.targetErrorThreshold = pow(10, -10);
		hpConfig.nearnessWeighting.type = Config::NearnessWeighting::Type::Exponential;
		hpConfig.nearnessWeighting.strength = 3.0;
		hpConfig.continuity.enforce = true;
		hpConfig.continuity.strength = 8.0;
        hpConfig.threadCount = std::thread::hardware_concurrency() != 0 ? std::thread::hardware_concurrency() : 1;

		Octree hpOctree;
		hpOctree.Create(hpConfig, SphereFunc);

		const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.5, -0.5, -0.5), Eigen::Vector3d(0.5, 0.5, 0.5));
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


	bool UnitTests::TestOctreeSerialisation()
	{
		auto SphereFunc = [](const Eigen::Vector3d& pt_) -> double
		{
			return (pt_ - Eigen::Vector3d(0.25, 0, 0)).norm() - 0.5;
		};

		Config hpConfig;
		hpConfig.targetErrorThreshold = pow(10, -10);
		hpConfig.nearnessWeighting.type = Config::NearnessWeighting::Type::Exponential;
		hpConfig.nearnessWeighting.strength = 3.0;
		hpConfig.continuity.enforce = true;
		hpConfig.continuity.strength = 8.0;
        hpConfig.threadCount = std::thread::hardware_concurrency() != 0 ? std::thread::hardware_concurrency() : 1;

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
			const f64 octS = hpOctree.Query(sample);
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
        auto SphereFunc = [](const Eigen::Vector3d& pt_) -> double
        {
            return (pt_ - Eigen::Vector3d(0.25, 0, 0)).norm() - 0.5;
        };

        Config hpConfig;
        hpConfig.targetErrorThreshold = pow(10, -10);
        hpConfig.nearnessWeighting.type = Config::NearnessWeighting::Type::Exponential;
        hpConfig.nearnessWeighting.strength = 3.0;
        hpConfig.continuity.enforce = false;
        hpConfig.threadCount = std::thread::hardware_concurrency() != 0 ? std::thread::hardware_concurrency() : 1;

        Octree hpOctree;
        hpOctree.Create(hpConfig, SphereFunc);

        Octree otherOctree = hpOctree;

        const Eigen::AlignedBox3d box(Eigen::Vector3d(-0.5, -0.5, -0.5), Eigen::Vector3d(0.5, 0.5, 0.5));
        for (usize i = 0; i < 1000000; ++i)
        {
            const Eigen::Vector3d sample(box.sample());
            const f64 octS = otherOctree.Query(sample);
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
            const f64 octS = otherOctree.Query(sample);
            const f64 trueS = SphereFunc(sample);

            if (abs(octS - trueS) > 0.01)
            {
                return false;
            }
        }

        return true;
    }
}