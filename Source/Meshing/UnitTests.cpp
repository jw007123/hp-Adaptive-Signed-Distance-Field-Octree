#include "UnitTests.h"

namespace Meshing
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
                case Test::Obj:
                {
                    testPassed = TestObjParsing();
                    break;
                }

                case Test::Meshing:
                {
                    testPassed = TestMeshCreation();
                    break;
                }

                case Test::Octree:
                {
                    testPassed = TestNNOctreeQuerying();
                    break;
                }

                case Test::BVHB:
                {
                    testPassed = TestBVHBuilding();
                    break;
                }

                case Test::BVHQ:
                {
                    testPassed = TestBVHQuerying();
                    break;
                }

				default:
					assert(0);
			}

			if (testPassed)
			{
				printf("%s: Passed\n\n", TestStrings[i]);
			}
			else
			{
				printf("%s: Failed\n\n", TestStrings[i]);
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


    bool UnitTests::TestObjParsing()
    {
        ObjParser objParser;
        objParser.Load("C:\\VS2017Projects\\SDF\\hp-Adaptive-Signed-Distance-Field-Octree\\Resources\\Ramesses.obj");

        return objParser.HasValidData();
    }


    bool UnitTests::TestMeshCreation()
    {
        Mesh objMesh;
        objMesh.CreateFromObj("C:\\VS2017Projects\\SDF\\hp-Adaptive-Signed-Distance-Field-Octree\\Resources\\Ramesses.obj");

        return objMesh.HasValidState();
    }


    bool UnitTests::TestNNOctreeQuerying()
    {
        const Eigen::AlignedBox3f octRoot(Eigen::Vector3f::Zero(), Eigen::Vector3f::Ones());

        NNOctree octree;
        octree.Create(octRoot);

        std::vector<Eigen::Vector3f> rndPts;
        for (u32 i = 0; i < 100000; ++i)
        {
            rndPts.push_back(octRoot.sample());
            octree.InsertPoint(rndPts[i], i);
        }

        for (u32 i = 0; i < 100000; ++i)
        {
            // 10.0f guarantees should always find point if octree is working
            const std::pair<Eigen::Vector3f, u32> ret = octree.NearestPoint(rndPts[i], 10.0f);
            if (ret.second != i)
            {
                return false;
            }

            if (!octree.RemovePoint(rndPts[i]))
            {
                return false;
            }
        }

        return true;
    }


    bool UnitTests::TestBVHBuilding()
    {
        Mesh objMesh;
        objMesh.CreateFromObj("C:\\VS2017Projects\\SDF\\hp-Adaptive-Signed-Distance-Field-Octree\\Resources\\Ramesses.obj");
        if (!objMesh.HasValidState())
        {
            return false;
        }

        BVH objBVH;
        objBVH.Create(objMesh);
        if (!objBVH.HasValidState())
        {
            return false;
        }

        return true;
    }


    bool UnitTests::TestBVHQuerying()
    {
        Mesh objMesh;
        objMesh.CreateFromObj("C:\\VS2017Projects\\SDF\\hp-Adaptive-Signed-Distance-Field-Octree\\Resources\\Ramesses.obj");
        if (!objMesh.HasValidState())
        {
            return false;
        }

        BVH objBVH;
        objBVH.Create(objMesh);
        if (!objBVH.HasValidState())
        {
            return false;
        }

        const Eigen::AlignedBox3f meshRoot = objMesh.CalculateMeshAABB();
        for (u32 i = 0; i < 50; ++i)
        {
            const Eigen::Vector3f sample = meshRoot.sample();
            const f32 d1                 = objMesh.SignedDistanceAtPt(sample);
            const f32 d2                 = objMesh.SignedDistanceAtPt(sample, objBVH);

            if ((d1 - d2) * (d1 - d2) > EPSILON_F32)
            {
                return false;
            }
        }

        return true;
    }
}