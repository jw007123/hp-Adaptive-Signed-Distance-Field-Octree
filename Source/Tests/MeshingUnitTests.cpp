#include "MeshingUnitTests.h"

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
        std::function<bool()> TestFuncs[Test::Num] =
        {
            std::bind(&UnitTests::TestObjParsing,       this),
            std::bind(&UnitTests::TestMeshCreation,     this),
            std::bind(&UnitTests::TestNNOctreeQuerying, this),
            std::bind(&UnitTests::TestBVHBuilding,      this),
            std::bind(&UnitTests::TestBVHQuerying,      this)
        };

        usize testsPassed = 0;
        for (usize i = 0; i < Test::Num; ++i)
        {
            const bool testPassed = TestFuncs[i]();
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
            printf("All tests passed!\n\n");
        }
        else
        {
            printf("Some tests failed!\n\n");
        }

        return (testsPassed == Test::Num);
    }


    bool UnitTests::TestObjParsing()
    {
        ObjParser objParser;
        return objParser.Load(".\\Resources\\Ramesses.obj");
    }


    bool UnitTests::TestMeshCreation()
    {
        Mesh objMesh;
        return objMesh.CreateFromObj(".\\Resources\\Ramesses.obj");
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
        if (!objMesh.CreateFromObj(".\\Resources\\Ramesses.obj"))
        {
            return false;
        }

        BVH objBVH;
        if (!objBVH.Create(objMesh))
        {
            return false;
        }

        return true;
    }


    bool UnitTests::TestBVHQuerying()
    {
        Mesh objMesh;
        if (!objMesh.CreateFromObj(".\\Resources\\Ramesses.obj"))
        {
            return false;
        }

        BVH objBVH;
        if (!objBVH.Create(objMesh))
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