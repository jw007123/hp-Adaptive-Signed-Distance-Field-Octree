#pragma once

#include "ObjParser.h"
#include "BVH.h"
#include "Mesh.h"

namespace Meshing
{
	class UnitTests
	{
	public:
		UnitTests();
		~UnitTests();

		bool Run();

	private:
		enum Test : u8
		{
			Obj     = 0,
            Meshing = 1,
            Octree  = 2,
            BVHB    = 3,
            BVHQ    = 4,
            Num
		};

		const char* TestStrings[Test::Num] = 
		{
			"Obj Parsing",
            "Mesh Creation",
            "NNOctree Querying",
            "BVH Building",
            "BVH Querying"
		};

        bool TestObjParsing();

        bool TestMeshCreation();

        bool TestNNOctreeQuerying();

        bool TestBVHBuilding();

        bool TestBVHQuerying();
	};
}