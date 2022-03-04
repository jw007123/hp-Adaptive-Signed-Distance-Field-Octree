#pragma once

#include "Octree.h"

namespace SDF
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
			Creation      = 0,
			Continuity    = 1,
			Serialisation = 2,
            Copying       = 3,
            SDFOperations = 4,
			Num
		};

		const char* TestStrings[Test::Num] = 
		{
			"Octree Creation: ",
			"Octree Continuity: ",
			"Octree Serialisation: ",
            "Octree Copying: ",
            "Octree SDF Operations: "
		};

		bool TestOctreeCreation();

		bool TestOctreeContinuity();
		
		bool TestOctreeSerialisation();

        bool TestOctreeCopying();

        bool TestOctreeSDFOperations();
	};
}