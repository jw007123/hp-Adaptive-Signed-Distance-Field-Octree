#pragma once

#include "HP/Octree.h"

namespace SDF
{
    class UnitTests
    {
    public:
	    bool Run();

    private:
	    enum Test : u8
	    {
		    Creation       = 0,
		    Continuity     = 1,
		    Serialisation  = 2,
            Copying        = 3,
            SDFOperations  = 4,
            CustomDomains  = 5,
		    Num
	    };

	    const char* TestStrings[Test::Num] = 
	    {
		    "Octree Creation",
		    "Octree Continuity",
		    "Octree Serialisation",
            "Octree Copying",
            "Octree SDF Operations",
            "Octree Custom Domains",
	    };

	    bool TestOctreeCreation();
	    bool TestOctreeContinuity();
	    bool TestOctreeSerialisation();
        bool TestOctreeCopying();
        bool TestOctreeSDFOperations();
        bool TestOctreeCustomDomains();
    };
}
