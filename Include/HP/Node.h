#pragma once

#include "Eigen/Geometry"

#include "Utility/Literals.h"
#include "HP/Consts.h"

namespace SDF
{
    struct Node
    {
	    Node();

	    // Nodes are allocated in blocks of 8, so the ith child's ID is childIdx + i
	    u32 childIdx;

	    // Standard AABB store
	    Eigen::AlignedBox3f aabb;

	    // Stores LegendreCoeffientCount[degree] doubles
	    struct Basis
	    {
		    union
		    {
			    f64*  coeffs = nullptr;
			    usize coeffsStart;
		    };
		    u8 degree;
	    } basis;

	    // Since we're storing an odd number of bytes with basis.degree, we may as well store depth as it'll be padded in anyway
	    u8 depth;
    };
}
