#pragma once

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "Literals.h"
#include "Utility.h"

namespace SDF
{
	struct Node
	{
		Node();
		~Node();

		// Nodes are allocated in blocks of 8, so the ith child's ID is childIdx + i
		u32 childIdx;

		// Standard AABB store
		Eigen::AlignedBox3f aabb;

		// Stores LegendreCoeffientCount[degree] doubles
		struct Basis
		{
			union
			{
				f64* coeffs;
				usize coeffsStart;
			};
			u8 degree;
		} basis;

		// Since we're storing an odd number of bytes with basis.degree, we may as well store depth as it'll be padded in anyway
		u8 depth;
	};
}