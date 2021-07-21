#pragma once

#include "Literals.h"

namespace SDF
{
	struct Config
	{
		struct NearnessWeighting
		{
			enum Type : u8
			{
				None = 0,
				Polynomial = 1,
				Exponential = 2
			} type;

			f64 strength;
		} nearnessWeighting;

		struct Continuity
		{
			bool enforce;
			f64 strength;
		} continuity;

		f64 targetErrorThreshold;
		usize threadCount;

		/// Checks that the config has valid parameters
		void IsValid() const;
	};
}
