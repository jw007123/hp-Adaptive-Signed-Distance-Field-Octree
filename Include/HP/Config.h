#pragma once

#include <thread>
#include <cassert>

#include "Eigen/Geometry"

#include "Utility/Literals.h"

namespace SDF
{
    struct Config
    {
        /// Reasonable defaults
        Config();

        struct NearnessWeighting
        {
	        enum Type : u8
	        {
		        None        = 0,
		        Polynomial  = 1,
		        Exponential = 2
	        } type;

	        f64 strength;
        } nearnessWeighting;

        struct Continuity
        {
	        bool enforce;
	        f64  strength;
        } continuity;

        bool enableLogging;
	    f64  targetErrorThreshold;
	    u32  threadCount;

        Eigen::AlignedBox3f root;

	    /// Checks that the config has valid parameters
	    void IsValid() const;
    };
}
