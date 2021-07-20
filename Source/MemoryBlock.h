#pragma once

#include <malloc.h>

#include "Literals.h"

namespace SDF
{
	struct MemoryBlock
	{
		usize size;
		void* ptr;
	};
}