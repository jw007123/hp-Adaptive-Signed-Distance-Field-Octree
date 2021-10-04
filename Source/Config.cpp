#include "Config.h"

namespace SDF
{
	void Config::IsValid() const
	{
		assert(targetErrorThreshold > 0.0);
		assert(threadCount > 0);

		if (nearnessWeighting.type != Config::NearnessWeighting::None)
		{
			assert(nearnessWeighting.strength > 0.0);
		}

		if (continuity.enforce)
		{
			assert(continuity.strength > 0.0);
		}
	}
}