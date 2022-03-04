#include "Config.h"

namespace SDF
{
    Config::Config()
    {
        targetErrorThreshold       = pow(10, -10);
        nearnessWeighting.type     = NearnessWeighting::Type::Exponential;
        nearnessWeighting.strength = 3.0;
        continuity.enforce         = true;
        continuity.strength        = 8.0;
        threadCount                = std::thread::hardware_concurrency() != 0 ? std::thread::hardware_concurrency() : 1;
        root                       = Eigen::AlignedBox3f(Eigen::Vector3f(-0.5, -0.5, -0.5), Eigen::Vector3f(0.5, 0.5, 0.5));
    }


    Config::~Config()
    {

    }


	void Config::IsValid() const
	{
		assert(targetErrorThreshold > 0.0);
		assert(threadCount > 0);
        assert(root.volume() > 0.0);
        assert(root.sizes()(0) == root.sizes()(1) && root.sizes()(1) == root.sizes()(2));
        
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