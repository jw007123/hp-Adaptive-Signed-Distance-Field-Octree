#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "Config.cpp"
#include "Node.cpp"
#include "Octree.cpp"
#include "ContinuityThreadPool.cpp"
#include "BuildThreadPool.cpp"

int main()
{
	auto SphereFunc = [](const Eigen::Vector3d& pt_) -> double
	{
		return pt_.norm() - 0.25;
	};

	SDF::Config hpConfig;
	hpConfig.targetErrorThreshold = pow(10, -6);
	hpConfig.nearnessWeighting.type = SDF::Config::NearnessWeighting::Type::Exponential;
	hpConfig.nearnessWeighting.strength = 3.0;
	hpConfig.continuity.enforce = true;
	hpConfig.continuity.strength = 8.0;
	hpConfig.threadCount = 12;

	SDF::Octree hpOctree;
	hpOctree.Create(hpConfig, SphereFunc);
	hpOctree.OutputFunctionSlice("SDF at Z = 0", 0.0);
	hpOctree.OutputTreeSlice("Octree at Z = 0", 0.0);

	return 0;
}