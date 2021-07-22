#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "Config.cpp"
#include "Node.cpp"
#include "Octree.cpp"
#include "ContinuityThreadPool.cpp"
#include "BuildThreadPool.cpp"

int main()
{
	auto Figure1Func = [](const Eigen::Vector3d& pt_) -> double
	{
		double minVal = std::numeric_limits<double>::max();
		int idx = 0;
		for (double x = -1; x < 1; x += 0.2)
		{
			for (double y = -1; y < 1; y += 0.2)
			{
				const Eigen::Vector3d ptT = pt_ + Eigen::Vector3d(x, y, 0.0);
				const int idxM3 = idx % 3;
				double val = 0.0;

				if (idxM3 == 0)
				{
					// Sphere
					val = ptT.norm() - 0.07;
				}
				else if (idxM3 == 1)
				{
					// Cube
					const Eigen::Vector3d q = (ptT).cwiseAbs() - Eigen::Vector3d::Ones() * 0.07;
					val = q.cwiseMax(0.0).norm() + std::min<double>(std::max<double>(q.x(), std::max<double>(q.y(), q.z())), 0.0);
				}
				else
				{
					// Torus
					const Eigen::Vector2d b(Eigen::Vector2d(ptT.x(), ptT.y()).norm() - 0.06, ptT.z());
					val = b.norm() - 0.0125;
				}

				minVal = std::min<double>(val, minVal);
				idx++;
			}
		}

		return minVal;
	};

	SDF::Config hpConfig;
	hpConfig.targetErrorThreshold = pow(10, -6);
	hpConfig.nearnessWeighting.type = SDF::Config::NearnessWeighting::Type::Exponential;
	hpConfig.nearnessWeighting.strength = 3.0;
	hpConfig.continuity.enforce = false;
	hpConfig.continuity.strength = 8.0;
	hpConfig.threadCount = 12;

	SDF::Octree hpOctree;
	hpOctree.Create(hpConfig, Figure1Func);
	hpOctree.OutputFunctionSlice("SDF at Z = 0", 0.0);

	return 0;
}