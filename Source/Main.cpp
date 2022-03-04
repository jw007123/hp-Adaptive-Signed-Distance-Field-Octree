#include "Ray.cpp"
#include "Config.cpp"
#include "Node.cpp"
#include "Octree.cpp"
#include "ContinuityThreadPool.cpp"
#include "BuildThreadPool.cpp"
#include "UnitTests.cpp"
#include "Benchmarks.cpp"

int main()
{
	SDF::UnitTests unitTests;
    if (!unitTests.Run())
    {
        std::this_thread::sleep_for(std::chrono::seconds(5));
        return -1;
    }

    SDF::Benchmarks benchmarks;
    benchmarks.Run();

    std::this_thread::sleep_for(std::chrono::seconds(5));

	return 0;
}