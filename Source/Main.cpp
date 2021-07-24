#include "Config.cpp"
#include "Node.cpp"
#include "Octree.cpp"
#include "ContinuityThreadPool.cpp"
#include "BuildThreadPool.cpp"
#include "UnitTests.cpp"

#include <iostream>

int main()
{
	SDF::UnitTests unitTests;
	unitTests.Run();

	return 0;
}