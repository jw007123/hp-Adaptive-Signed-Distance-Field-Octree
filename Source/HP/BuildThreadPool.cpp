#include "HP/BuildThreadPool.h"
#include "HP/Octree.h"

namespace SDF
{
	BuildThreadPool::BuildThreadPool()
	{

	}


	BuildThreadPool::~BuildThreadPool()
	{

	}


	void BuildThreadPool::StartThreads(const InitialData* const initialData_, Octree* const octree_)
	{
		assert(initialData_->config.threadCount > 0);

		threadPool.reserve(initialData_->config.threadCount);
		for (usize i = 0; i < initialData_->config.threadCount; ++i)
		{
			threadPool.push_back(std::thread(&Octree::TickBuildThread, octree_, initialData_, i));
		}
	}


	void BuildThreadPool::StopThreads()
	{
		for (u32 i = 0; i < threadPool.size(); ++i)
		{
			threadPool[i].join();
		}

		threadPool.clear();
	}
}