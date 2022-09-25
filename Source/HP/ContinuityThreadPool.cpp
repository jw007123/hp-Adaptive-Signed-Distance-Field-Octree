#include "HP/ContinuityThreadPool.h"
#include "HP/Octree.h"

namespace SDF
{
    ContinuityThreadPool::ContinuityThreadPool()
    {

    }


    ContinuityThreadPool::~ContinuityThreadPool()
    {

    }


    void ContinuityThreadPool::StartThreads(InitialData* const initialData_, Octree* const octree_)
    {
	    assert(initialData_->config.threadCount > 0);

	    threadPool.reserve(initialData_->config.threadCount);
	    for (usize i = 0; i < initialData_->config.threadCount; ++i)
	    {
		    threadPool.push_back(std::thread(&Octree::TickContinuityThread, octree_, initialData_, i));
	    }
    }


    void ContinuityThreadPool::StopThreads()
    {
	    for (usize i = 0; i < threadPool.size(); ++i)
	    {
		    threadPool[i].join();
	    }

	    threadPool.clear();
    }
}