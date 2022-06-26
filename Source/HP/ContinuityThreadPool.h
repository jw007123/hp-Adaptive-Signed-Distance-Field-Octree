#pragma once

#include <queue>
#include <mutex>

#include <malloc.h>
#include <atomic>

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/Sparse"

#include "Literals.h"
#include "Config.h"

namespace SDF
{
	class Octree;

	class ContinuityThreadPool
	{
	public:
		struct Input
		{
			std::pair<usize, usize> nodeIdxs;
			u8                      dim;
		};

		struct InitialData
		{
			Config config;

			std::atomic<bool>* shutdownAtom;

			std::queue<Input>* inputQueue;
			std::mutex*        inputQueueMutex;

			std::vector<Eigen::Triplet<f64>>* matTriplets;
		};

		ContinuityThreadPool();
		~ContinuityThreadPool();

		void StartThreads(InitialData* const initialData_, Octree* const octree_);

		void StopThreads();

	private:
		std::vector<std::thread> threadPool;
	};
}