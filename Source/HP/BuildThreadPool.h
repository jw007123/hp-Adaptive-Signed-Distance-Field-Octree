#pragma once

#include <functional>
#include <queue>
#include <mutex>
#include <atomic>

#include <malloc.h>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "Literals.h"

#include "Node.h"
#include "Config.h"

namespace SDF
{
	class Octree;

	class BuildThreadPool
	{
	public:
		struct Input
		{
			Node                node;
			std::pair<u32, f64> nodeIdxAndErr;
		};

		struct Output
		{
			Node::Basis nodeBasis;

			enum BasisType : u8
			{
				P = 0,
				H = 1
			} basisType;

			// If basisType == H...
			u8 childIdx;

            u32 nodeIdx;
			f64 initialErr;
			f64 newErr;
		};

		struct InitialData
		{
			Config config;

			std::atomic<bool>* shutdownAtom;

			std::queue<Input>* inputQueue;
			std::mutex*        inputQueueMutex;

			std::queue<Output>* outputQueue;
			std::mutex*         outputQueueMutex;
		};

		BuildThreadPool();
		~BuildThreadPool();

		void StartThreads(const InitialData* const initialData_, Octree* const octree_);
		void StopThreads();

	private:
		std::vector<std::thread> threadPool;
	};
}