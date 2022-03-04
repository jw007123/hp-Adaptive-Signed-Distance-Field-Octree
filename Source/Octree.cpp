#include "Octree.h"

namespace SDF
{
	Octree::Octree()
	{
		coeffStore = nullptr;
	}


	Octree::~Octree()
	{
		if (coeffStore)
		{
			free(coeffStore);
		}
	}


    Octree::Octree(const Octree& other_)
    {
        // Easy to copy
        config             = other_.config;
        nodes              = other_.nodes;
        configRootCentre   = other_.configRootCentre;
        configRootInvSizes = other_.configRootInvSizes;

        // Count coeffs
        usize nCoeffs = 0;
        for (usize i = 0; i < nodes.size(); ++i)
        {
            if (nodes[i].basis.degree != (BASIS_MAX_DEGREE + 1))
            {
                nCoeffs += LegendreCoeffientCount[nodes[i].basis.degree];
            }
        }

        // Copy over
        coeffStore = (f64*)malloc(sizeof(f64) * nCoeffs);
        memcpy(coeffStore, other_.coeffStore, sizeof(f64) * nCoeffs);
    }


    Octree& Octree::operator=(const Octree& other_)
    {
        Clear();

        // Easy to copy
        config             = other_.config;
        nodes              = other_.nodes;
        configRootCentre   = other_.configRootCentre;
        configRootInvSizes = other_.configRootInvSizes;
        
        // Count coeffs
        usize nCoeffs = 0;
        for (usize i = 0; i < nodes.size(); ++i)
        {
            if (nodes[i].basis.degree != (BASIS_MAX_DEGREE + 1))
            {
                nCoeffs += LegendreCoeffientCount[nodes[i].basis.degree];
            }
        }

        // Copy over
        coeffStore = (f64*)malloc(sizeof(f64) * nCoeffs);
        memcpy(coeffStore, other_.coeffStore, sizeof(f64) * nCoeffs);

        return *this;
    }


    Octree::Octree(Octree&& other_)
    {
        config     = other_.config;
        nodes      = std::move(other_.nodes);
        coeffStore = other_.coeffStore;

        configRootCentre   = other_.configRootCentre;
        configRootInvSizes = other_.configRootInvSizes;

        other_.coeffStore = nullptr;
    }


    Octree& Octree::operator=(Octree&& other_)
    {
        Clear();

        config     = other_.config;
        nodes      = std::move(other_.nodes);
        coeffStore = other_.coeffStore;

        configRootCentre   = other_.configRootCentre;
        configRootInvSizes = other_.configRootInvSizes;

        other_.coeffStore = nullptr;

        return *this;
    }


	void Octree::UniformlyRefine(f64& totalCoeffError_)
	{
		// Perform a coarse refinement of 4 levels and fit a degree 2 polynomial at each level
		constexpr usize coarseDepth  = 4;
		constexpr usize coarseDegree = 2;

		struct Inserter
		{
			usize parentIdxs[coarseDepth] = { 0 };
			usize childNos[coarseDepth]   = { 0 };
			usize depth                   = 1;
		} inserter;

		while (1)
		{
			const usize curParentIdx = inserter.parentIdxs[inserter.depth - 1];
			const usize curNodeIdx   = nodes[curParentIdx].childIdx + inserter.childNos[inserter.depth - 1];

			if (inserter.depth < coarseDepth)
			{
				if (nodes[curNodeIdx].childIdx == -1)
				{
					// If not at max depth and children don't exist, subdivide node and then go into it
					Subdivide(curNodeIdx);

					inserter.parentIdxs[inserter.depth] = curNodeIdx;
					inserter.childNos[inserter.depth]   = 0;
					inserter.depth++;
				}
				else
				{
					if (inserter.childNos[inserter.depth - 1] < 7)
					{
						// If children exist, try and go to next child
						inserter.childNos[inserter.depth - 1]++;
					}
					else
					{
						// Otherwise, go up again
						if (inserter.depth == 1)
						{
							break;
						}
						else
						{
							inserter.depth--;
						}
					}
				}
			}
			else
			{
				// Alloc basis mem
				nodes[curNodeIdx].basis.coeffs = (f64*)malloc(sizeof(f64) * LegendreCoeffientCount[coarseDegree]);
				nodes[curNodeIdx].basis.degree = 0;

				// Otherwise, fit polynomial to node
				const f64 newError = FitPolynomial(nodes[curNodeIdx].basis, nodes[curNodeIdx].aabb, coarseDegree, coarseDepth);
				totalCoeffError_  += newError;

				// Add to queue
				std::pair<usize, f64> newNodeIdxAndErr = { curNodeIdx, newError };
				nodeQueue.push(newNodeIdxAndErr);

				if (inserter.childNos[inserter.depth - 1] < 7)
				{
					// If children exist, try and go to next child
					inserter.childNos[inserter.depth - 1]++;
				}
				else
				{
					// Otherwise, go up again
					inserter.depth--;
				}
			}
		}
	}


	void Octree::RunBuildThreadPool(f64& totalCoeffError_)
	{
		// Create structs
		BuildThreadPool::InitialData threadData;
		threadData.config           = config;
		threadData.inputQueue       = new (std::queue<BuildThreadPool::Input>);
		threadData.inputQueueMutex  = new (std::mutex);
		threadData.outputQueue      = new (std::queue<BuildThreadPool::Output>);
		threadData.outputQueueMutex = new (std::mutex);
		threadData.shutdownAtom     = new (std::atomic<bool>);
		threadData.shutdownAtom->store(false);

		// Start threads
		BuildThreadPool threadPool;
		threadPool.StartThreads(&threadData, this);
		
		// Thread IO
		bool finished = false;
		while (!finished)
		{
			// Check if we can finish
			finished = totalCoeffError_ < config.targetErrorThreshold || nodeQueue.empty();
			if (finished)
			{
				// Stop threads in loop so no need to duplicate outputQueue code outside loop
				threadData.shutdownAtom->store(true);
				threadPool.StopThreads();
			}

			// Ensure that the threads have enough work
			threadData.inputQueueMutex->lock();
			{
				usize curQueueSize = threadData.inputQueue->size();
				for (usize i = curQueueSize; i < 10 * config.threadCount && !finished; ++i)
				{
					// Pop node off and send all required info in Input struct
					std::pair<usize, f64> nodeIdxAndErr = nodeQueue.top();
					nodeQueue.pop();

					BuildThreadPool::Input newInput;
					newInput.nodeIdxAndErr = nodeIdxAndErr;
					newInput.node          = nodes[nodeIdxAndErr.first];
					threadData.inputQueue->push(newInput);
				}
			}
			threadData.inputQueueMutex->unlock();

			// Receive completed jobs from threads
			threadData.outputQueueMutex->lock();
			{
				while (!threadData.outputQueue->empty())
				{
					BuildThreadPool::Output newOutput = threadData.outputQueue->front();
					threadData.outputQueue->pop();

					usize newNodeIdx = -1;
					switch (newOutput.basisType)
					{
						case BuildThreadPool::Output::BasisType::P:
						{
							newNodeIdx = newOutput.nodeIdx;
							free(nodes[newNodeIdx].basis.coeffs);
							totalCoeffError_ += (newOutput.newErr - newOutput.initialErr);

							break;
						}

						case BuildThreadPool::Output::BasisType::H:
						{
							// Split node if haven't already done so
							if (nodes[newOutput.nodeIdx].childIdx == -1)
							{
								free(nodes[newOutput.nodeIdx].basis.coeffs);
								nodes[newOutput.nodeIdx].basis.degree = (BASIS_MAX_DEGREE + 1);

								Subdivide(newOutput.nodeIdx);

								totalCoeffError_ -= newOutput.initialErr;
							}

							newNodeIdx        = nodes[newOutput.nodeIdx].childIdx + newOutput.childIdx;
							totalCoeffError_ += newOutput.newErr;

							break;
						}

						default:
							assert(0);
					}

					// Copy new basis. Threads don't dealloc so just take ptr
					nodes[newNodeIdx].basis = newOutput.nodeBasis;

					// Place node back into nodes queue
					const std::pair<usize, f64> newNodeIdxAndErr = { newNodeIdx,  newOutput.newErr };
					nodeQueue.push(newNodeIdxAndErr);
				}
			}
			threadData.outputQueueMutex->unlock();

			std::this_thread::sleep_for(std::chrono::milliseconds(5));
		}

		delete (threadData.inputQueue);
		delete (threadData.inputQueueMutex);
		delete (threadData.outputQueue);
		delete (threadData.outputQueueMutex);
		delete (threadData.shutdownAtom);
	}


	void Octree::Create(const Config& config_, std::function<f64(const Eigen::Vector3d& pt_)> F_)
	{
		// Reset state
		Clear();

		// Sanity
		config_.IsValid();
		config = config_;

        // Setup internal domain to be over unit cube
        configRootCentre                 = config.root.center().cast<f64>();
        configRootInvSizes               = config.root.sizes().cwiseInverse().cast<f64>();
        const Eigen::Vector3d rootBounds = config.root.sizes().cast<f64>();
        F = [F_, centre = configRootCentre, rootBounds](const Eigen::Vector3d& pt_) -> f64
        {
            return F_(pt_.cwiseProduct(rootBounds) + centre);
        };
		
		// Create root and then obtain a coarse F approximation
		f64 totalCoeffError = 0.0;
		CreateRoot();
		UniformlyRefine(totalCoeffError);

		// Main algo loop
		if (totalCoeffError > config.targetErrorThreshold)
		{
			RunBuildThreadPool(totalCoeffError);
		}

		// Optimise mem layout
		const usize nCoeffs = ReallocCoeffs();

		// Perform continuity post-process
		if (config.continuity.enforce)
		{
			PerformContinuityPostProcess(nCoeffs);
		}

		// Cleanup
		procMap.clear();
		while (!nodeQueue.empty())
		{
			nodeQueue.pop();
		}
	}


    void Octree::UnionSDF(std::function<f64(const Eigen::Vector3d& pt_)> F_)
    {
        /*
            Note: Requires keeping old tree around until operation has completed, but a lot simpler
            than doing operation over each leaf node and then recursively merging redundant leaves.

            Another caveat is that continuity corrections will change areas unaffected by the 
            union operation even if this isn't needed. This will have the effect of moving our approximation
            further and further away from the original SDF in these areas. Could be fixed by flagging
            nodes that have/haven't changed and removing shared faces of the latter from the optimisation.
        */
        Octree oldTree = std::move(*this);
        
        auto UnionF = [&oldTree, F_](const Eigen::Vector3d& pt_) -> f64
        {
            return std::min(oldTree.Query(pt_), F_(pt_));
        };

        Create(oldTree.config, UnionF);
    }


    void Octree::SubtractSDF(std::function<f64(const Eigen::Vector3d& pt_)> F_)
    {
        Octree oldTree = std::move(*this);

        auto SubtractF = [&oldTree, F_](const Eigen::Vector3d& pt_) -> f64
        {
            return std::max(oldTree.Query(pt_) * -1.0, F_(pt_));
        };

        Create(oldTree.config, SubtractF);
    }


    void Octree::IntersectSDF(std::function<f64(const Eigen::Vector3d& pt_)> F_)
    {
        Octree oldTree = std::move(*this);

        auto IntersectSDF = [&oldTree, F_](const Eigen::Vector3d& pt_) -> f64
        {
            return std::max(oldTree.Query(pt_), F_(pt_));
        };

        Create(oldTree.config, IntersectSDF);
    }


	void Octree::FromMemoryBlock(MemoryBlock octBlock_)
	{
		assert(octBlock_.size && octBlock_.ptr);
		Clear();

		const usize nCoeffs = *((usize*)octBlock_.ptr);
		coeffStore          = (f64*)malloc(sizeof(f64) * nCoeffs);
		memcpy(coeffStore, (u8*)octBlock_.ptr + sizeof(usize), sizeof(f64) * nCoeffs);

		const usize nNodes = *((usize*)((u8*)octBlock_.ptr + sizeof(usize) + sizeof(f64) * nCoeffs));
		nodes.resize(nNodes);
		memcpy(nodes.data(), (u8*)octBlock_.ptr + sizeof(usize) + sizeof(f64) * nCoeffs + sizeof(usize), sizeof(Node) * nNodes);

		config = *(Config*)((u8*)octBlock_.ptr + (sizeof(usize) + sizeof(f64) * nCoeffs + sizeof(usize) + sizeof(Node) * nNodes));

        // Recalculate
        configRootCentre   = config.root.center().cast<f64>();
        configRootInvSizes = config.root.sizes().cwiseInverse().cast<f64>();
	}


	MemoryBlock Octree::ToMemoryBlock() const
	{
		MemoryBlock octBlock = { 0, nullptr };

		usize nCoeffs = 0;
		for (usize i = 0; i < nodes.size(); ++i)
		{
			if (nodes[i].basis.degree != (BASIS_MAX_DEGREE + 1))
			{
				nCoeffs += LegendreCoeffientCount[nodes[i].basis.degree];
			}
		}

		// Calculate block size
		octBlock.size += sizeof(usize);
		octBlock.size += sizeof(f64) * nCoeffs;
		octBlock.size += sizeof(usize);
		octBlock.size += sizeof(Node) * nodes.size();
		octBlock.size += sizeof(Config);

		// Alloc block and fill
		octBlock.ptr = malloc(octBlock.size);

		*((usize*)octBlock.ptr) = nCoeffs;
		memcpy((u8*)octBlock.ptr + sizeof(usize), coeffStore, sizeof(f64) * nCoeffs);

		*(usize*)((u8*)octBlock.ptr + (sizeof(usize) + sizeof(f64) * nCoeffs)) = nodes.size();
		memcpy((u8*)octBlock.ptr + (sizeof(usize) + sizeof(f64) * nCoeffs + sizeof(usize)), nodes.data(), sizeof(Node) * nodes.size());

		*(Config*)((u8*)octBlock.ptr + (sizeof(usize) + sizeof(f64) * nCoeffs + sizeof(usize) + sizeof(Node) * nodes.size())) = config;

		return octBlock;
	}


	void Octree::Clear()
	{
		if (coeffStore)
		{
			free(coeffStore);
		}

		nodes.clear();
		while (!nodeQueue.empty())
		{
			nodeQueue.pop();
		}
	}


	usize Octree::ReallocCoeffs()
	{
		// Determine size of block required
		usize nCoeffs = 0;
		for (usize i = 0; i < nodes.size(); ++i)
		{
			if (nodes[i].basis.degree != (BASIS_MAX_DEGREE + 1))
			{
				nCoeffs += LegendreCoeffientCount[nodes[i].basis.degree];
			}
		}

		// Alloc
		if (coeffStore)
		{
			free(coeffStore);
		}
		coeffStore = (f64*)malloc(sizeof(f64) * nCoeffs);

		// Point all coeffs to new store
		usize curCoeffIdx = 0;
		for (usize i = 0; i < nodes.size(); ++i)
		{
			if (nodes[i].basis.degree != (BASIS_MAX_DEGREE + 1))
			{
				// Move mem
				memcpy(coeffStore + curCoeffIdx, nodes[i].basis.coeffs, sizeof(f64) * LegendreCoeffientCount[nodes[i].basis.degree]);
				free(nodes[i].basis.coeffs);

				// Change address
				nodes[i].basis.coeffsStart = curCoeffIdx;
				curCoeffIdx += LegendreCoeffientCount[nodes[i].basis.degree];
			}
		}

		return nCoeffs;
	}


	void Octree::TickBuildThread(const BuildThreadPool::InitialData* const threadData_)
	{
		while (1)
		{
			// Check if we need to exit
			if (threadData_->shutdownAtom->load())
			{
				break;
			}

			// Find new job
			BuildThreadPool::Input newJob;
			bool haveJob = false;
			threadData_->inputQueueMutex->lock();
			{
				haveJob = !threadData_->inputQueue->empty();
				if (haveJob)
				{
					newJob = threadData_->inputQueue->front();
					threadData_->inputQueue->pop();
				}
			}
			threadData_->inputQueueMutex->unlock();

			// Nothing to do
			if (!haveJob)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(5));
				continue;
			}

			// Estimate improvement via H and P-refinement
			Node::Basis hBases[8];
			f64 hBasesErrors[8];
			Node::Basis pBasis;
			f64 pBasisError;
			const f64 hImprovement = EstimateHImprovement(newJob, hBases, hBasesErrors);
			const f64 pImprovement = EstimatePImprovement(newJob, pBasis, pBasisError);

			// Decide whether to refine in H or P
			const usize nodeDepth   = newJob.node.depth;
			const usize basisDegree = newJob.node.basis.degree;
			const bool refineP      = basisDegree < BASIS_MAX_DEGREE && (nodeDepth == TREE_MAX_DEPTH || pImprovement > hImprovement);
			const bool refineH      = nodeDepth < TREE_MAX_DEPTH && !refineP;

			// Push to output queue
			threadData_->outputQueueMutex->lock();
			{
				if (refineP)
				{
					for (u8 i = 0; i < 8; ++i)
					{
						free(hBases[i].coeffs);
					}

					BuildThreadPool::Output newOutput;
					newOutput.initialErr = newJob.nodeIdxAndErr.second;
					newOutput.nodeIdx    = newJob.nodeIdxAndErr.first;
					newOutput.basisType  = BuildThreadPool::Output::BasisType::P;
					newOutput.newErr     = pBasisError;
					newOutput.nodeBasis  = pBasis;

					threadData_->outputQueue->push(newOutput);
				}

				if (refineH)
				{
					free(pBasis.coeffs);

					for (u8 i = 0; i < 8; ++i)
					{
						BuildThreadPool::Output newOutput;
						newOutput.initialErr = newJob.nodeIdxAndErr.second;
						newOutput.nodeIdx    = newJob.nodeIdxAndErr.first;
						newOutput.basisType  = BuildThreadPool::Output::BasisType::H;
						newOutput.childIdx   = i;
						newOutput.newErr     = hBasesErrors[i];
						newOutput.nodeBasis  = hBases[i];

						threadData_->outputQueue->push(newOutput);
					}
				}
			}
			threadData_->outputQueueMutex->unlock();
		}
	}


	f64 Octree::Query(const Eigen::Vector3d& pt_) const
	{
        // Move to unit cube
        const Eigen::Vector3d pt = (pt_ - configRootCentre).cwiseProduct(configRootInvSizes);

		// Not in volume
		if (!nodes[0].aabb.contains(pt.cast<f32>()))
		{
			return std::numeric_limits<f64>::max();
		}

		// Start at root
		usize curNodeIdx = 0;
		while (1)
		{
			const Eigen::Vector3f& aabbMin = nodes[curNodeIdx].aabb.min();
			const Eigen::Vector3f& aabbMax = nodes[curNodeIdx].aabb.max();
			const f32 curNodeAABBHalf      = (aabbMax.x() - aabbMin.x()) * 0.5f;

			const usize xIdx = (pt.x() >= (aabbMin.x() + curNodeAABBHalf));
			const usize yIdx = (pt.y() >= (aabbMin.y() + curNodeAABBHalf)) << 1;
			const usize zIdx = (pt.z() >= (aabbMin.z() + curNodeAABBHalf)) << 2;

			const usize childIdx = nodes[curNodeIdx].childIdx + xIdx + yIdx + zIdx;
			const Node& curChild = nodes[childIdx];
			if (curChild.basis.degree != (BASIS_MAX_DEGREE + 1))
			{
				// Leaf
				Node::Basis basis;
				basis.degree = curChild.basis.degree;
				basis.coeffs = coeffStore + curChild.basis.coeffsStart;

				return FApprox(basis, curChild.aabb, pt, curChild.depth);
			}
			else
			{
				// Not leaf
				curNodeIdx = childIdx;
			}
		}
	}


    bool Octree::QueryRay(const Ray& ray_, const f64 tMax_, f64& t_) const
    {
        constexpr usize MAX_STEPS = 200;
        constexpr f64 eps         = 0.0001;
        constexpr f64 minStep     = 0.0001;

        // Move to unit cube
        Ray ray((ray_.origin - configRootCentre).cwiseProduct(configRootInvSizes), ray_.direction);

        Eigen::Vector3d intMin = ray.origin;
        Eigen::Vector3d intMax;

        // Either misses root or we find first intersection point and clamp to that
        if (!nodes[0].aabb.contains(ray.origin.cast<f32>()) && !ray.IntersectAABB(nodes[0].aabb.cast<f64>(), intMin, intMax))
        {
            return false;
        }

        f64 d = 0.0;

        for (usize i = 0; i < MAX_STEPS; ++i)
        {
            const f64 v = Query(intMin + d * ray.direction);

            if (v < eps)
            {
                t_ = v;

                return true;
            }

            // Assume max of 5% err in approximated SDF and move by at least minStep 
            d += v * 0.95 + minStep;

            if (d > tMax_)
            {
                return false;
            }
        }

        return false;
    }
	
	
	f64 Octree::QueryWithGradient(const Eigen::Vector3d& pt_, Eigen::Vector3d& unitGradient_) const
	{
        // Move to unit cube
        const Eigen::Vector3d pt = (pt_ - configRootCentre).cwiseProduct(configRootInvSizes);

		// Not in volume
		if (!nodes[0].aabb.contains(pt.cast<f32>()))
		{
			return std::numeric_limits<f64>::max();
		}

		// Start at root
		usize curNodeIdx = 0;
		while (1)
		{
			const Eigen::Vector3f& aabbMin = nodes[curNodeIdx].aabb.min();
			const Eigen::Vector3f& aabbMax = nodes[curNodeIdx].aabb.max();
			const f32 curNodeAABBHalf      = (aabbMax.x() - aabbMin.x()) * 0.5f;

			const usize xIdx = (pt.x() >= (aabbMin.x() + curNodeAABBHalf));
			const usize yIdx = (pt.y() >= (aabbMin.y() + curNodeAABBHalf)) << 1;
			const usize zIdx = (pt.z() >= (aabbMin.z() + curNodeAABBHalf)) << 2;

			const usize childIdx = nodes[curNodeIdx].childIdx + xIdx + yIdx + zIdx;
			const Node& curChild = nodes[childIdx];
			if (curChild.basis.degree != (BASIS_MAX_DEGREE + 1))
			{
				// Leaf
				Node::Basis basis;
				basis.degree = curChild.basis.degree;
				basis.coeffs = coeffStore + curChild.basis.coeffsStart;

				return FApproxWithGradient(basis, curChild.aabb, pt, curChild.depth, unitGradient_);
			}
			else
			{
				// Not leaf
				curNodeIdx = childIdx;
			}
		}
	}


	void Octree::CreateRoot()
	{
		nodes.push_back({});

		Node& rootNode = nodes.back();
		rootNode.depth = 0;
		rootNode.aabb  = Eigen::AlignedBox3f(Eigen::Vector3f(-0.5, -0.5, -0.5), Eigen::Vector3f(0.5, 0.5, 0.5));

		Subdivide(0);
	}


	f64 Octree::EstimateHImprovement(const BuildThreadPool::Input& inputData_, Node::Basis* hBases_, f64* hBasesErrors_)
	{
		// Calculate maximum possible error over all children
		f64 maxNewErr = 0.0;
		for (usize i = 0; i < 8; ++i)
		{
			// Indicate to FitPolynomial that we need to calculate all coeffs
			hBases_[i].degree = 0;
			hBases_[i].coeffs = (f64*)malloc(sizeof(f64) * LegendreCoeffientCount[inputData_.node.basis.degree]);

			hBasesErrors_[i] = FitPolynomial(hBases_[i], CornerAABB(inputData_.node.aabb, i), inputData_.node.basis.degree, inputData_.node.depth + 1);
			maxNewErr        = std::max<f64>(maxNewErr, hBasesErrors_[i]);
		}

		// Equation (9)
		return (1.0 / (7.0 * LegendreCoeffientCount[inputData_.node.basis.degree])) * (inputData_.nodeIdxAndErr.second - 8.0 * maxNewErr);
	}


	f64 Octree::EstimatePImprovement(const BuildThreadPool::Input& inputData_, Node::Basis& pBasis_, f64& pBasisError_)
	{
		const u8 basisDegree = inputData_.node.basis.degree;
		const u8 nodeDepth   = inputData_.node.depth;

		// Copy current node basis into P
		pBasis_.coeffs = (f64*)malloc(sizeof(f64) * LegendreCoeffientCount[basisDegree + 1]);
		memcpy(pBasis_.coeffs, inputData_.node.basis.coeffs, sizeof(f64) * LegendreCoeffientCount[basisDegree]);
		pBasis_.degree = basisDegree;

		// Calculate new error and fit to P
		pBasisError_ = FitPolynomial(pBasis_, inputData_.node.aabb, basisDegree + 1, nodeDepth);

		// Equation (8)
		return (1.0 / (LegendreCoeffientCount[basisDegree + 1] - LegendreCoeffientCount[basisDegree])) * (inputData_.nodeIdxAndErr.second - 8.0 * pBasisError_);
	}


	f64 Octree::FApprox(const Node::Basis& basis_, const Eigen::AlignedBox3f& aabb_, const Eigen::Vector3d& pt_, const usize depth_) const
	{
		// Move pt_ to unit cube
		const Eigen::Vector3d unitPt = (pt_ - aabb_.center().cast<f64>()) * (2 << depth_);

        // Create lookup table for pt_
        f64 LpXLookup[BASIS_MAX_DEGREE][3];
        for (usize i = 0; i < 3; ++i)
        {
            // Constant
            LpXLookup[0][i] = NormalisedLengths[0][depth_];

            // Initial values for recurrence
            f64 LjMinus2 = 0.0;
            f64 LjMinus1 = 1.0;
            f64 Lj       = 1.0;

            // Determine remaining values
            for (usize j = 1; j <= basis_.degree; ++j)
            {
                Lj       = LegendreCoefficent[j][0] * unitPt(i) * LjMinus1 - LegendreCoefficent[j][1] * LjMinus2;
                LjMinus2 = LjMinus1;
                LjMinus1 = Lj;

                LpXLookup[j][i] = Lj * NormalisedLengths[j][depth_];
            }
        }

		// Sum up basis coeffs
		f64 fApprox = 0.0;
		for (usize i = 0; i < LegendreCoeffientCount[basis_.degree]; ++i)
		{
			f64 Lp = 1.0;
			for (usize j = 0; j < 3; ++j)
			{
				Lp *= LpXLookup[BasisIndexValues[i][j]][j];
			}

			fApprox += basis_.coeffs[i] * Lp;
		}

		return fApprox;
	}
	
	
	f64 Octree::FApproxWithGradient(const Node::Basis& basis_, const Eigen::AlignedBox3f& aabb_, const Eigen::Vector3d& pt_, const usize depth_, Eigen::Vector3d& unitGradient_) const
	{
		// Move pt_ to unit cube
		const Eigen::Vector3d unitPt = (pt_ - aabb_.center().cast<f64>()) * (2 << depth_);

		// CD eps
		const f64 eps = 0.0001;

        // Create lookup table for pt_
        f64 LpXLookup[BASIS_MAX_DEGREE][3][3];
        for (usize i = 0; i < 3; ++i)
        {
            // Constants
            LpXLookup[0][i][0] = NormalisedLengths[0][depth_];
            LpXLookup[0][i][1] = NormalisedLengths[0][depth_];
            LpXLookup[0][i][2] = NormalisedLengths[0][depth_];

            // Initial values for recurrences
            f64 Lj0Minus2 = 0.0;
            f64 Lj0Minus1 = 1.0;
            f64 Lj0       = 1.0;

            f64 Lj1Minus2 = 0.0;
            f64 Lj1Minus1 = 1.0;
            f64 Lj1       = 1.0;

            f64 Lj2Minus2 = 0.0;
            f64 Lj2Minus1 = 1.0;
            f64 Lj2       = 1.0;

            // Determine remaining values
            for (usize j = 1; j <= basis_.degree; ++j)
            {
                Lj0       = LegendreCoefficent[j][0] * unitPt(i) * Lj0Minus1 - LegendreCoefficent[j][1] * Lj0Minus2;
                Lj0Minus2 = Lj0Minus1;
                Lj0Minus1 = Lj0;

                Lj1       = LegendreCoefficent[j][0] * (unitPt(i) + eps) * Lj1Minus1 - LegendreCoefficent[j][1] * Lj1Minus2;
                Lj1Minus2 = Lj1Minus1;
                Lj1Minus1 = Lj1;

                Lj2       = LegendreCoefficent[j][0] * (unitPt(i) - eps) * Lj2Minus1 - LegendreCoefficent[j][1] * Lj2Minus2;
                Lj2Minus2 = Lj2Minus1;
                Lj2Minus1 = Lj2;

                LpXLookup[j][i][0] = Lj0 * NormalisedLengths[j][depth_];
                LpXLookup[j][i][1] = Lj1 * NormalisedLengths[j][depth_];
                LpXLookup[j][i][2] = Lj2 * NormalisedLengths[j][depth_];
            }
        }

		// Calculate gradient
		for (usize k = 0; k < 3; ++k)
		{
			f64 fApproxP1 = 0.0;
			f64 fApproxM1 = 0.0;

			for (usize i = 0; i < LegendreCoeffientCount[basis_.degree]; ++i)
			{
				fApproxP1 += basis_.coeffs[i] * LpXLookup[BasisIndexValues[i][k]][k][1];
				fApproxM1 += basis_.coeffs[i] * LpXLookup[BasisIndexValues[i][k]][k][2];
			}

			unitGradient_(k) = (fApproxP1 - fApproxM1) / (2.0 * eps);
		}
		unitGradient_.normalize();

		// Return actual value
		f64 fApprox = 0.0;
		for (usize i = 0; i < LegendreCoeffientCount[basis_.degree]; ++i)
		{
			f64 Lp = 1.0;
			for (usize j = 0; j < 3; ++j)
			{
				Lp *= LpXLookup[BasisIndexValues[i][j]][j][0];
			}

			fApprox += basis_.coeffs[i] * Lp;
		}

		return fApprox;
	}
		

	f64 Octree::LpX(const usize p_, const f64 x_)
	{
		// Set initial values for recurrence
		f64 LiMinus2 = 0.0;
		f64 LiMinus1 = 1.0;
		f64 Li       = 1.0;

		// Apply procedure p_ times
		for (usize i = 1; i <= p_; ++i)
		{
			Li       = LegendreCoefficent[i][0] * x_ * LiMinus1 - LegendreCoefficent[i][1] * LiMinus2;
			LiMinus2 = LiMinus1;
			LiMinus1 = Li;
		}

		return Li;
	}


	f64 Octree::FitPolynomial(Node::Basis& basis_, const Eigen::AlignedBox3f& aabb_, const u8 degree_, const usize depth_)
	{
		assert(basis_.degree < degree_);

		// If we already have a basis, we can re-use coeffs up to degree_ - 1
		const usize startingIdx = basis_.degree > 0 ? LegendreCoeffientCount[basis_.degree] : 0;
		const usize endingIdx   = LegendreCoeffientCount[degree_];

		// Determine GQ steps to perform
		const usize GQStart = SumToN[4 * degree_];
		const usize GQEnd   = SumToN[4 * degree_ + 1];

		// Determine mappings from aabb_ to unit cube
		const Eigen::Vector3d aabbScale  = aabb_.sizes().cast<f64>() * 0.5;
		const Eigen::Vector3d aabbCentre = aabb_.center().cast<f64>();

		// Set coeffs to 0
		memset(basis_.coeffs + startingIdx, 0, (endingIdx - startingIdx) * sizeof(f64));

		// Perform GQ and calc new basis
		for (usize i = GQStart; i < GQEnd; ++i)
		{
			for (usize j = GQStart; j < GQEnd; ++j)
			{
				for (usize k = GQStart; k < GQEnd; ++k)
				{
					// Generate sample
					const Eigen::Vector3d unitSample(LegendreRoots[i], LegendreRoots[j], LegendreRoots[k]);
					const Eigen::Vector3d unitWeights(LegendreWeights[i], LegendreWeights[j], LegendreWeights[k]);
					const Eigen::Vector3d aabbSample = unitSample.cwiseProduct(aabbScale) + aabbCentre;

					// Evaluate F at this position
					const f64 weightsMult    = unitWeights.prod();
					const f64 normScalesMult = aabbScale.prod();
					const f64 FaabSample     = normScalesMult * weightsMult * F(aabbSample);

					// Use this sample to continue calculating coeffs
					for (usize c = startingIdx; c < endingIdx; ++c)
					{
						f64 Lp = 1.0;
						for (usize p = 0; p < 3; ++p)
						{
							Lp *= LpX(BasisIndexValues[c][p], unitSample(p));
							Lp *= NormalisedLengths[BasisIndexValues[c][p]][depth_];
						}

						basis_.coeffs[c] += Lp * FaabSample;
					}
				}
			}
		}

		// Update degree
		basis_.degree = degree_;

		// Determine new error using (6)
		f64 newError = 0.0;
		for (usize i = 0; i < endingIdx; ++i)
		{
			if ((BasisIndexValues[i][0] + BasisIndexValues[i][1] + BasisIndexValues[i][2]) == degree_)
			{
				newError += basis_.coeffs[i] * basis_.coeffs[i];
			}
		}

		switch (config.nearnessWeighting.type)
		{
			case Config::NearnessWeighting::None:
			{
				return newError;
			}

			case Config::NearnessWeighting::Polynomial:
			{
				return newError * CalculatePolyWeighting(basis_, aabb_, depth_);
			}

			case Config::NearnessWeighting::Exponential:
			{
				return newError * CalculateExpWeighting(basis_, aabb_, depth_);
			}

			default:
				assert(0);
		}

        return newError;
	}


	Eigen::AlignedBox3f Octree::CornerAABB(const Eigen::AlignedBox3f& aabb_, const usize i_)
	{
		Eigen::AlignedBox3f cornerAABB = aabb_;

		if (i_ & 1)
		{
			cornerAABB.min().x() = (aabb_.max().x() + aabb_.min().x()) * 0.5f;
		}
		else
		{
			cornerAABB.max().x() = (aabb_.max().x() + aabb_.min().x()) * 0.5f;
		}

		if (i_ & 2)
		{
			cornerAABB.min().y() = (aabb_.max().y() + aabb_.min().y()) * 0.5f;
		}
		else
		{
			cornerAABB.max().y() = (aabb_.max().y() + aabb_.min().y()) * 0.5f;
		}

		if (i_ & 4)
		{
			cornerAABB.min().z() = (aabb_.max().z() + aabb_.min().z()) * 0.5f;
		}
		else
		{
			cornerAABB.max().z() = (aabb_.max().z() + aabb_.min().z()) * 0.5f;
		}


		return cornerAABB;
	}


	void Octree::Subdivide(const usize nodeIdx_)
	{
		nodes[nodeIdx_].childIdx = (u32)nodes.size();

		for (usize i = 0; i < 8; ++i)
		{
			nodes.push_back({});

			// Create AABB
			Node& newNode = nodes.back();
			newNode.aabb  = CornerAABB(nodes[nodeIdx_].aabb, i);
			newNode.depth = nodes[nodeIdx_].depth + 1;
		}
	}


#if HAS_STB
	void Octree::OutputFunctionSlice(const char* fName_, const f64 c_)
	{
		// Obtain mem and set size
		const usize nSamples = 2048;
		f64* sdfVals         = (f64*)malloc(sizeof(f64) * nSamples * nSamples);
		u8* imageBytes       = (u8*)malloc(nSamples * nSamples * 3);

		// Fill in pixels
		std::pair<f64, f64> minMaxPosVals = { std::numeric_limits<f64>::max(), 0.0 };
		std::pair<f64, f64> minMaxNegVals = { 0.0, std::numeric_limits<f64>::max() * -1.0 };
		for (usize i = 0; i < nSamples; ++i)
		{
			for (usize j = 0; j < nSamples; ++j)
			{
				// Generate sample
				Eigen::Vector3d sampleLoc = nodes[0].aabb.min().cast<f64>();
				sampleLoc.x()            += ((j + 0.5) / nSamples);
				sampleLoc.y()            += ((i + 0.5) / nSamples);
				sampleLoc.z()             = c_;

				// Update 
				const f64 sdfVal = Query(sampleLoc);
				if (sdfVal > EPSILON_F32)
				{
					minMaxPosVals.first  = std::min<f64>(sdfVal, minMaxPosVals.first);
					minMaxPosVals.second = std::max<f64>(sdfVal, minMaxPosVals.second);
				}
				else
				{
					minMaxNegVals.first  = std::min<f64>(sdfVal, minMaxNegVals.first);
					minMaxNegVals.second = std::max<f64>(sdfVal, minMaxNegVals.second);
				}

				sdfVals[i * nSamples + j] = sdfVal;
			}
		}

		// Rescale to obtain better image constrant
		for (usize i = 0; i < nSamples; ++i)
		{
			for (usize j = 0; j < nSamples; ++j)
			{
				const f32 unscaledVal = (f32)sdfVals[i * nSamples + j];
				if (unscaledVal > 0.0f)
				{
					const u8 scaledByte = (u8)(255 * (unscaledVal - minMaxPosVals.second) / (minMaxPosVals.first - minMaxPosVals.second));

					imageBytes[3 * (i * nSamples + j)]     = 0;
					imageBytes[3 * (i * nSamples + j) + 1] = scaledByte;
					imageBytes[3 * (i * nSamples + j) + 2] = 0;
				}
				else
				{
					const u8 scaledByte = (u8)(255 * (unscaledVal - minMaxNegVals.first) / (minMaxNegVals.second - minMaxNegVals.first));

					imageBytes[3 * (i * nSamples + j)]     = 0;
					imageBytes[3 * (i * nSamples + j) + 1] = 0;
					imageBytes[3 * (i * nSamples + j) + 2] = scaledByte;
				}
			}
		}

		// Write to file
		char fNameWithExt[PATH_MAX];
		sprintf_s(fNameWithExt, PATH_MAX, "%s.bmp", fName_);
		const bool success = stbi_write_bmp(fNameWithExt, nSamples, nSamples, 3, imageBytes);
		assert(success);

		// Cleanup
		free(imageBytes);
		free(sdfVals);
	}
#endif


	f64 Octree::CalculatePolyWeighting(const Node::Basis& basis_, const Eigen::AlignedBox3f& aabb_, const usize depth_)
	{
		// Determine constants
		const f64 V = aabb_.volume();
		const f64 d = sqrt(3.0);

		// Approximate F over aabb_ with MC
		const usize nSamples = 100;
		f64 fIntegral = 0.0;
		for (usize i = 0; i < nSamples; ++i)
		{
			fIntegral += FApprox(basis_, aabb_, aabb_.sample().cast<f64>(), depth_);
		}
		fIntegral *= (aabb_.volume() / nSamples);
		fIntegral  = std::abs<f64>(fIntegral);

		// Calc weighting and then clamp to [0, 1]
		const f64 k = std::pow(1.0 - fIntegral / (V * d), config.nearnessWeighting.strength);
		return std::min<f64>(1.0, std::max<f64>(k, 0.0));
	}


	f64 Octree::CalculateExpWeighting(const Node::Basis& basis_, const Eigen::AlignedBox3f& aabb_, const usize depth_)
	{
		// Determine constants
		const f64 V = aabb_.volume();
		const f64 d = sqrt(3.0);

		// Approximate F over aabb_ with MC
		const usize nSamples = 100;
		f64 fIntegral = 0.0;
		for (usize i = 0; i < nSamples; ++i)
		{
			fIntegral += FApprox(basis_, aabb_, aabb_.sample().cast<f64>(), depth_);
		}
		fIntegral *= (aabb_.volume() / nSamples);
		fIntegral  = std::abs<f64>(fIntegral);

		// Calc weighting
		return std::exp(-1.0 * config.nearnessWeighting.strength * fIntegral / (V * d));
	}


	void Octree::EvaluateSharedFaceIntegralNumerically(const usize nodeIdxA_, const usize nodeIdxB_, const u8 dim_, std::vector<Eigen::Triplet<f64>>& matTriplets_)
	{
        using StorageIndex = Eigen::SparseMatrix<f64>::StorageIndex;

		const Node& nodeA = nodes[nodeIdxA_];
		const Node& nodeB = nodes[nodeIdxB_];

		// Determine the shared face and integration constants
		const Eigen::AlignedBox3f sharedFace   = Eigen::AlignedBox3f(nodeA.aabb).clamp(nodeB.aabb);
		const Eigen::Vector3d sharedFaceScale  = sharedFace.sizes().cast<f64>() * 0.5;
		const Eigen::Vector3d sharedFaceCentre = sharedFace.center().cast<f64>();

		// Determine GQ steps to perform
		const usize maxDegree = std::max<usize>(nodeA.basis.degree, nodeB.basis.degree);
		const usize GQStart   = SumToN[maxDegree];
		const usize GQEnd     = SumToN[maxDegree + 1];

		// Determine shared face scale
		const usize depthDiff = nodeA.depth > nodeB.depth ? (nodeA.depth - nodeB.depth) : (nodeB.depth - nodeA.depth);
		const f64 invDist = 1.0 / std::pow<f64>(2.0, depthDiff);

		// Determine shared face translation
		Eigen::Vector3d invTranslation(Eigen::Vector3d::Zero());
		if (nodeA.depth > nodeB.depth)
		{
			invTranslation(Mod3Lookup[dim_][1]) = (nodeA.aabb.center()(Mod3Lookup[dim_][1]) - nodeB.aabb.center()(Mod3Lookup[dim_][1])) / (nodeA.aabb.sizes()(Mod3Lookup[dim_][1]) * 0.5);
			invTranslation(Mod3Lookup[dim_][2]) = (nodeA.aabb.center()(Mod3Lookup[dim_][2]) - nodeB.aabb.center()(Mod3Lookup[dim_][2])) / (nodeA.aabb.sizes()(Mod3Lookup[dim_][2]) * 0.5);
		}
		else
		{
			invTranslation(Mod3Lookup[dim_][1]) = (nodeB.aabb.center()(Mod3Lookup[dim_][1]) - nodeA.aabb.center()(Mod3Lookup[dim_][1])) / (nodeB.aabb.sizes()(Mod3Lookup[dim_][1]) * 0.5);
			invTranslation(Mod3Lookup[dim_][2]) = (nodeB.aabb.center()(Mod3Lookup[dim_][2]) - nodeA.aabb.center()(Mod3Lookup[dim_][2])) / (nodeB.aabb.sizes()(Mod3Lookup[dim_][2]) * 0.5);
		}
		invTranslation *= invDist;

		// Pl * Pl in (13)
		for (usize i = 0; i < LegendreCoeffientCount[nodeA.basis.degree]; ++i)
		{
			for (usize j = 0; j < LegendreCoeffientCount[nodeA.basis.degree]; ++j)
			{
				// GL quadrature as in FitPolynomial
				f64 integral = 0.0;
				for (usize x = GQStart; x < GQEnd; ++x)
				{
					for (usize y = GQStart; y < GQEnd; ++y)
					{
						// Determine sample on shared face
						Eigen::Vector3d aUnitSample;
						aUnitSample(dim_)                = 1.0;
						aUnitSample(Mod3Lookup[dim_][1]) = LegendreRoots[x];
						aUnitSample(Mod3Lookup[dim_][2]) = LegendreRoots[y];

						// Move to smaller shared face
						if (nodeB.depth > nodeA.depth)
						{
							aUnitSample(Mod3Lookup[dim_][1]) = aUnitSample(Mod3Lookup[dim_][1]) * invDist + invTranslation(Mod3Lookup[dim_][1]);
							aUnitSample(Mod3Lookup[dim_][2]) = aUnitSample(Mod3Lookup[dim_][2]) * invDist + invTranslation(Mod3Lookup[dim_][2]);
						}

						// Calculate integral step
						f64 areaVal = LegendreWeights[x] * LegendreWeights[y];
						for (usize k = 0; k < 3; ++k)
						{
							areaVal *= LpX(BasisIndexValues[i][k], aUnitSample(k));
							areaVal *= LpX(BasisIndexValues[j][k], aUnitSample(k));
						}
						integral += areaVal;
					}
				}

				f64 basisWeights = 1.0;
				for (usize k = 0; k < 3; ++k)
				{
					basisWeights *= NormalisedLengths[BasisIndexValues[i][k]][nodeA.depth];
					basisWeights *= NormalisedLengths[BasisIndexValues[j][k]][nodeA.depth];
				}
				integral *= sharedFaceScale(Mod3Lookup[dim_][1]) * sharedFaceScale(Mod3Lookup[dim_][2]) * basisWeights;

				// Finally...
				if (std::abs<f32>((f32)integral) > EPSILON_F32)
				{
					matTriplets_.push_back(Eigen::Triplet<f64>((StorageIndex)(nodeA.basis.coeffsStart + i), 
                                                               (StorageIndex)(nodeA.basis.coeffsStart + j), integral));
				}
			}
		}

		// -Pl * Pr and -Pr * Pl in (13)
		for (usize i = 0; i < LegendreCoeffientCount[nodeA.basis.degree]; ++i)
		{
			for (usize j = 0; j < LegendreCoeffientCount[nodeB.basis.degree]; ++j)
			{
				// GL quadrature as in FitPolynomial
				f64 integral = 0.0;
				for (usize x = GQStart; x < GQEnd; ++x)
				{
					for (usize y = GQStart; y < GQEnd; ++y)
					{
						// Determine samples on shared face
						Eigen::Vector3d aUnitSample;
						aUnitSample(dim_)                = 1.0;
						aUnitSample(Mod3Lookup[dim_][1]) = LegendreRoots[x];
						aUnitSample(Mod3Lookup[dim_][2]) = LegendreRoots[y];

						Eigen::Vector3d bUnitSample;
						bUnitSample(dim_)                = -1.0;
						bUnitSample(Mod3Lookup[dim_][1]) = LegendreRoots[x];
						bUnitSample(Mod3Lookup[dim_][2]) = LegendreRoots[y];

						// Move to smaller shared face
						if (nodeB.depth > nodeA.depth)
						{
							aUnitSample(Mod3Lookup[dim_][1]) = aUnitSample(Mod3Lookup[dim_][1]) * invDist + invTranslation(Mod3Lookup[dim_][1]);
							aUnitSample(Mod3Lookup[dim_][2]) = aUnitSample(Mod3Lookup[dim_][2]) * invDist + invTranslation(Mod3Lookup[dim_][2]);
						}
						else if (nodeA.depth > nodeB.depth)
						{
							bUnitSample(Mod3Lookup[dim_][1]) = bUnitSample(Mod3Lookup[dim_][1]) * invDist + invTranslation(Mod3Lookup[dim_][1]);
							bUnitSample(Mod3Lookup[dim_][2]) = bUnitSample(Mod3Lookup[dim_][2]) * invDist + invTranslation(Mod3Lookup[dim_][2]);
						}

						// Calculate integral step
						f64 areaVal = LegendreWeights[x] * LegendreWeights[y];
						for (usize k = 0; k < 3; ++k)
						{
							areaVal *= LpX(BasisIndexValues[i][k], aUnitSample(k));
							areaVal *= LpX(BasisIndexValues[j][k], bUnitSample(k));
						}
						integral += areaVal;
					}
				}

				f64 basisWeights = 1.0;
				for (usize k = 0; k < 3; ++k)
				{
					basisWeights *= NormalisedLengths[BasisIndexValues[i][k]][nodeA.depth];
					basisWeights *= NormalisedLengths[BasisIndexValues[j][k]][nodeB.depth];
				}
				integral *= sharedFaceScale(Mod3Lookup[dim_][1]) * sharedFaceScale(Mod3Lookup[dim_][2]) * basisWeights * -1.0;

				if (std::abs<f32>((f32)integral) > EPSILON_F32)
				{
					matTriplets_.push_back(Eigen::Triplet<f64>((StorageIndex)(nodeA.basis.coeffsStart + i),
                                                               (StorageIndex)(nodeB.basis.coeffsStart + j), integral));
					matTriplets_.push_back(Eigen::Triplet<f64>((StorageIndex)(nodeB.basis.coeffsStart + j), 
                                                               (StorageIndex)(nodeA.basis.coeffsStart + i), integral));
				}
			}
		}

		// Pr * Pr in (13)
		for (usize i = 0; i < LegendreCoeffientCount[nodeB.basis.degree]; ++i)
		{
			for (usize j = 0; j < LegendreCoeffientCount[nodeB.basis.degree]; ++j)
			{
				// GL quadrature as in FitPolynomial
				f64 integral = 0.0;
				for (usize x = GQStart; x < GQEnd; ++x)
				{
					for (usize y = GQStart; y < GQEnd; ++y)
					{
						// Determine sample on shared face
						Eigen::Vector3d bUnitSample;
						bUnitSample(dim_)                = -1.0;
						bUnitSample(Mod3Lookup[dim_][1]) = LegendreRoots[x];
						bUnitSample(Mod3Lookup[dim_][2]) = LegendreRoots[y];

						if (nodeA.depth > nodeB.depth)
						{
							bUnitSample(Mod3Lookup[dim_][1]) = bUnitSample(Mod3Lookup[dim_][1]) * invDist + invTranslation(Mod3Lookup[dim_][1]);
							bUnitSample(Mod3Lookup[dim_][2]) = bUnitSample(Mod3Lookup[dim_][2]) * invDist + invTranslation(Mod3Lookup[dim_][2]);
						}

						// Calculate integral step
						f64 areaVal = LegendreWeights[x] * LegendreWeights[y];
						for (usize k = 0; k < 3; ++k)
						{
							areaVal *= LpX(BasisIndexValues[i][k], bUnitSample(k));
							areaVal *= LpX(BasisIndexValues[j][k], bUnitSample(k));
						}
						integral += areaVal;
					}
				}

				f64 basisWeights = 1.0;
				for (usize k = 0; k < 3; ++k)
				{
					basisWeights *= NormalisedLengths[BasisIndexValues[i][k]][nodeB.depth];
					basisWeights *= NormalisedLengths[BasisIndexValues[j][k]][nodeB.depth];
				}
				integral *= sharedFaceScale(Mod3Lookup[dim_][1]) * sharedFaceScale(Mod3Lookup[dim_][2]) * basisWeights;
				
				if (std::abs<f32>((f32)integral) > EPSILON_F32)
				{
					matTriplets_.push_back(Eigen::Triplet<f64>((StorageIndex)(nodeB.basis.coeffsStart + i), 
                                                               (StorageIndex)(nodeB.basis.coeffsStart + j), integral));
				}
			}
		}
	}


	void Octree::EvaluateSharedFaceIntegralAnalytically(const usize nodeIdxA_, const usize nodeIdxB_, const u8 dim_, std::vector<Eigen::Triplet<f64>>& matTriplets_)
	{
        using StorageIndex = Eigen::SparseMatrix<f64>::StorageIndex;

		const Node& nodeA = nodes[nodeIdxA_];
		const Node& nodeB = nodes[nodeIdxB_];

		// Pl * Pl in (13)
		for (usize i = 0; i < LegendreCoeffientCount[nodeA.basis.degree]; ++i)
		{
			for (usize j = 0; j < LegendreCoeffientCount[nodeA.basis.degree]; ++j)
			{
				const bool kdIJ1 = (BasisIndexValues[i][Mod3Lookup[dim_][1]] == BasisIndexValues[j][Mod3Lookup[dim_][1]]);
				const bool kdIJ2 = (BasisIndexValues[i][Mod3Lookup[dim_][2]] == BasisIndexValues[j][Mod3Lookup[dim_][2]]);
				if (!kdIJ1 || !kdIJ2)
				{
					// Coeff in block = 0
					continue;
				}

				f64 integral = 1.0;
				integral    *= LpX(BasisIndexValues[i][dim_], 1.0);
				integral    *= NormalisedLengths[BasisIndexValues[i][dim_]][nodeA.depth];
				integral    *= LpX(BasisIndexValues[j][dim_], 1.0);
				integral    *= NormalisedLengths[BasisIndexValues[j][dim_]][nodeA.depth];

                matTriplets_.push_back(Eigen::Triplet<f64>((StorageIndex)(nodeA.basis.coeffsStart + i),
                                                           (StorageIndex)(nodeA.basis.coeffsStart + j), integral));
			}
		}

		// -Pl * Pr and -Pr * Pl in (13)
		for (usize i = 0; i < LegendreCoeffientCount[nodeA.basis.degree]; ++i)
		{
			for (usize j = 0; j < LegendreCoeffientCount[nodeB.basis.degree]; ++j)
			{
				const bool kdIJ1 = (BasisIndexValues[i][Mod3Lookup[dim_][1]] == BasisIndexValues[j][Mod3Lookup[dim_][1]]);
				const bool kdIJ2 = (BasisIndexValues[i][Mod3Lookup[dim_][2]] == BasisIndexValues[j][Mod3Lookup[dim_][2]]);
				if (!kdIJ1 || !kdIJ2)
				{
					// Coeff in block = 0
					continue;
				}

				f64 integral = -1.0;
				integral    *= LpX(BasisIndexValues[i][dim_], 1.0);
				integral    *= NormalisedLengths[BasisIndexValues[i][dim_]][nodeA.depth];
				integral    *= LpX(BasisIndexValues[j][dim_], -1.0);
				integral    *= NormalisedLengths[BasisIndexValues[j][dim_]][nodeB.depth];

                matTriplets_.push_back(Eigen::Triplet<f64>((StorageIndex)(nodeA.basis.coeffsStart + i),
                                                           (StorageIndex)(nodeB.basis.coeffsStart + j), integral));
                matTriplets_.push_back(Eigen::Triplet<f64>((StorageIndex)(nodeB.basis.coeffsStart + j),
                                                           (StorageIndex)(nodeA.basis.coeffsStart + i), integral));
			}
		}

		// Pr * Pr in (13)
		for (usize i = 0; i < LegendreCoeffientCount[nodeB.basis.degree]; ++i)
		{
			for (usize j = 0; j < LegendreCoeffientCount[nodeB.basis.degree]; ++j)
			{
				const bool kdIJ1 = (BasisIndexValues[i][Mod3Lookup[dim_][1]] == BasisIndexValues[j][Mod3Lookup[dim_][1]]);
				const bool kdIJ2 = (BasisIndexValues[i][Mod3Lookup[dim_][2]] == BasisIndexValues[j][Mod3Lookup[dim_][2]]);
				if (!kdIJ1 || !kdIJ2)
				{
					// Coeff in block = 0
					continue;
				}

				f64 integral = 1.0;
				integral    *= LpX(BasisIndexValues[i][dim_], -1.0);
				integral    *= NormalisedLengths[BasisIndexValues[i][dim_]][nodeB.depth];
				integral    *= LpX(BasisIndexValues[j][dim_], -1.0);
				integral    *= NormalisedLengths[BasisIndexValues[j][dim_]][nodeB.depth];

                matTriplets_.push_back(Eigen::Triplet<f64>((StorageIndex)(nodeB.basis.coeffsStart + i),
                                                           (StorageIndex)(nodeB.basis.coeffsStart + j), integral));
			}
		}
	}


	void Octree::NodeProc(const usize nodeIdx_, std::queue<ContinuityThreadPool::Input>& jobQueue_)
	{
		const Node& node = nodes[nodeIdx_];
		const bool hasChildren = node.childIdx != -1;

		if (hasChildren)
		{
			// Call NodeProc on all the children
			for (usize i = 0; i < 8; ++i)
			{
				NodeProc(node.childIdx + i, jobQueue_);
			}

			// Call FaceProc[X/Y/Z] on shared faces
			FaceProcX(node.childIdx + 0, node.childIdx + 1, jobQueue_);
			FaceProcX(node.childIdx + 2, node.childIdx + 3, jobQueue_);
			FaceProcX(node.childIdx + 4, node.childIdx + 5, jobQueue_);
			FaceProcX(node.childIdx + 6, node.childIdx + 7, jobQueue_);

			FaceProcY(node.childIdx + 0, node.childIdx + 2, jobQueue_);
			FaceProcY(node.childIdx + 1, node.childIdx + 3, jobQueue_);
			FaceProcY(node.childIdx + 4, node.childIdx + 6, jobQueue_);
			FaceProcY(node.childIdx + 5, node.childIdx + 7, jobQueue_);

			FaceProcZ(node.childIdx + 0, node.childIdx + 4, jobQueue_);
			FaceProcZ(node.childIdx + 1, node.childIdx + 5, jobQueue_);
			FaceProcZ(node.childIdx + 2, node.childIdx + 6, jobQueue_);
			FaceProcZ(node.childIdx + 3, node.childIdx + 7, jobQueue_);
		}
	}


	void Octree::FaceProcX(const usize nodeIdxA_, const usize nodeIdxB_, std::queue<ContinuityThreadPool::Input>& jobQueue_)
	{
		const Node& nodeA       = nodes[nodeIdxA_];
		const Node& nodeB       = nodes[nodeIdxB_];
		const bool aHasChildren = nodeA.childIdx != -1;
		const bool bHasChildren = nodeB.childIdx != -1;

		// Continue down if we aren't yet at leaves. Otherwise, calculate the integral between the leaf faces
		if (aHasChildren || bHasChildren)
		{
			FaceProcX(aHasChildren ? nodeA.childIdx + 1: nodeIdxA_, bHasChildren ? nodeB.childIdx + 0: nodeIdxB_, jobQueue_);
			FaceProcX(aHasChildren ? nodeA.childIdx + 3: nodeIdxA_, bHasChildren ? nodeB.childIdx + 2: nodeIdxB_, jobQueue_);
			FaceProcX(aHasChildren ? nodeA.childIdx + 5: nodeIdxA_, bHasChildren ? nodeB.childIdx + 4: nodeIdxB_, jobQueue_);
			FaceProcX(aHasChildren ? nodeA.childIdx + 7: nodeIdxA_, bHasChildren ? nodeB.childIdx + 6: nodeIdxB_, jobQueue_);
		}
		else
		{
			// Order nodes based on dim
			const usize nodeIdxA = nodes[nodeIdxA_].aabb.min()(0) < nodes[nodeIdxB_].aabb.min()(0) ? nodeIdxA_ : nodeIdxB_;
			const usize nodeIdxB = nodes[nodeIdxA_].aabb.min()(0) < nodes[nodeIdxB_].aabb.min()(0) ? nodeIdxB_ : nodeIdxA_;

			// Determine if we've done this integral before
			if (procMap.find(std::pair<usize, usize>(nodeIdxA, nodeIdxB)) != procMap.end())
			{
				return;
			}
			else
			{
				procMap.emplace(std::pair<usize, usize>(nodeIdxA, nodeIdxB), true);
			}

			// New job is unique
			ContinuityThreadPool::Input newJob;
			newJob.nodeIdxs = std::pair<usize, usize>(nodeIdxA, nodeIdxB);
			newJob.dim      = 0;
			jobQueue_.push(newJob);
		}
	}


	void Octree::FaceProcY(const usize nodeIdxA_, const usize nodeIdxB_, std::queue<ContinuityThreadPool::Input>& jobQueue_)
	{
		const Node& nodeA       = nodes[nodeIdxA_];
		const Node& nodeB       = nodes[nodeIdxB_];
		const bool aHasChildren = nodeA.childIdx != -1;
		const bool bHasChildren = nodeB.childIdx != -1;

		// Continue down if we aren't yet at leaves. Otherwise, calculate the integral between the leaf faces
		if (aHasChildren || bHasChildren)
		{
			FaceProcY(aHasChildren ? nodeA.childIdx + 2 : nodeIdxA_, bHasChildren ? nodeB.childIdx + 0 : nodeIdxB_, jobQueue_);
			FaceProcY(aHasChildren ? nodeA.childIdx + 3 : nodeIdxA_, bHasChildren ? nodeB.childIdx + 1 : nodeIdxB_, jobQueue_);
			FaceProcY(aHasChildren ? nodeA.childIdx + 6 : nodeIdxA_, bHasChildren ? nodeB.childIdx + 4 : nodeIdxB_, jobQueue_);
			FaceProcY(aHasChildren ? nodeA.childIdx + 7 : nodeIdxA_, bHasChildren ? nodeB.childIdx + 5 : nodeIdxB_, jobQueue_);
		}
		else
		{
			// Order nodes based on dim
			const usize nodeIdxA = nodes[nodeIdxA_].aabb.min()(1) < nodes[nodeIdxB_].aabb.min()(1) ? nodeIdxA_ : nodeIdxB_;
			const usize nodeIdxB = nodes[nodeIdxA_].aabb.min()(1) < nodes[nodeIdxB_].aabb.min()(1) ? nodeIdxB_ : nodeIdxA_;

			// Determine if we've done this integral before
			if (procMap.find(std::pair<usize, usize>(nodeIdxA, nodeIdxB)) != procMap.end())
			{
				return;
			}
			else
			{
				procMap.emplace(std::pair<usize, usize>(nodeIdxA, nodeIdxB), true);
			}

			// New job is unique
			ContinuityThreadPool::Input newJob;
			newJob.nodeIdxs = std::pair<usize, usize>(nodeIdxA, nodeIdxB);
			newJob.dim      = 1;
			jobQueue_.push(newJob);
		}
	}


	void Octree::FaceProcZ(const usize nodeIdxA_, const usize nodeIdxB_, std::queue<ContinuityThreadPool::Input>& jobQueue_)
	{
		const Node& nodeA       = nodes[nodeIdxA_];
		const Node& nodeB       = nodes[nodeIdxB_];
		const bool aHasChildren = nodeA.childIdx != -1;
		const bool bHasChildren = nodeB.childIdx != -1;

		// Continue down if we aren't yet at leaves. Otherwise, calculate the integral between the leaf faces
		if (aHasChildren || bHasChildren)
		{
			FaceProcZ(aHasChildren ? nodeA.childIdx + 4 : nodeIdxA_, bHasChildren ? nodeB.childIdx + 0 : nodeIdxB_, jobQueue_);
			FaceProcZ(aHasChildren ? nodeA.childIdx + 5 : nodeIdxA_, bHasChildren ? nodeB.childIdx + 1 : nodeIdxB_, jobQueue_);
			FaceProcZ(aHasChildren ? nodeA.childIdx + 6 : nodeIdxA_, bHasChildren ? nodeB.childIdx + 2 : nodeIdxB_, jobQueue_);
			FaceProcZ(aHasChildren ? nodeA.childIdx + 7 : nodeIdxA_, bHasChildren ? nodeB.childIdx + 3 : nodeIdxB_, jobQueue_);
		}
		else
		{
			// Order nodes based on dim
			const usize nodeIdxA = nodes[nodeIdxA_].aabb.min()(2) < nodes[nodeIdxB_].aabb.min()(2) ? nodeIdxA_ : nodeIdxB_;
			const usize nodeIdxB = nodes[nodeIdxA_].aabb.min()(2) < nodes[nodeIdxB_].aabb.min()(2) ? nodeIdxB_ : nodeIdxA_;

			// Determine if we've done this integral before
			if (procMap.find(std::pair<usize, usize>(nodeIdxA, nodeIdxB)) != procMap.end())
			{
				return;
			}
			else
			{
				procMap.emplace(std::pair<usize, usize>(nodeIdxA, nodeIdxB), true);
			}

			// New job is unique
			ContinuityThreadPool::Input newJob;
			newJob.nodeIdxs = std::pair<usize, usize>(nodeIdxA, nodeIdxB);
			newJob.dim      = 2;
			jobQueue_.push(newJob);
		}
	}


	void Octree::TickContinuityThread(const ContinuityThreadPool::InitialData* const threadData_, const usize threadIdx_)
	{
		while (1)
		{
			// Check if we need to exit
			if (threadData_->shutdownAtom->load())
			{
				break;
			}

			// Find new job
			ContinuityThreadPool::Input newJob;
			bool haveJob = false;
			threadData_->inputQueueMutex->lock();
			{
				haveJob = !threadData_->inputQueue->empty();
				if (haveJob)
				{
					newJob = threadData_->inputQueue->front();
					threadData_->inputQueue->pop();
				}
			}
			threadData_->inputQueueMutex->unlock();

			// Nothing to do
			if (!haveJob)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(25));
				continue;
			}

			// Consts
			const usize nodeIdxA = newJob.nodeIdxs.first;
			const usize nodeIdxB = newJob.nodeIdxs.second;

			// Perform integration
			if (nodes[nodeIdxA].depth == nodes[nodeIdxB].depth)
			{
				EvaluateSharedFaceIntegralAnalytically(nodeIdxA, nodeIdxB, newJob.dim, (threadData_->matTriplets)[threadIdx_]);
			}
			else
			{
				EvaluateSharedFaceIntegralNumerically(nodeIdxA, nodeIdxB, newJob.dim, (threadData_->matTriplets)[threadIdx_]);
			}
		}
	}


	void Octree::RunContinuityThreadPool(std::vector<Eigen::Triplet<f64>>& integralMatrixTriplets_)
	{
		// Create structs
		ContinuityThreadPool::InitialData threadData;
		threadData.config          = config;
		threadData.inputQueue      = new (std::queue<ContinuityThreadPool::Input>);
		threadData.inputQueueMutex = new (std::mutex);
		threadData.shutdownAtom    = new (std::atomic<bool>);
		threadData.matTriplets     = new std::vector<Eigen::Triplet<f64>>[config.threadCount];
        threadData.shutdownAtom->store(false);

		// Go through octree and allocate jobs to threads. Do this before threads are up to avoid pointless locking/unlocking
		for (usize i = 0; i < nodes.size(); ++i)
		{
			NodeProc(i, *threadData.inputQueue);
		}

		// Start threads
		ContinuityThreadPool threadPool;
		threadPool.StartThreads(&threadData, this);

		// Wait for threads to finish
		bool finished = false;
		while (!finished)
		{
			std::this_thread::sleep_for(std::chrono::milliseconds(25));

			threadData.inputQueueMutex->lock();
			{
				finished = threadData.inputQueue->empty();
			}
			threadData.inputQueueMutex->unlock();
		}

		// Stop threads
		threadData.shutdownAtom->store(true);
		threadPool.StopThreads();

		// Copy results and force a free whilst doing so due to memory limits
		for (usize i = 0; i < config.threadCount; ++i)
		{
			integralMatrixTriplets_.reserve(integralMatrixTriplets_.size() + threadData.matTriplets[i].size());
			integralMatrixTriplets_.insert(integralMatrixTriplets_.end(), threadData.matTriplets[i].begin(), threadData.matTriplets[i].end());
			std::vector<Eigen::Triplet<f64>>().swap(threadData.matTriplets[i]);
		}

		// Cleanup
		delete   (threadData.inputQueue);
		delete   (threadData.inputQueueMutex);
		delete   (threadData.shutdownAtom);
		delete[] (threadData.matTriplets);
	}


	void Octree::PerformContinuityPostProcess(const usize nCoeffs_)
	{
		// Determine sparse mat entries
		std::vector<Eigen::Triplet<f64>> integralMatrixTriplets;
		RunContinuityThreadPool(integralMatrixTriplets);

		// Regularise
		for (usize i = 0; i < nCoeffs_; ++i)
		{
			integralMatrixTriplets.push_back(Eigen::Triplet<f64>((Eigen::SparseMatrix<f64>::StorageIndex)i, 
                                                                 (Eigen::SparseMatrix<f64>::StorageIndex)i, 
                                                                  config.continuity.strength));
		}

		// Create sparse mat from entries
		Eigen::SparseMatrix<f64, Eigen::RowMajor> integralMatrix;
		integralMatrix.resize(nCoeffs_, nCoeffs_);
		integralMatrix.setFromTriplets(integralMatrixTriplets.begin(), integralMatrixTriplets.end());
		std::vector<Eigen::Triplet<f64>>().swap(integralMatrixTriplets);

		// Create vector for old coeffs
		Eigen::VectorXd oldCoeffs;
		oldCoeffs.resize(nCoeffs_);
		memcpy(oldCoeffs.data(), coeffStore, sizeof(f64) * nCoeffs_);
		oldCoeffs *= config.continuity.strength;

		// If OpenMP enabled, leverage additional threads
		Eigen::initParallel();
		Eigen::setNbThreads((i32)config.threadCount);

		// Solve system and copy new coeffs over
		Eigen::ConjugateGradient<Eigen::SparseMatrix<f64, Eigen::RowMajor>, Eigen::Lower | Eigen::Upper,
        Eigen::IncompleteCholesky<f64, Eigen::Lower | Eigen::Upper, Eigen::NaturalOrdering<int>>> extSolver(integralMatrix);

		extSolver.setTolerance(EPSILON_F32);
		const Eigen::VectorXd newCoeffs = extSolver.solveWithGuess(oldCoeffs, oldCoeffs);
		memcpy(coeffStore, newCoeffs.data(), sizeof(f64) * nCoeffs_);

		// No side effects!
		Eigen::setNbThreads(1);
	}
}