#include "Node.h"

namespace SDF
{
	Node::Node()
	{
		childIdx = -1;

		aabb.setEmpty();

		basis.coeffs = nullptr;
		basis.degree = 0;

		depth = -1;
	}


	Node::~Node()
	{

	}
}