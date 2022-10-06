#include "HP/Node.h"

namespace SDF
{
    Node::Node()
    {
	    childIdx = -1;

	    aabb.setEmpty();

	    basis.coeffs = nullptr;
	    basis.degree = BASIS_MAX_DEGREE + 1;

	    depth = TREE_MAX_DEPTH + 1;
    }
}
