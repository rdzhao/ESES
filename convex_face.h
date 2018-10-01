#ifndef _CONVEX_FACE_
#define _CONVEX_FACE_

#include "types.h"

class convex_face{
public:
	//convex face is defined by its boundary circles
	std::list<atom_cycleIter> boundary_cycles;

	convex_face(){ }


};
#endif