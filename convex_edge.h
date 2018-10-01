/*********************
class defined for convex arc/edge on 
atoms (composing the saddle face
and convex faces)
***********************/
#ifndef _CONVEX_EDGE_
#define _CONVEX_EDGE_

#include "Vector3d.h"

class convex_edge{
	public:
	CVector3d pos1, pos2; //starting and ending vertices for convex edge
	int atom_index; //index of the adjacent atom

	convex_edge(CVector3d p1, CVector3d p2, int index);
};

#endif