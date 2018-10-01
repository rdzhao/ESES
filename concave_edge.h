/************************
class defined for concave edge
defined on probe sphere (which 
composes the concave face
and saddle face)

************************/
#ifndef _CONCAVE_EDGE_
#define _CONCAVE_EDGE_

#include "types.h"
#include "Vector3d.h"
class concave_edge{
public:

	CVector3d pos1, pos2; //starting and ending point 
	int atom1,atom2; //index of corresponding atom
	CVector3d n_ijk; //concave arc plane normal vector
	Concave_SphereIter current_concave;
	
	concave_edge();
	concave_edge(CVector3d p1,CVector3d p2,int n1,int n2);
	void initial_concave_edge(CVector3d probe_center,CVector3d torus_center,CVector3d torus_axis);

};



#endif