#include "convex_edge.h"

convex_edge::convex_edge(CVector3d p1, CVector3d p2, int index){
	pos1 = p1;
	pos2 = p2;
	atom_index = index;
}