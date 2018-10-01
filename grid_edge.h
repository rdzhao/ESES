#ifndef _GRID_EDGE
#define _GRID_EDGE

#include <vector>
#include "Vector3d.h"
#include "types.h"

class grid_edge{
public:
	grid_pointIter point1,point2;  //pointers of the two ends
	grid_edge(grid_pointIter p1, grid_pointIter p2);
	~grid_edge();

};
#endif