#ifndef _GRID_POINT
#define _GRID_POINT

#include <vector>
#include "../utility/Vector3d.h"
#include "../utility/Types.h"

class grid_point{
public:
	int grid_x,grid_y,grid_z;  //index in the grid
	CVector3d grid_pos;
	bool is_processed;
	bool is_irregular; //itself is outside and its neighboring grid point is inside or it is inside and its neighboring grid point is outside
	bool is_inside; // mark whether it is inside the surface
	bool is_outside_aug_atom; //mark whether it is outside the augmented atoms
//	bool is_inside_atom; //mark whether it is inside the atoms
//	bool is_inside_concave_sphere;
//	bool is_inside_tori; //mark whether it is inside the tori
	bool is_inside_VS_or_dualTet; //mark whether it is in the Visibility sphere of torus or the dual tet for concave face
	bool is_regular;

	std::vector<AtomIter> possible_atom_interior; //stores the atoms which possibly includes current grid point
	std::vector<TorusIter> possible_torus_interior;
	std::vector<Concave_SphereIter> possible_concavesphere_interior;

	grid_point();
	grid_point(int X, int Y, int Z);
	~grid_point();

	void initialize_grid_point(double x_step,double y_step,double z_step, double x_min, double y_min, double z_min);

};
#endif