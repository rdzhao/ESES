#include "grid_point.h"

grid_point::grid_point(){

}

grid_point::grid_point(int X, int Y, int Z){
	grid_x = X;
	grid_y = Y;
	grid_z = Z;

	is_inside = false; //default	
	is_outside_aug_atom = true;
//	is_inside_atom = false;
//	is_inside_concave_sphere = false;
	is_inside_VS_or_dualTet = false;
//	is_inside_tori = false;
	is_processed = false;

	is_irregular = false;
	is_regular = false;
}

grid_point::~grid_point(){
	possible_atom_interior.clear();
	possible_torus_interior.clear();
	possible_concavesphere_interior.clear();
}
void grid_point::initialize_grid_point(double x_step,double y_step,double z_step, double min_x, double min_y, double min_z){
	grid_pos = CVector3d(grid_x * x_step + min_x,grid_y * y_step + min_y, grid_z* z_step + min_z);
}