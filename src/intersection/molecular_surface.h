/******************************
April/16/2014 

This header file is for computing the molecular surface based on paper

"Analytical Molecular Surface Calculation" 
   By Michael L. Connolly 1983

Format for Input:
// ".xyzr" file containing the atom information
//probe raidus
//grid size
//extension of the bounding box


Algorithm of the paper
  
(1) Torus construction:
construct torus for each pair of neighboring atoms i,j (when i>j, in case of construct the same torus twice)
(equations in Table 2, note the case when the square root of 'torus radius r_ij' is imaginary 
<when d_ij>r_i+r_j+2r_p or d_ij < fabs(r_i - r_j) >
and then no torus will be constructed)

(2)Probe placement/ Mark the torus
for each torus constructed for atom i,j,

a. when the denominator of 'Base point b_ijk' is zero (i.e, the centers of atom i,j,k are
co-linear), there is no probe placement

b.when the case a is checked, we also need to check the square root of 'Probe height h_ijk'. If it is imaginary:
<1> for EVERY neighboring atom l of atom i,j (or there is no neighboring atom of atom i,j), atom l is too far away
from atom i,j to collide, then the torus of atom i,j is marked as free torus
<2> for ANY neighboring atom l of atom i,j, if atom l is located between atom i, j, this torus is (partially) buried in the atom k

(3) Construct concave sphere and edges (for i < j < k)
a. for a torus which is not blocked and non-free,
when the quantity under 'Probe height h_ijk' is positive, we need to check the collision condition for each of the two possible
probe position. (collision detection with all the neighboring atoms of atom i,j,k)

b. when it survives the collision detection, we will have three concave edges (appended to the torus between that pair of atoms) and 
one concave face.

(4) Construct saddle face
a. Free torus is a special case where there is no concave edges and its convex edges are complete contact circles

b.For other torus:
<1> reorder the concave edges such that they are in clockwiswe order when viewed along the torus axis from atom i to j;
<2> construct the convex edges/arcs (which on atoms) such that when looking along the torus axis from atom i to j, the circle
on atom i is COUNTER-CLOCKWISE and on atom j is CLOCKWISE
<3> construct the saddle faces by choosing a pair of adjacent concave edges and corresponding convex edges so that the orientation
of the resulting rectangle is CLOCKWISE when viewed from outside of the molecule. When each saddle face constructed, a pair of
convex arc is generated and appended to the corresponding atom.

(5) Construct convex faces
a.special case: when the atom's entire surface is accessible

b.for others:
<1> group the convex edges into cycles with orientation
<2> check for each pair of cycles and store their interior and exterior info
<3> construct the convex faces given the info <2>

*********************************/
#ifndef _MOLECULAR_SURFACE_
#define _MOLECULAR_SURFACE_

#include <iostream>
#include <iomanip> 
#include <fstream> 
#include <algorithm>
#include <string>
#include <sstream>
#include <list>
#include <vector>
#include <cmath>
#include <time.h>
#include <omp.h>

#include "../utility/Types.h"
#include "../utility/Vector3d.h"
#include "../utility/Jenkins_Traub.h"

#include "../surface/atom_sphere.h"
#include "../surface/torus.h"
#include "../surface/concave_edge.h"
#include "../surface/concave_sphere.h"
#include "../surface/saddle_face.h"

#include "../grid/grid_point.h"
#include "../grid/grid_edge.h"


#undef max
#undef min
#define  TEST 0
#define  NEIGHBOR_LENGTH 50
//#define BLOCK_SIZE 127

class intersection_info{
public:
    std::vector<CVector3d> m_grid_inside; // this is int, but use CVector3d for convenience
    std::vector<CVector3d> m_grid_outside; // this is int, but use CVector3d for convenience
    std::vector<CVector3d> m_intersection;
    std::vector<CVector3d> m_normal;
    std::vector<std::vector<int>> m_adj_atoms; 
};

class molecular_surface{
public:
	double m_probe_radius;
	int m_N_atoms;  //# of atoms
	std::vector<Atom *> m_atoms_vec;
	std::vector<Torus *> m_torus_vec;
	std::vector<Concave_Sphere *> m_concave_spheres;
	std::vector<concave_edge *> m_concave_edges;
	std::vector<saddle_face *> m_saddle_faces;

	std::vector<grid_point *> m_grid_points; 
	std::vector<grid_edge *> m_grid_edges;

	//vector<bool> m_grid_status;
	//vector<bool> m_grid_fix_volume;
	//vector<bool> m_grid_regularity;

	double m_min_x,m_min_y,m_min_z;  //bounding box
	double m_max_x,m_max_y,m_max_z;
	//double t_min_x, t_min_y, t_min_z;  //total bounding box
	//double t_max_x, t_max_y, t_max_z;
	double m_grid_resolution;
	int m_x_num,m_y_num,m_z_num; //number of grid points in each dimension
	//int t_x_num, t_y_num, t_z_num; //total number of grid points in each dimension
	double m_x_step,m_y_step,m_z_step;
	//int block_dim_x, block_dim_y, block_dim_z;
	//int a, b, c; //block index

	std::vector<CVector3d> m_intersection_points; //store the output: location of intersection points
	std::vector<CVector3d> m_intersection_normals; //store the normal direction of the intersection points
	std::vector<grid_pointIter> m_intersection_interior;
	std::vector<grid_pointIter> m_intersection_outside;
	std::vector<cell_unit> m_intersection_cell;
	std::vector<int> m_intersection_types; // 0: atom; 1:torus; 2:concave

	//CVector3d m_mesh_min,m_mesh_max; //store the bounding box of mesh
	//double m_extend_bounding;
	

	double m_surface_area,m_surface_volume;
	double m_area_unit, m_volume_unit; //area unit: h^2; volume unit: h^3
	std::vector<double> m_partition_area; //store the area on each atom;
	std::map<int, int> m_local_to_global_atom_idx;
	double m_area_from_partition;
	double m_another_area_from_partiion;


	molecular_surface();
	~molecular_surface();


CVector3d grid_pos(int X,int Y,int Z);
//return the global index of the grid point
int grid_index(int X,int Y, int Z);

//construct the grid points and grid edge
void initialize_grid();


/**********************constructing the molecular surface starts**********************************/
/**********************constructing the molecular surface starts**********************************/
/**********************constructing the molecular surface starts**********************************/

//initialize the molecular surface
void initialize_molecular_surface();



/*********************
construct the saddle face by choosing
the right concave edges and convex edges

(meaning the current concave edge should 
pointing from vi to vj, while the next concave edge
should pointing from vj to vi)
*********************/
void construct_saddle_each_torus(TorusIter current_torus);

void construct_saddle_faces();

//fix the duplex concave spheres (when the probe is touching more than three atoms at the same time).
void check_duplex_concave();

void store_duplex_concaves(Concave_SphereIter current_concave,std::vector<Concave_SphereIter> &full_duplex_concaves);

//construct the torus
void construct_torus();


//check whether atom_i and atom_j are neighboring atom (i.e, their distance  fabs(r_j-r_i)< dij < r_i + r_j + 2*r_p)
bool is_neighbor_atom(AtomIter atom_i, AtomIter atom_j);

/*********************************************
	mark if it is a free torus or not.  A torus (t_ij) is free
	when for EVERY mutual neighbor atom k of atom i,j,
	there is no intersection with atom k

********************************************/
void mark_torus();


/***************************************************
collision detection when constructing concave sphere
probe_1, probe_2 : possible location of probe sphere center
i,j,k : index for atom i,j,k
is_collide_1, is_collide_2 : return the collision detection result

Check whether the location of probe sphere will collide with the mutual neighbors of atom i, j,k

***************************************************/
void collision_detection(CVector3d probe_1,CVector3d probe_2,AtomIter atom_i,AtomIter atom_j, AtomIter atom_k,bool &is_collide_1, bool &is_collide_2);

void construct_single_concave(CVector3d probe_center,AtomIter atom_i,AtomIter atom_j,AtomIter atom_k,TorusIter torus_ij,TorusIter torus_jk, TorusIter torus_ki);

void probe_placement();

/**********************constructing the molecular surface ends**********************************/
/**********************constructing the molecular surface ends**********************************/
/**********************constructing the molecular surface ends**********************************/


//check whether current gird point is the atom sphere or not
bool is_inside_atom(AtomIter current_atom,CVector3d grid_loc);

bool is_inside_concave_sphere(CVector3d grid_pos,CVector3d sphere_center);

/********************************************************************************
check whether test_point is inside the tetrahedron composed by p1, p2, p3, p4

Here we use the method from :
http://steve.hollasch.net/cgindex/geometry/ptintet.html
**********************************************************************************/
bool is_inside_tetrahedron(CVector3d p1,CVector3d p2, CVector3d p3, CVector3d p4, CVector3d test_point);

bool is_f_InsideTorus_positive(double x,double y,double z,double R,double r);

bool is_f_quartic_equation_positive(double x,double y,double z,double R,double r);

bool is_located_in_selfintersecting_part(double x,double y,double z,double R,double r);

bool is_located_within_lemon_part(double x,double y,double z,double R,double r);

bool is_located_on_lemon_part(double x,double y,double z,double R,double r);

bool is_covered_by_saddles(CVector3d grid_point,TorusIter current_torus);

bool is_inside_Visibility_Sphere_old(TorusIter current_torus,CVector3d point_pos);

bool is_inside_Visibility_Sphere(TorusIter current_torus,CVector3d point_pos);

bool validate_atom_intersection(AtomIter atom,CVector3d pos);

/**********************labeling the grid points starts**********************************/
/**********************labeling the grid points starts**********************************/
/**********************labeling the grid points starts**********************************/
void process_each_atom_newer(AtomIter current_atom);

void process_each_concave_sphere_newer(Concave_SphereIter current_concave);

void process_each_torus_newer(TorusIter current_torus);

//if return true, then this grid point is outside the molecular surface
bool is_inside_accessible_torus_old(TorusIter current_torus,grid_pointIter current_point);

//if return true, then this grid point is outside the molecular surface
bool is_inside_accessible_torus(TorusIter current_torus,grid_pointIter current_point);

bool is_inside_accessible_torus_with_output(TorusIter current_torus,grid_pointIter current_point);

void check_regular();

void check_the_grid_points();

//detect the other unprocessed grid points by counting the number of intersection points
void process_others_ray_counting();

void process_ray_racing_each_point(grid_pointIter current_point);

void process_each_augmented_atom(AtomIter current_atom);


int compute_surface_intersection_num(grid_pointIter current_gridpoint,grid_pointIter processed_ngb);

int compute_intersection_single_atom(grid_pointIter point_inside, grid_pointIter point_outside,AtomIter atom,std::vector<double> &t);

int compute_atom_intersection_num(grid_pointIter point_inside, grid_pointIter point_outside);

int compute_torus_intersection_num(grid_pointIter point_inside,grid_pointIter point_outside);

//----------------------------------------------------------------------------
bool solveQuartic(double a, double b, double c, double d, double e, double *root);

/**********************************************
compute line and torus intersection.

Coefficients are assembled based on
"http://www.emeyex.com/site/projects/raytorus.pdf"

note here pos_interior and pos_outside are in the PARAMETRIZED domain 


return: 
t as the intersection point of p + t *td; 
inter_normal: the normal of the intersection point on the PARAMETRIZED domain
*********************************************/
int compute_line_torus_intersection(TorusIter current_torus, CVector3d grid_pos_interior,CVector3d grid_pos_outside,CVector3d pos_interior, 
									CVector3d pos_outside,std::vector<double> &t,std::vector<CVector3d> &inter_normal);

int compute_concave_intersection_num(grid_pointIter point_inside,grid_pointIter point_outside);

/**********************labeling the grid points ends**********************************/
/**********************labeling the grid points ends**********************************/
/**********************labeling the grid points ends**********************************/



/********************compute intersection starts**************************************/
/********************compute intersection starts**************************************/
/********************compute intersection starts**************************************/
int process_grid_intersection_new(grid_edgeIter current_edge);

void process_grid_atom_intersection_new(grid_pointIter point_inside, grid_pointIter point_outside,std::vector<CVector3d> &intersect_points,
										std::vector<CVector3d> &intersect_normal,std::vector<int> &intersect_type,std::vector<cell_unit> &intersect_cell);

void process_grid_torus_intersection_new(grid_pointIter point_inside,grid_pointIter point_outside,TorusIter current_torus, std::vector<CVector3d> &intersect_points,
										 std::vector<CVector3d> &intersect_normal,std::vector<int> &intersect_type,std::vector<cell_unit> &intersect_cell);

void process_grid_edge_concave_new(grid_pointIter point_inside,grid_pointIter point_outside,std::vector<CVector3d> &intersect_points,
								   std::vector<CVector3d> &intersect_normalvec,std::vector<int> &intersect_type,std::vector<cell_unit> &intersect_cell);

/********************compute intersection ends**************************************/
/********************compute intersection ends**************************************/
/********************compute intersection ends**************************************/

void clean_memory();

/**********************
	here we check whether it is inside the torus or outside 
	 the parametrized expression for torus is:
	 (R - sqrt(x*x + y*y) ) ^2 + z^2 = r^2

	 To be in the 'interior' of this torus patch, it should require
	 (1) [R - sqrt(x*x+y*y)]^2 + z^2 - r^2 >0
	 (2) within the visibility sphere (Fig 7. of paper "Interactive Visualization of Molecular Surface Dynamics")
************************************************************/
bool is_inside_torus(TorusIter current_torus,grid_pointIter current_point);

/**********************
	here we check whether it is inside the torus or outside 
	 the parametrized expression for torus is:
	 (R - sqrt(x*x + y*y) ) ^2 + z^2 = r^2

	 To be in the 'interior' of this torus patch, it should require
	 (1) [R - sqrt(x*x+y*y)]^2 + z^2 - r^2 >0
	 (2) within the visibility sphere (Fig 7. of paper "Interactive Visualization of Molecular Surface Dynamics")
************************************************************/
bool is_inside_torus_new(TorusIter current_torus,grid_pointIter current_point);

bool is_outside_concave_sphere(TorusIter current_torus,CVector3d current_point);

//find the closest neighboring grid point which is already processed
bool search_the_processed_neighbor(grid_pointIter current_point, grid_pointIter &processed_neighbor,int &found_k,int &mode);

bool solveQuadraticOther(double a, double b, double c, double &root);

//----------------------------------------------------------------------------
bool solveQuadratic(double a, double b, double c, double &root);

bool solveCubic(double a, double b, double c, double d, double &root);

double f_value(double a ,double b, double c, double d, double e, double x);

double f_derivative_value(double a ,double b, double c, double d, double e, double x);

bool solve_newton_method(double a ,double b, double c, double d, double e,double x_0, double &t);

bool is_inside_spherical_triangle(Concave_SphereIter current_concave,CVector3d test_point);

void accum_surface_area(grid_pointIter point1,grid_pointIter point2,CVector3d intersect_normal);

int partition_area(int a, int b, int c, int block_size);

int process_each_grid_edge(grid_pointIter point_inside, grid_pointIter point_outside,std::vector<CVector3d> &intersect_points,
			   std::vector<CVector3d> &intersect_normal,std::vector<int> &intersect_type,std::vector<cell_unit> &intersect_cell);

void write_intersections(intersection_info& intersection, int block_x, int block_y, int block_z, int block_size);
 
//given the status of the each grid point (inside and outside), check the missing intersections
void test_intersection();
};

#endif
