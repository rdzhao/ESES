/*********************
Class for torus surface when probe sphere
touches with two atom spheres

*******************/
#ifndef _TORUS_
#define _TORUS_

#include <vector>
#include "../utility/Vector3d.h"
#include "../utility/Types.h"
#include "../utility/Matrix44.h"

class Torus{
public:
	Torus();
	Torus(AtomIter atomi,AtomIter atomj,double probe_radius);
	~Torus();

	AtomIter atom_i, atom_j; //pointer of corresponding atom
	double torus_r;  //radius of the torus
	CVector3d torus_center; //center of the torus  
	CVector3d torus_axis;  //axis of the torus, direction pointing from the center of atom i to center of atom j
	bool is_free;  //mark the condition of torus (free and non-free)
	bool is_blocked; //mark whether this torus in entirely inside one atom
	
	CVector3d circle_center_i, circle_center_j; //center of the contact circle on atom i and j
	double phi_i, phi_j;
	CVector3d rotation_axis;  //for visualization 
	double rotation_angle;     //for visualization 

	std::vector<AtomIter> mutual_atom_list; //pointer of the mutual neighbors of atom i, j
	std::vector<concave_edgeIter>concave_edge_list;  //store the concave edge list for torus
	std::vector<saddle_faceIter> saddle_face_list; //store the corresponding saddle faces	
	
	CMatrix44 para_to_torus, torus_to_para;
	CVector3d vs_sphere_center; double vs_sphere_r; //Fig 7. of paper "Interactive Visualization of Molecular Surface Dynamics"


	struct order_elem{
		double angle;
		int index;
	};

	
	
	void remove_duplicate_concave_edge();
	void reorder_concave_edge(); //reorder the concave edges such that it is clockwise when viewing from atom i to atom j
	double angle_between_vectors(CVector3d a, CVector3d b,CVector3d n);
	void selectionSort(struct order_elem data[], int lenD);
	CMatrix44 return_rotation_matrix(double c,double s, CVector3d rotation_axis);
	void generate_true_self_intersecting_boundary(double probe_radius);


};

#endif