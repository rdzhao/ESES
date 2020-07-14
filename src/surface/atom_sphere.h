/*********************
class of atom sphere
********************/
#ifndef _ATOM_
#define _ATOM_

#include <cmath>
#include <list>
#include <vector>
#include <map>
#include "../utility/Vector3d.h"
#include "../utility/Types.h"

class Atom{
public:

	CVector3d center;
	double r;
	int index; //index of atom
	std::vector<AtomIter> neighbor_list; //pointer of the neighboring atom
	bool is_inside; //true if this atom is completely inside another atom
	//std::vector<TorusIter> adj_torus_list; //the array size is number of atoms, store the pointer to the corresponding torus
	std::map<int, TorusIter> adj_torus_list;
	std::vector<CVector3d> boundary_cicle_centers;
	std::vector<CVector3d> boundary_circle_normals;

	bool is_chosen; //for rendering 

	Atom();
	Atom(CVector3d center_pos,double radius,int n);
	~Atom();
	void initialize_atom(int num_atoms);
	
	


};

#endif