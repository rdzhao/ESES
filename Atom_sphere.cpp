#include "Atom_sphere.h" 
#include "Torus.h"

Atom::Atom(){
	center = CVector3d(0.0,0.0,0.0);
	r = 0.0;

}

Atom::Atom(CVector3d center_pos,double radius,int n){
	center = center_pos;
	r = radius;
	index = n;
	is_inside = false;
	
}

Atom::~Atom(){
	
}

void Atom::initialize_atom(int num_atoms){
	//adj_torus_list.resize(num_atoms);
	is_chosen = false;
}