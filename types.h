#ifndef _TYPES_H_
#define _TYPES_H_

#include <list>
#include <vector>

#define ZERO_TOL 1e-12
#define PI 3.14159265358979323846264

class Torus;
class Atom;
class concave_edge;
class saddle_face;
class Concave_Sphere;
class grid_point;
class grid_edge;


//typedef std::list<Atom>::iterator AtomIter;
typedef Torus* TorusIter;
typedef Atom* AtomIter; 
typedef Concave_Sphere* Concave_SphereIter;
typedef saddle_face* saddle_faceIter;
typedef concave_edge* concave_edgeIter;
typedef grid_point* grid_pointIter;
typedef grid_edge* grid_edgeIter;

typedef std::vector<Atom *>::iterator AtomPointerIter;


struct cell_unit 
{
	AtomIter current_atom;
	TorusIter current_torus;
	Concave_SphereIter current_concave;
};

#endif