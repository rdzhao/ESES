#ifndef _TYPES_H_
#define _TYPES_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <list>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <cmath>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <limits>
#include <omp.h>

#define ZERO_TOL 1e-12
#define PI 3.14159265358979323846264

const int BLOCK_SIZE = 127;

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

struct cell_unit{
	AtomIter current_atom;
	TorusIter current_torus;
	Concave_SphereIter current_concave;
};

enum AtomType{ATOM_C=0, ATOM_H, ATOM_O, ATOM_N, ATOM_P, ATOM_S, ATOM_K};
enum PatchType{PATCH_CONVEX=0, PATCH_SADDLE, PATCH_CONCAVE};

const std::unordered_map<char, AtomType> atmap = {
    {'C', ATOM_C},
    {'H', ATOM_H},
    {'O', ATOM_O},
    {'N', ATOM_N},
    {'P', ATOM_P},
    {'S', ATOM_S}  
};

const std::string CONSOLE_RED("\033[0;31m");
const std::string CONSOLE_GREEN("\033[1;32m");
const std::string CONSOLE_YELLOW("\033[1;33m");
const std::string CONSOLE_CYAN("\033[0;36m");
const std::string CONSOLE_MAGENTA("\033[0;35m");
const std::string CONSOLE_WHITE("\033[0m");

class SVector2I{
public:
    int x;
    int y;

    SVector2I(int xx, int yy) : x(xx), y(yy){}

    const bool operator==(const SVector2I& v){
        return (x==v.x && y==v.y);
    }
};

class SVector3I{
public:
    int x;
    int y;
    int z;

    SVector3I(int xx, int yy, int zz) : x(xx), y(yy), z(zz){}

    const bool operator==(const SVector3I& v){
        return (x==v.x && y==v.y && z==v.z);
    }
};

class SVector2ICompare{
public:
    const bool operator()(const SVector2I& v1, const SVector2I& v2) const {
        if(v1.x != v2.x) return v1.x < v2.x;
        else return v1.y < v2.y;
    }
};

class SVector3ICompare{
public:
    const bool operator()(const SVector3I& v1, const SVector3I& v2) const {
        if(v1.x != v2.x) return v1.x < v2.x;
        else if(v1.y != v2.y) return v1.y < v2.y;
        else return v1.z < v2.z;
    }
};

#endif