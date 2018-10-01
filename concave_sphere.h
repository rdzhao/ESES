#ifndef _CONCAVE_SPHERE_
#define _CONCAVE_SPHERE_

#include "Vector3d.h"
#include "types.h"


class Concave_Sphere{
public:
	int index_i,index_j,index_k;  //index of corresponding atoms
	CVector3d center_i, center_j, center_k;  //center of the corresponding atoms
	double r_i,r_j,r_k;
	concave_edgeIter edge_ij,edge_jk,edge_ki;
	CVector3d probe_center;


	//for probe touching more than three atoms at the same time
	std::vector<Concave_SphereIter> duplex_concave_faces;
	bool has_check_duplex;

	Concave_Sphere();
	Concave_Sphere(concave_edgeIter e1, concave_edgeIter e2,concave_edgeIter e3, int i,int j, int k,CVector3d a_i,CVector3d a_j, CVector3d a_k,double radius_i,
		double radius_j,double radius_k,CVector3d probe);

};
#endif