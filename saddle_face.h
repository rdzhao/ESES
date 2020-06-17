#ifndef _SADDLE_FACE_
#define _SADDLE_FACE_

#include "types.h"
#include "Vector3d.h"
#include "Matrix44.h"


class saddle_face{
public:
	
	concave_edgeIter concave_edge_1, concave_edge_2;
	CVector3d boundary_i_1,boundary_i_2; //boundary unit vector of the saddle face on the touching cycle of atom i
	CVector3d boundary_cross; //boundary_i_1.cross(boundary_i_2)
	bool boundary_check;//if true, then the true saddle face is located between boundary_i_1 and boundary_i_2, otherwise it is located in the
	//complimentary part

	
	double theta_s;
	double angle_z,angle_x; //angle need to rotate to align z and x direction
	CMatrix44 para_to_saddle;
	
	saddle_face();
	saddle_face(concave_edgeIter e_1,concave_edgeIter e_2,TorusIter current_torus,double probe_radius);
	~saddle_face();

	/****************************
	the orientation of edges in this saddle face:
	convex edge: e_i  (vi_2---> vi_1)
	concave edge: e_1 (vi_1 --> vj_1)
	convex edge: e_j (vj_1 --> vj_2)
	concave edge: e_2 (vj_2 --> vi_2)

	********************************/
	
	double return_rotation_angle(TorusIter current_torus,CVector3d previous_direction);
	CMatrix44 return_rotation_matrix(double rotation_angle, CVector3d rotation_axis);
	void compute_transformation(concave_edgeIter e_1,concave_edgeIter e_2,TorusIter current_torus);

	//void return_rotation_matrix(double rotation_angle, CVector3d rotation_axis, double **Q);

};
#endif