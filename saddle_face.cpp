#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream> 
#include "saddle_face.h"
#include "concave_edge.h"
#include "Torus.h"

saddle_face::saddle_face(){

}

saddle_face::~saddle_face(){

}

saddle_face::saddle_face(concave_edgeIter e_1,concave_edgeIter e_2,TorusIter current_torus,double probe_radius){

	concave_edge_1 = e_1;// the concave edge pointing from the point at atom_i to atom_j
	concave_edge_2 = e_2; // the concave edge pointing from the point at atom_j to atom_i

	//compute the transformation to the standard 'torus'
	compute_transformation(e_1,e_2,current_torus);

	boundary_i_1 = e_1->pos1 - current_torus->circle_center_i;
	boundary_i_2 = e_2->pos2 - current_torus->circle_center_i;
	boundary_i_1.Normalize();
	boundary_i_2.Normalize();
	boundary_cross = (boundary_i_1.Cross(boundary_i_2));

	if(boundary_cross.Dot(current_torus->torus_axis) < 0)
		boundary_check = false;
	else
		boundary_check = true;
	
}


void saddle_face::compute_transformation(concave_edgeIter e_1,concave_edgeIter e_2,TorusIter current_torus){
	CVector3d n_ijk = e_1->n_ijk;
	CVector3d n_ijl = e_2->n_ijk;
	double check_val = (n_ijk.Cross(n_ijl)).Dot(current_torus->torus_axis);

	//the vector representing the vector pointing from the atom i center to the point on the atom i 
	CVector3d new_X = e_1->pos1 - current_torus->torus_center;
	new_X = new_X - (new_X.Dot(current_torus->torus_axis))*current_torus->torus_axis;
	new_X.Normalize();

	angle_x = return_rotation_angle(current_torus,new_X);

	CVector3d new_X2 = e_2->pos1 - current_torus->torus_center;
	new_X2 = new_X2 - (new_X2.Dot(current_torus->torus_axis))*current_torus->torus_axis;
	new_X2.Normalize();




	double angle_x2 = return_rotation_angle(current_torus,new_X2);



	angle_z = -current_torus->rotation_angle;
	theta_s = angle_x2 - angle_x;

	//if((angle_x >0 && angle_x2<0) || (angle_x <0 && angle_x2>0)){
	if(fabs(theta_s)>PI)
		check_val *= -1;

	if(check_val <0){
		if(theta_s > 0)
			theta_s -= 2*PI;
		else
			theta_s += 2*PI;
	}


	//special case i.e concave_edge_1 and concave_edge_2 are the opposite edge
	if((n_ijk.Cross(n_ijl)).Length() < ZERO_TOL)
		theta_s = 0;

	//compute the transformation matrix
	CMatrix44 rotation_matrix1 = return_rotation_matrix(angle_z,current_torus->rotation_axis);
	CMatrix44 rotation_matrix2 = return_rotation_matrix(angle_x,CVector3d(0.0,0.0,1.0));
	para_to_saddle = rotation_matrix1.MultRight(rotation_matrix2);


}

CMatrix44 saddle_face::return_rotation_matrix(double rotation_angle, CVector3d rotation_axis){
	double x = rotation_axis[0];
	double y = rotation_axis[1];
	double z = rotation_axis[2];
	double c = cos(rotation_angle);
	double s = sin(rotation_angle);

	double Q[4][4];

	Q[0][0] = x*x*(1-c)+c;    Q[0][1] = x*y*(1-c)-z*s; Q[0][2] = x*z*(1-c)+y*s; Q[0][3] = 0.0;
	Q[1][0] = y*x*(1-c)+z*s; Q[1][1] = y*y*(1-c)+c;  Q[1][2] = y*z*(1-c)-x*s;  Q[1][3] = 0.0;
	Q[2][0] = x*z*(1-c)-y*s; Q[2][1] = y*z*(1-c)+x*s; Q[2][2] = z*z*(1-c)+c;   Q[2][3] = 0.0;
	Q[3][0] = 0.0;  Q[3][1] = 0.0; Q[3][2] = 0.0; Q[3][3] = 1.0;

	CMatrix44 matrix = CMatrix44(Q[0][0],Q[0][1],Q[0][2],Q[0][3],
		Q[1][0],Q[1][1],Q[1][2],Q[1][3],
		Q[2][0],Q[2][1],Q[2][2],Q[2][3],
		Q[3][0],Q[3][1],Q[3][2],Q[3][3]);

	return matrix;
}


double saddle_face::return_rotation_angle(TorusIter current_torus,CVector3d previous_direction){
	
	CMatrix44 Q = current_torus->torus_to_para;

	CVector3d current_direction = CVector3d(0.0,0.0,0.0);

	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			current_direction[i] += previous_direction[j] * Q.Get(i,j);

	current_direction.Normalize();
	CVector3d X_direction = CVector3d(1.0,0.0,0.0);
	CVector3d rotation_X = current_direction.Cross(X_direction);
	double angle_X = acos(current_direction.Dot(X_direction));
	rotation_X.Normalize();

	if(rotation_X.Dot(CVector3d(0.0,0.0,1.0)) > 0)
		return -angle_X;
	else
		return angle_X;

}