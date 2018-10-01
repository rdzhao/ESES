#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h> 
#include "types.h"
#include "Torus.h"
#include "Atom_sphere.h"
#include "concave_edge.h"
#include "concave_sphere.h"
#include "saddle_face.h"

Torus::Torus() : torus_r(0.0)
{
}

Torus::~Torus(){
	mutual_atom_list.clear();
	concave_edge_list.clear();
	saddle_face_list.clear();
}

//initialize the torus: the adjacent atoms, torus parameters contact circles
Torus::Torus(AtomIter atomi,AtomIter atomj,double probe_radius){

	atom_i = atomi;  atom_j = atomj;
	
	CVector3d a_i = atom_i->center;   
	CVector3d a_j = atom_j->center;
	double r_i = atom_i->r;
	double r_j = atom_j->r;
	double d_ij =(a_i-a_j).Length();

	//compute torus center, axis and radius
	torus_center = (a_i + a_j)/2.0 + (a_j - a_i)/2.0*((r_i+probe_radius)*(r_i+probe_radius) - (r_j+probe_radius)*(r_j+probe_radius))/(d_ij*d_ij);
	torus_r = sqrt((r_i+r_j+2.0*probe_radius)*(r_i+r_j+2.0*probe_radius) - d_ij*d_ij) * sqrt(d_ij*d_ij - (r_i-r_j)*(r_i-r_j))/(2.0*d_ij);
	torus_axis = a_j -a_i;
	torus_axis.Normalize();
	is_free = true; //initially set every torus as free torus
	is_blocked = false;


	//compute contact circles
	circle_center_i = (r_i * torus_center + probe_radius * a_i)/(r_i+probe_radius);
	circle_center_j = (r_j * torus_center + probe_radius * a_j)/(r_j+probe_radius);
	double circle_radius_i = torus_r * r_i/(probe_radius+r_i);
	double circle_radius_j = torus_r * r_j/(probe_radius+r_j);
	atom_i->boundary_cicle_centers.push_back(circle_center_i);
	atom_i->boundary_circle_normals.push_back(-torus_axis);
	atom_j->boundary_cicle_centers.push_back(circle_center_j);
	atom_j->boundary_circle_normals.push_back(torus_axis);
	

	double circle_displacement_i = torus_axis.Dot(circle_center_i - a_i);
	double circle_displacement_j = torus_axis.Dot(circle_center_j - a_j);
	phi_i = atan(circle_displacement_i/circle_radius_i);
	phi_j = atan(circle_displacement_j/circle_radius_j);
	double right_phi_i = -phi_j;
	double right_phi_j = -phi_i;
	phi_i = right_phi_i;
	phi_j = right_phi_j;


	//compute rotation and transformation from the parametrized torus centered at origin
	CVector3d z_direction = CVector3d(0.0,0.0,1.0);
	rotation_axis = torus_axis.Cross(z_direction);
	rotation_angle = atan2(rotation_axis.Length(),torus_axis.Dot(z_direction));
	double x_val = torus_axis.Dot(z_direction);
	double y_val = rotation_axis.Length();
	double t_val = sqrt(x_val*x_val + y_val*y_val);
	double cos_val = x_val/t_val;
	double sin_val = -y_val/t_val;
	rotation_axis.Normalize();


	/*******************************
	compute the bounding box and transformation matrix 

	transform para to torus: first para to torus then translation
	transform torus to para: first translation then para
	***********************************/	 
	 
	para_to_torus = return_rotation_matrix(cos_val,sin_val,rotation_axis);;
	torus_to_para = para_to_torus.Transpose();
	
	/**********************
	compute the visibility sphere given the paper
	'Interactive Visualization of Molecular Surface Dynamics'
	we use the Eq(4)
	************************/
	CVector3d arbitrary = torus_axis;
	arbitrary[0] += rand()%10+1.0;
	arbitrary[1] += rand()%10+1.0;
	arbitrary[2] += rand()%10+1.0;
	CVector3d arb_tangent = arbitrary.Cross(torus_axis);
	arb_tangent.Normalize();
	CVector3d p = torus_center + torus_r * arb_tangent;
	CVector3d x_tmp = p - a_i; 
	x_tmp.Normalize();
	CVector3d x = r_i*x_tmp + a_i;
	vs_sphere_center = ((p-a_i).Length())/((p-a_i).Length() + (p-a_j).Length()) * (a_j - a_i) + a_i;
	vs_sphere_r = (vs_sphere_center - x).Length();


	////compute singular points
	//if(torus_r < probe_radius){
	//	CVector3d tmp_p1 = CVector3d(0,0,sqrt(probe_radius*probe_radius-torus_r*torus_r)/torus_r);
	//	CVector3d tmp_p2 = -tmp_p1;

	//	singluar_p1 = para_to_torus.MultMatVec(tmp_p1) + torus_center;
	//	singluar_p2 = para_to_torus.MultMatVec(tmp_p2) + torus_center;
	//	
	//}


}


CMatrix44 Torus::return_rotation_matrix(double c,double s, CVector3d rotation_axis){
	double x = rotation_axis[0];
	double y = rotation_axis[1];
	double z = rotation_axis[2];
	

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

double Torus::angle_between_vectors(CVector3d a, CVector3d b,CVector3d n) 
// rotate from a to b in anti-clock
{
	CVector3d axb = a.Cross(b);
	double angle = atan2( axb.Length(), a.Dot(b) );
	if (axb.Dot(n)< 0.0) angle = -angle;
	return angle;
}



void Torus::selectionSort(struct order_elem data[], int lenD)
{
	int j = 0;
	struct order_elem tmp;
	for(int i=0;i<lenD;i++){
		j = i;
		for(int k = i;k<lenD;k++){
			if(data[j].angle<data[k].angle){
				j = k;
			}
		}
		tmp = data[i];
		data[i] = data[j];
		data[j] = tmp;
	}
}


//remove the concave edge generated from the duplex concave face or the concave edges which are opposite of each other
void Torus::remove_duplicate_concave_edge(){
	if(concave_edge_list.size()==0)
		return;

	//first check concave edges that are opposite
	std::vector<concave_edgeIter> tmp_concave_edge;
	for(int i=0;i<concave_edge_list.size();i++){
		bool exsit_opposite = false;
		for(int j=0;j<concave_edge_list.size();j++){
			if((concave_edge_list[i]->pos1-concave_edge_list[j]->pos2).Length() < ZERO_TOL &&
				(concave_edge_list[i]->pos2-concave_edge_list[j]->pos1).Length() < ZERO_TOL){
					exsit_opposite = true;
					concave_edge_list[i]->current_concave->duplex_concave_faces.push_back(concave_edge_list[j]->current_concave);
					break;
			}
		}
		if(!exsit_opposite)
			tmp_concave_edge.push_back(concave_edge_list[i]);
	}

	concave_edge_list = tmp_concave_edge;
	tmp_concave_edge.clear();
	//then check the duplex concave edges
	for(int i=0;i<concave_edge_list.size();i++){
		bool exsit_duplex = false;
		for(int j=i+1;j<concave_edge_list.size();j++){
			if((concave_edge_list[i]->pos1-concave_edge_list[j]->pos1).Length() < ZERO_TOL &&
				(concave_edge_list[i]->pos2-concave_edge_list[j]->pos2).Length() < ZERO_TOL){
					exsit_duplex = true;
					concave_edge_list[i]->current_concave->duplex_concave_faces.push_back(concave_edge_list[j]->current_concave);
					concave_edge_list[j]->current_concave->duplex_concave_faces.push_back(concave_edge_list[i]->current_concave);
					break;
			}
		}
		if(!exsit_duplex)
			tmp_concave_edge.push_back(concave_edge_list[i]);
	}

	concave_edge_list.clear();
	concave_edge_list = tmp_concave_edge;
}

/**********************************
reorder the concave edges such that it is clockwise
when viewing from atom i to atom j.

We tackle this problem by computing the
vertices on the contact circle of atom i,
and order the vertices such that
they are in the right order (angle in decreasing order)

***********************************/
void Torus::reorder_concave_edge(){

	remove_duplicate_concave_edge();

	int num_edge = concave_edge_list.size();

	if(num_edge==0)
		return;

	struct order_elem *compare_list = new struct order_elem[num_edge];

	CVector3d standard_direction;
	
	for(int i=0;i<num_edge;i++){
		CVector3d contact_pos_i;
		if(concave_edge_list[i]->atom1 == atom_i->index)
			contact_pos_i = concave_edge_list[i]->pos1;
		else
			contact_pos_i = concave_edge_list[i]->pos2;

		CVector3d current_direction = contact_pos_i - circle_center_i;
		current_direction.Normalize();

		if(i==0)
			standard_direction = current_direction;
		
		compare_list[i].angle = angle_between_vectors(current_direction,standard_direction,torus_axis);
		compare_list[i].index = i;
	}

	selectionSort(compare_list,num_edge);

	//save
	std::vector<concave_edgeIter>new_list;
	new_list = concave_edge_list;

	for(int i=0; i<new_list.size(); i++)
		new_list[i] = concave_edge_list[compare_list[i].index];

	concave_edge_list = new_list;
	new_list.clear();

	delete []compare_list;


}


void Torus::generate_true_self_intersecting_boundary(double probe_radius){
	if(is_free || is_blocked)
		return;
	if(torus_r >= probe_radius)
		return;
}

