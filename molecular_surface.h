/******************************
April/16/2014 

This header file is for computing the molecular surface based on paper

"Analytical Molecular Surface Calculation" 
   By Michael L. Connolly 1983

Format for Input:
//*.xyzr file containing the atom information
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
#include "Vector3d.h"
#include "Atom_sphere.h"
#include "Torus.h"
#include "concave_edge.h"
#include "concave_sphere.h"
#include "saddle_face.h"
#include "grid_point.h"
#include "grid_edge.h"
#include "types.h"

#include "Jenkins_Traub.h"


#undef max
#undef min
#define  TEST 0
#define  NEIGHBOR_LENGTH 50
//#define BLOCK_SIZE 127

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


	molecular_surface(){	
	}

	~molecular_surface(){
	}


CVector3d grid_pos(int X,int Y,int Z){
	return CVector3d(X * m_x_step + m_min_x,Y * m_y_step + m_min_y,Z * m_z_step + m_min_z);
}

//return the global index of the grid point
int grid_index(int X,int Y, int Z){
	return (X + m_x_num*Y + m_x_num*m_y_num*Z);
}

//construct the grid points and grid edge
void initialize_grid(){
	m_grid_points.resize(m_x_num*m_y_num*m_z_num);
	//construct the grid edge
	for(int i=0; i<m_x_num; i++)
		for(int j=0; j<m_y_num; j++)
			for(int k=0; k<m_z_num; k++){
				grid_point*  new_point = new grid_point(i,j,k);
				new_point->initialize_grid_point(m_x_step,m_y_step,m_z_step,m_min_x,m_min_y,m_min_z);
				//std::cout << new_point->grid_pos.x()<<" "<< new_point->grid_pos.y() << " "<<new_point->grid_pos.z() <<std::endl;
				m_grid_points[grid_index(i,j,k)] = new_point;
			}

			for(int i=0;i<m_x_num;i++)
				for(int j=0; j<m_y_num; j++)
					for(int k=0;k<m_z_num; k++){
						grid_pointIter point_1 = m_grid_points[grid_index(i,j,k)];

						if(i<m_x_num-1){
							grid_pointIter point_2 = m_grid_points[grid_index(i+1,j,k)];
							grid_edgeIter new_edge = new grid_edge(point_1,point_2);
							m_grid_edges.push_back(new_edge);
						}
						if(j<m_y_num-1){
							grid_pointIter point_2 = m_grid_points[grid_index(i,j+1,k)];
							grid_edgeIter new_edge = new grid_edge(point_1,point_2);
							m_grid_edges.push_back(new_edge);

						}

						if(k<m_z_num-1){
							grid_pointIter point_2 = m_grid_points[grid_index(i,j,k+1)];
							grid_edgeIter new_edge = new grid_edge(point_1,point_2);
							m_grid_edges.push_back(new_edge);
						}
					}

}


/**********************constructing the molecular surface starts**********************************/
/**********************constructing the molecular surface starts**********************************/
/**********************constructing the molecular surface starts**********************************/

//initialize the molecular surface
void initialize_molecular_surface(){
	construct_torus();   //construct torus
	cout<<"construct torus complete.."<<endl;
	mark_torus(); //mark the torus as free and non-free, blocked or not
	cout<<"mark torus complete.."<<endl;
	probe_placement(); //construct the concave face and concave edge
	cout<<"probe placement complete.."<<endl;
	check_duplex_concave();
	cout<<"check duplex complete.."<<endl;
	construct_saddle_faces();
	cout<<"construct saddle complete.."<<endl;
	
}



/*********************
construct the saddle face by choosing
the right concave edges and convex edges

(meaning the current concave edge should 
pointing from vi to vj, while the next concave edge
should pointing from vj to vi)
*********************/
void construct_saddle_each_torus(TorusIter current_torus){

	int nsize = (current_torus->concave_edge_list).size();

	for(int i=0; i<nsize;i++){
		concave_edgeIter concave_e = current_torus->concave_edge_list[i];
		concave_edgeIter concave_next = current_torus->concave_edge_list[(i+1)%nsize];
		if(concave_e->atom1 == current_torus->atom_i->index && concave_next->atom1 == current_torus->atom_j->index){
			//first construct the convex edge
			saddle_faceIter current_saddle = new saddle_face(concave_e,concave_next,current_torus,m_probe_radius);
			m_saddle_faces.push_back(current_saddle);
			
			(current_torus->saddle_face_list).push_back(current_saddle);
		}
	}

}

void construct_saddle_faces(){
	for(int i=0;i<m_torus_vec.size();i++){
		TorusIter current_torus = m_torus_vec[i];
		//if this torus has intersections (free torus is not included)
		if(current_torus->concave_edge_list.size() >0){
		
			construct_saddle_each_torus(current_torus);
			continue;
		}
	}
}

//fix the duplex concave spheres (when the probe is touching more than three atoms at the same time).
void check_duplex_concave(){
	for (int i=0;i<m_torus_vec.size();i++){
		TorusIter current_torus = m_torus_vec[i];
		current_torus->mutual_atom_list.clear();
		if(!current_torus->is_free && current_torus->concave_edge_list.size()==0){
			current_torus->is_blocked = true;
			continue;
		}

		if(current_torus->concave_edge_list.size()==1)
			continue;

		if(current_torus->concave_edge_list.size()>0){
			current_torus->reorder_concave_edge();
		}
	}

	std::vector<Concave_SphereIter> tmp_duplex_concaves;
	/*std::vector<Concave_SphereIter> gloabal_store;
	std::vector<std::vector<int>> global_duplex_ids;*/
	for(int s_id=0;s_id<m_concave_spheres.size();s_id++){
		Concave_SphereIter current_concave = m_concave_spheres[s_id];
		if(current_concave->duplex_concave_faces.size()==0)
			continue;
		if(current_concave->has_check_duplex)
			continue;

		tmp_duplex_concaves.clear();
		store_duplex_concaves(current_concave,tmp_duplex_concaves);
		std::vector<int> duplex_atom_ids;
		for(int i=0;i<tmp_duplex_concaves.size();i++){
			tmp_duplex_concaves[i]->duplex_concave_faces = tmp_duplex_concaves;
			tmp_duplex_concaves[i]->has_check_duplex = true;
		}

	

	}
}

void store_duplex_concaves(Concave_SphereIter current_concave,std::vector<Concave_SphereIter> &full_duplex_concaves){
	for(int i=0;i<current_concave->duplex_concave_faces.size();i++){
		Concave_SphereIter duplex_concave = current_concave->duplex_concave_faces[i];
		bool is_checked = false;
		for(int k=0;k<full_duplex_concaves.size();k++){
			if(full_duplex_concaves[k]==duplex_concave){
				is_checked = true;
				break;
			}
		}

		if(is_checked)
			continue;
		full_duplex_concaves.push_back(duplex_concave);
		store_duplex_concaves(duplex_concave,full_duplex_concaves);
	}
}

//construct the torus
void construct_torus(){
	//cout << "Checking neighbour info.." << endl;
	//first check the neighboring information of the atoms
	for(AtomPointerIter atom_i=m_atoms_vec.begin(); atom_i!=m_atoms_vec.end();atom_i++){
		(*atom_i)->initialize_atom(m_N_atoms);
		AtomPointerIter atom_j = atom_i;
		atom_j++;
		for(; atom_j!=m_atoms_vec.end(); atom_j++){
			if(is_neighbor_atom(*atom_i,*atom_j)){
				(*atom_i)->neighbor_list.push_back(*atom_j);
				(*atom_j)->neighbor_list.push_back(*atom_i);
			}
		}
	}
	
	//cout << "Constructing torus.." << endl;
	//construct the torus
	for(AtomPointerIter atom_i=m_atoms_vec.begin(); atom_i!=m_atoms_vec.end();atom_i++){
		//inside torus
		if((*atom_i)->is_inside)
			continue;
		for(int i=0;i<(*atom_i)->neighbor_list.size();i++){
			AtomIter atom_j = (*atom_i)->neighbor_list[i];
			if(atom_j->is_inside)
				continue;
			//avoid duplex torus
			if((*atom_i)->index > atom_j->index)
				continue;

			Torus* torus_ij = new Torus((*atom_i),atom_j,m_probe_radius);
			m_torus_vec.push_back(torus_ij);
		}
	}
	//cout << "Initalize torus relation with atoms.." << endl;
	//initialize the torus mutual ngb
	for(int t_id=0;t_id<m_torus_vec.size();t_id++){
		TorusIter current_torus = m_torus_vec[t_id];
		//find the common index in the neighbors of atom i,j
		int index_i = current_torus->atom_i->index;
		int index_j = current_torus->atom_j->index;

		AtomIter atom_i = current_torus->atom_i;
		AtomIter atom_j = current_torus->atom_j;

		atom_i->adj_torus_list[index_j] = current_torus;
		atom_j->adj_torus_list[index_i] = current_torus;

		for(int m=0; m<atom_i->neighbor_list.size(); m++){
			if(atom_i->neighbor_list[m]->is_inside)
				continue;
			int current_index = atom_i->neighbor_list[m]->index;
			for(int l=0; l<atom_j->neighbor_list.size(); l++)
				if(current_index== atom_j->neighbor_list[l]->index)
					current_torus->mutual_atom_list.push_back(atom_i->neighbor_list[m]);
		}
	}

}


//check whether atom_i and atom_j are neighboring atom (i.e, their distance  fabs(r_j-r_i)< dij < r_i + r_j + 2*r_p)
bool is_neighbor_atom(AtomIter atom_i, AtomIter atom_j){
	double r_i = atom_i->r;
	double r_j = atom_j->r;
	double atom_distance = (atom_i->center - atom_j->center).Length();

	if(atom_distance > r_i+r_j+2.0*m_probe_radius)
		return false;
	if(atom_distance < fabs(r_j - r_i)){
		if(r_i > r_j)
			atom_j->is_inside = true;
		else
			atom_i->is_inside = true;

		return false;
	}
	else
		return true;

}

/*********************************************
	mark if it is a free torus or not.  A torus (t_ij) is free
	when for EVERY mutual neighbor atom k of atom i,j,
	there is no intersection with atom k

********************************************/
void mark_torus(){
	for(int t_id=0;t_id<m_torus_vec.size();t_id++){
		TorusIter current_torus = m_torus_vec[t_id];
		//already marked as non-free
		if(!current_torus->is_free)
			continue;
		
		CVector3d u_ij = current_torus->torus_axis;
		CVector3d t_ij = current_torus->torus_center;
		double r_ij = current_torus->torus_r;
		int i = current_torus->atom_i->index;
		int j = current_torus->atom_j->index;
		CVector3d a_i = current_torus->atom_i->center;
		double r_i = current_torus->atom_i->r;
		for(int l=0; l<current_torus->mutual_atom_list.size(); l++){
			AtomIter atom_k = current_torus->mutual_atom_list[l];
			CVector3d a_k = atom_k->center;
			CVector3d u_ik = a_k - a_i;
			double r_k = atom_k->r;
			u_ik.Normalize();

			double tmp_value = u_ij.Dot(u_ik);
			double w_ijk = acos(tmp_value);


			//check whether the center of atom i,j,k are colinear
			if(fabs(w_ijk)<ZERO_TOL || fabs(w_ijk - PI)<ZERO_TOL)
				continue;   

			CVector3d u_ijk = u_ij.Cross(u_ik);
			u_ijk.Normalize();
			CVector3d u_tb = u_ijk.Cross(u_ij);
			CVector3d t_ik = (atom_k->adj_torus_list[i])->torus_center;
			CVector3d b_ijk = t_ij + u_tb*(u_ik.Dot(t_ik - t_ij))/sqrt(1-tmp_value*tmp_value);

			double check_square = (r_i + m_probe_radius)*(r_i + m_probe_radius) - (b_ijk - a_i).LengthSquared();

			if(check_square < 0){
				/****************************
				means there is no probe placement for atom i,j,k; 
				now we need to check whether it is because
				(1) torus (t_ij) circle lie entirely outside sphere atom k---no intersection

				OR

				(2) the torus lies entirely inside the sphere atom k --- not accessible

				use the Eq (A11)
				********************************/
				double right_val = sqrt((r_k + m_probe_radius)*(r_k + m_probe_radius) - r_ij*r_ij );
				double left_val = (t_ij - a_k).Length();

				if(left_val < right_val){
					//the torus lies entirely in the sphere atom k
					current_torus->is_free = false; 
					current_torus->is_blocked = true;
					break;
				}
			}

			else{
				//there is possible probe_placement, then the torus_ij,torus_jk, torus_ik is not free
				current_torus->is_free = false;
				(current_torus->atom_i->adj_torus_list[atom_k->index])->is_free = false; //t_ik
				(atom_k->adj_torus_list[j])->is_free = false; //t_jk
			}
		}

	}
}


/***************************************************
collision detection when constructing concave sphere
probe_1, probe_2 : possible location of probe sphere center
i,j,k : index for atom i,j,k
is_collide_1, is_collide_2 : return the collision detection result

Check whether the location of probe sphere will collide with the mutual neighbors of atom i, j,k

***************************************************/
void collision_detection(CVector3d probe_1,CVector3d probe_2,AtomIter atom_i,AtomIter atom_j, AtomIter atom_k,bool &is_collide_1, bool &is_collide_2){
	std::vector<AtomIter> mutual_ngb_ijk; //store the pointer of mutual neighbors of atom i,j,k
	TorusIter torus_ij = atom_i->adj_torus_list[atom_j->index];
	TorusIter torus_jk = atom_j->adj_torus_list[atom_k->index];

	for(int l=0; l<torus_ij->mutual_atom_list.size(); l++){
		int current_index = torus_ij->mutual_atom_list[l]->index;
		for(int m=0; m<torus_jk->mutual_atom_list.size(); m++)
			if(current_index == torus_jk->mutual_atom_list[m]->index)
				mutual_ngb_ijk.push_back(torus_ij->mutual_atom_list[l]);
	}

	is_collide_1 = is_collide_2 = false;

	for(int l=0; l<mutual_ngb_ijk.size(); l++){
		AtomIter tmp_atom = mutual_ngb_ijk[l];
		CVector3d atom_center = tmp_atom->center;
		double distance_1 = (probe_1 - atom_center).Length();
		double distance_2 = (probe_2 - atom_center).Length();
		double comp_r = tmp_atom->r+m_probe_radius-ZERO_TOL;

		if(!is_collide_1 && (distance_1 < comp_r))
			is_collide_1 = true;
		
		if(!is_collide_2 && (distance_2 < comp_r))
			is_collide_2 = true;

		if(is_collide_1 && is_collide_2)
			return;
		
	}


}

void construct_single_concave(CVector3d probe_center,AtomIter atom_i,AtomIter atom_j,AtomIter atom_k,TorusIter torus_ij,TorusIter torus_jk, TorusIter torus_ki){
	double r_i  = atom_i->r;
	double r_j  = atom_j->r;
	double r_k = atom_k->r;

	CVector3d a_i = atom_i->center;
	CVector3d a_j = atom_j->center;
	CVector3d a_k = atom_k->center;

	CVector3d vertex_pi = (r_i*probe_center + m_probe_radius*a_i)/(r_i + m_probe_radius);
	CVector3d vertex_pj = (r_j*probe_center + m_probe_radius*a_j)/(r_j + m_probe_radius);
	CVector3d vertex_pk= (r_k*probe_center+ m_probe_radius*a_k)/(r_k + m_probe_radius);

	//make sure the concave edge has counter-clockwise orientation from the center of the probe sphere
	CVector3d tmp = (vertex_pi + vertex_pj + vertex_pk)/3.0;
	CVector3d tmp_center = probe_center - tmp;
	CVector3d e_ij = vertex_pj - vertex_pi; CVector3d e_jk = vertex_pk - vertex_pj;
	CVector3d tmp_norm = e_ij.Cross(e_jk);

	if(tmp_norm.Dot(tmp_center) > 0){
		concave_edgeIter edge_ij  = new concave_edge(vertex_pi,vertex_pj,atom_i->index,atom_j->index);
		concave_edgeIter edge_jk = new concave_edge(vertex_pj,vertex_pk,atom_j->index,atom_k->index);
		concave_edgeIter edge_ki = new concave_edge(vertex_pk,vertex_pi,atom_k->index,atom_i->index);
		Concave_SphereIter current_concave = new Concave_Sphere(edge_ij,edge_jk,edge_ki,atom_i->index,atom_j->index,atom_k->index,a_i,a_j,a_k,
			atom_i->r,atom_j->r,atom_k->r,probe_center);

		m_concave_edges.push_back(edge_ij);
		m_concave_edges.push_back(edge_jk);
		m_concave_edges.push_back(edge_ki);
		m_concave_spheres.push_back(current_concave);

		edge_ij->initial_concave_edge(probe_center,torus_ij->torus_center,torus_ij->torus_axis);
		edge_jk->initial_concave_edge(probe_center,torus_jk->torus_center,torus_jk->torus_axis);
		edge_ki->initial_concave_edge(probe_center,torus_ki->torus_center,torus_ki->torus_axis);

		//store the concave edges
		torus_ij->concave_edge_list.push_back(edge_ij);
		torus_jk->concave_edge_list.push_back(edge_jk);
		torus_ki->concave_edge_list.push_back(edge_ki);

		edge_ij->current_concave = current_concave;
		edge_jk->current_concave = current_concave;
		edge_ki->current_concave = current_concave;
	}
	else{
		concave_edgeIter edge_ji  = new concave_edge(vertex_pj,vertex_pi,atom_j->index,atom_i->index);
		concave_edgeIter edge_kj = new concave_edge(vertex_pk,vertex_pj,atom_k->index,atom_j->index);
		concave_edgeIter edge_ik = new concave_edge(vertex_pi,vertex_pk,atom_i->index,atom_k->index); 
		Concave_SphereIter current_concave = new Concave_Sphere(edge_ji,edge_kj,edge_ik,atom_j->index,atom_k->index,atom_i->index,a_j,a_k,a_i,
			atom_j->r,atom_k->r,atom_i->r,probe_center);

		m_concave_edges.push_back(edge_ji);
		m_concave_edges.push_back(edge_kj);
		m_concave_edges.push_back(edge_ik);
		m_concave_spheres.push_back(current_concave);

		edge_ji->current_concave = current_concave;
		edge_ik->current_concave = current_concave;
		edge_kj->current_concave = current_concave;
	
		edge_ji->initial_concave_edge(probe_center,torus_ij->torus_center,torus_ij->torus_axis);
		edge_kj->initial_concave_edge(probe_center,torus_jk->torus_center,torus_jk->torus_axis);
		edge_ik->initial_concave_edge(probe_center,torus_ki->torus_center,torus_ki->torus_axis);

		torus_ij->concave_edge_list.push_back(edge_ji);
		torus_jk->concave_edge_list.push_back(edge_kj);
		torus_ki->concave_edge_list.push_back(edge_ik);
	}

}


void probe_placement(){
	for(int t_id=0;t_id<m_torus_vec.size();t_id++){
		TorusIter current_torus = m_torus_vec[t_id];
		if(current_torus->is_free || current_torus->is_blocked)
			continue;   //no concave sphere for torus which is free or blocked

		int i = current_torus->atom_i->index;
		int j = current_torus->atom_j->index;
		AtomIter atom_i = current_torus->atom_i;
		AtomIter atom_j = current_torus->atom_j;
		CVector3d u_ij = current_torus->torus_axis;
		CVector3d t_ij = current_torus->torus_center;
		double r_ij = current_torus->torus_r;
		CVector3d a_i =atom_i->center;
		double r_i = atom_i->r;


		for(int l=0; l<current_torus->mutual_atom_list.size(); l++){
			AtomIter atom_k = current_torus->mutual_atom_list[l];
			if((atom_k->adj_torus_list[i]->is_free || atom_k->adj_torus_list[i]->is_blocked) &&
				(atom_k->adj_torus_list[j]->is_free || atom_k->adj_torus_list[j]->is_blocked))
				continue;

			int k = atom_k->index;
			if(k < j || k < i)
				continue;

			
			CVector3d a_k = atom_k->center;
			CVector3d u_ik = a_k - a_i;
			double r_k = atom_k->r;
			u_ik.Normalize();

			double tmp_value = u_ij.Dot(u_ik);
			double w_ijk = acos(tmp_value);

			//check whether the center of atom i,j,k are colinear
			if(fabs(w_ijk)<ZERO_TOL || fabs(w_ijk - PI)<ZERO_TOL)
				continue;   

			CVector3d u_ijk = u_ij.Cross(u_ik);
			u_ijk.Normalize();
			CVector3d u_tb = u_ijk.Cross(u_ij);
			CVector3d t_ik = (atom_k->adj_torus_list[i])->torus_center;
			CVector3d b_ijk = t_ij + u_tb*(u_ik.Dot(t_ik - t_ij))/sqrt(1-tmp_value*tmp_value);

			double check_square = (r_i + m_probe_radius)*(r_i + m_probe_radius) - (b_ijk - a_i).LengthSquared();

		
			if(check_square>0){
				double h_ijk = sqrt(check_square);
				CVector3d probe_1 = b_ijk + h_ijk * u_ijk;
				CVector3d probe_2 = b_ijk - h_ijk * u_ijk;



				bool is_collide1,is_collide2;
				//collision detect
				collision_detection(probe_1,probe_2,atom_i,atom_j,atom_k,is_collide1,is_collide2);

				if(!is_collide1)
					construct_single_concave(probe_1,atom_i,atom_j,atom_k,current_torus,(atom_j->adj_torus_list[k]), (atom_k->adj_torus_list[i]));
				if(!is_collide2)
					construct_single_concave(probe_2,atom_i,atom_j,atom_k,current_torus,(atom_j->adj_torus_list[k]), (atom_k->adj_torus_list[i]));
			}
		}
	}


	for(int i=0;i<m_atoms_vec.size();i++)
		m_atoms_vec[i]->adj_torus_list.clear();
}

/**********************constructing the molecular surface ends**********************************/
/**********************constructing the molecular surface ends**********************************/
/**********************constructing the molecular surface ends**********************************/


//check whether current gird point is the atom sphere or not
bool is_inside_atom(AtomIter current_atom,CVector3d grid_loc){
	double radius = current_atom->r;
	double tmp = (grid_loc - current_atom->center).LengthSquared() - radius*radius;

	if(tmp>0)
		return false;
	else
		return true;

}

bool is_inside_concave_sphere(CVector3d grid_pos,CVector3d sphere_center){
	double tmp = (grid_pos-sphere_center).LengthSquared()-pow(m_probe_radius,2);
	if(tmp>0)
		return false;
	else
		return true;

}

/********************************************************************************
check whether test_point is inside the tetrahedron composed by p1, p2, p3, p4

Here we use the method from :
http://steve.hollasch.net/cgindex/geometry/ptintet.html
**********************************************************************************/
bool is_inside_tetrahedron(CVector3d p1,CVector3d p2, CVector3d p3, CVector3d p4, CVector3d test_point){
	CMatrix44 D0 = CMatrix44(p1[0],p1[1],p1[2],1.0,p2[0],p2[1],p2[2],1.0,
		p3[0],p3[1],p3[2],1.0,p4[0],p4[1],p4[2],1.0);

	CMatrix44 D1 = CMatrix44(test_point[0],test_point[1],test_point[2],1.0,p2[0],p2[1],p2[2],1.0,
		p3[0],p3[1],p3[2],1.0,p4[0],p4[1],p4[2],1.0);

	CMatrix44 D2 = CMatrix44(p1[0],p1[1],p1[2],1.0,test_point[0],test_point[1],test_point[2],1.0,
		p3[0],p3[1],p3[2],1.0,p4[0],p4[1],p4[2],1.0);

	CMatrix44 D3 = CMatrix44(p1[0],p1[1],p1[2],1.0,p2[0],p2[1],p2[2],1.0,
		test_point[0],test_point[1],test_point[2],1.0,p4[0],p4[1],p4[2],1.0);

	CMatrix44 D4 = CMatrix44(p1[0],p1[1],p1[2],1.0,p2[0],p2[1],p2[2],1.0,
		p3[0],p3[1],p3[2],1.0,test_point[0],test_point[1],test_point[2],1.0);

	double check_D0 = D0.Determinant();
	double check_D1 = D1.Determinant();
	double check_D2 = D2.Determinant();
	double check_D3 = D3.Determinant();
	double check_D4 = D4.Determinant();

	if(check_D0>=0 && check_D1>=0 && check_D2>=0 && check_D3>=0 && check_D4>=0)
		return true;
	if(check_D0<=-0 && check_D1<=-0 && check_D2<=-0 && check_D3<=-0 && check_D4<=-0)
		return true;

	return false;
}

bool is_f_InsideTorus_positive(double x,double y,double z,double R,double r){
	double tmp = R - sqrt(x*x+y*y);
	double f = tmp*tmp + z*z - r*r;
	if(f>0)
		return true;
	else
		return false;
}

bool is_f_quartic_equation_positive(double x,double y,double z,double R,double r){
	double f = 4*R*R*(z*z-r*r) + (pow(r*r+R*R-x*x-y*y-z*z,2));
	if(f>0)
		return true;
	else
		return false;
}

bool is_located_in_selfintersecting_part(double x,double y,double z,double R,double r){
	
	if((z*z <= r*r-R*R) && 4*R*R*(r*r-z*z) <= (pow(r*r+R*R-x*x-y*y-z*z,2)))
		return true;
	else
		return false;
}

bool is_located_within_lemon_part(double x,double y,double z,double R,double r){
	//current torus is a donut torus
	if(R>r)
		return false;

	double m = sqrt(x*x+y*y);
	if(m<=r-R && pow(z,2)<=(r*r-pow(m+R,2)))
		return true;
	else
		return false;

	//if((z*z <= r*r-R*R) && 4*R*R*(r*r-z*z) <= (pow(r*r+R*R-x*x-y*y-z*z,2)+1e-6))
	//	return true;
	//else
	//	return false;
}

bool is_located_on_lemon_part(double x,double y,double z,double R,double r){
	//current torus is a donut torus
	if(R>r)
		return false;

	double m = sqrt(x*x+y*y);
	if(m<=r-R+1e-6 && pow(z,2)<=(r*r-pow(m+R,2))+1e-6)
		return true;
	else
		return false;

	//if((z*z <= r*r-R*R) && 4*R*R*(r*r-z*z) <= (pow(r*r+R*R-x*x-y*y-z*z,2)+1e-6))
	//	return true;
	//else
	//	return false;
}

bool is_covered_by_saddles(CVector3d grid_point,TorusIter current_torus){

	//now check whether it is located in the saddle faces part
	CVector3d current_x = grid_point - current_torus->circle_center_i;
	current_x = current_x - (current_x.Dot(current_torus->torus_axis))*current_torus->torus_axis;
	current_x.Normalize();

	for(int i=0;i<current_torus->saddle_face_list.size();i++){
		saddle_faceIter current_saddle = current_torus->saddle_face_list[i];
		bool is_between = false;

		if((current_saddle->boundary_i_1.Cross(current_x)).Dot(current_saddle->boundary_cross)>0 &&
			(current_x.Cross(current_saddle->boundary_i_2)).Dot(current_saddle->boundary_cross)>0 )
			is_between = true;

		if(current_saddle->boundary_check && is_between)
			return true;
		if(!current_saddle->boundary_check && !is_between)
			return true;
	}

	return false;
}

bool is_inside_Visibility_Sphere_old(TorusIter current_torus,CVector3d point_pos){
	double check = (point_pos - current_torus->vs_sphere_center).LengthSquared() - pow(current_torus->vs_sphere_r,2);
	
	if(check>0)
		return false;

	if(current_torus->is_blocked){
		CVector3d para_pos = (current_torus->torus_to_para).MultMatVec(point_pos - current_torus->torus_center);
		double x = para_pos[0];
		double y = para_pos[1];
		double z = para_pos[2];

		if(is_f_InsideTorus_positive(x,y,z,current_torus->torus_r,m_probe_radius))
			return true;
		else
			return false;
	}

	if(current_torus->is_free || is_covered_by_saddles(point_pos,current_torus))
		return true;
	else
		return false;
}

bool is_inside_Visibility_Sphere(TorusIter current_torus,CVector3d point_pos){
	double check = (point_pos - current_torus->vs_sphere_center).LengthSquared() - pow(current_torus->vs_sphere_r,2);

	if(check>0)
		return false;

	if(current_torus->is_blocked){
		CVector3d para_pos = (current_torus->torus_to_para).MultMatVec(point_pos - current_torus->torus_center);
		double x = para_pos[0];
		double y = para_pos[1];
		double z = para_pos[2];

		if(!is_f_quartic_equation_positive(x,y,z,current_torus->torus_r,m_probe_radius))
			return false;
		else{
			if(is_located_within_lemon_part(x,y,z,current_torus->torus_r,m_probe_radius))
				return false;
			else
				return true;
		}
	}

	if(current_torus->is_free || is_covered_by_saddles(point_pos,current_torus))
		return true;
	else
		return false;
}


bool validate_atom_intersection(AtomIter atom,CVector3d pos){
	for(int k=0;k<atom->boundary_cicle_centers.size(); k++){
		CVector3d current_center = atom->boundary_cicle_centers[k];
		CVector3d current_normal = atom->boundary_circle_normals[k];
		if((pos-current_center).Dot(current_normal) < 0 )
			return false;
	}


	//check it is not inside the neighboring atoms
	for(int m=0;m<atom->neighbor_list.size();m++){
		AtomIter neighbor_atom = atom->neighbor_list[m];
		if(is_inside_atom(neighbor_atom,pos))
			return false;
	}

	return true;

}

/**********************labeling the grid points starts**********************************/
/**********************labeling the grid points starts**********************************/
/**********************labeling the grid points starts**********************************/



void process_each_atom_newer(AtomIter current_atom){
	if(current_atom->is_inside)
		return;

	int tmp_min_x = int(floor((current_atom->center[0] - current_atom->r - m_probe_radius - m_min_x)/m_x_step));
	int tmp_min_y = int(floor((current_atom->center[1] - current_atom->r - m_probe_radius - m_min_y)/m_y_step));
	int tmp_min_z = int(floor((current_atom->center[2] - current_atom->r - m_probe_radius - m_min_z)/m_z_step));

	int tmp_max_x = int(ceil((current_atom->center[0] + current_atom->r + m_probe_radius - m_min_x)/m_x_step));
	int tmp_max_y = int(ceil((current_atom->center[1] + current_atom->r + m_probe_radius - m_min_y)/m_y_step));
	int tmp_max_z = int(ceil((current_atom->center[2] + current_atom->r + m_probe_radius - m_min_z)/m_z_step));

	//sanity check
	int min_x_id = (tmp_min_x>0 ? tmp_min_x:0);
	int min_y_id = (tmp_min_y>0 ? tmp_min_y:0);
	int min_z_id = (tmp_min_z>0 ? tmp_min_z:0);

	int max_x_id = (tmp_max_x<m_x_num-1 ? tmp_max_x : m_x_num-1);
	int max_y_id = (tmp_max_y<m_y_num-1 ? tmp_max_y : m_y_num-1);
	int max_z_id = (tmp_max_z<m_z_num-1 ? tmp_max_z : m_z_num-1);
	for(int x_i=min_x_id; x_i<=max_x_id; x_i++)
		for(int y_i=min_y_id; y_i<=max_y_id; y_i++)
			for(int z_i=min_z_id; z_i<=max_z_id; z_i++){
				grid_pointIter current_point = m_grid_points[grid_index(x_i,y_i,z_i)];
				double dist_square = (current_point->grid_pos - current_atom->center).LengthSquared();

				if(current_point->is_outside_aug_atom){
					double check_distance_aug = dist_square - (current_atom->r+m_probe_radius)*(current_atom->r+m_probe_radius);
					if(check_distance_aug > 0)
						continue;
					current_point->is_outside_aug_atom = false;  //it is located within the augmented sphere
				}

				//stores the possible atoms which include current grid point
				double check_distance_possible = dist_square - (current_atom->r+m_grid_resolution)*(current_atom->r+m_grid_resolution);
				if(check_distance_possible > 0)
					continue;
				current_point->possible_atom_interior.push_back(current_atom);
				
				if(current_point->is_processed)
					continue;

				double check_distance_interior = dist_square - (current_atom->r)*(current_atom->r);
				if(check_distance_interior > 0)
					continue;				
			//	current_point->is_inside_atom = true;
				current_point->is_processed = true;
				current_point->is_inside = true;
			}
}


void process_each_concave_sphere_newer(Concave_SphereIter current_concave){

	int tmp_min_x = int(floor((current_concave->probe_center[0]-m_probe_radius-m_min_x)/m_x_step));
	int tmp_min_y = int(floor((current_concave->probe_center[1]-m_probe_radius-m_min_y)/m_y_step));
	int tmp_min_z = int(floor((current_concave->probe_center[2]-m_probe_radius-m_min_z)/m_z_step));

	int tmp_max_x = int(ceil((current_concave->probe_center[0]+m_probe_radius-m_min_x)/m_x_step));
	int tmp_max_y = int(ceil((current_concave->probe_center[1]+m_probe_radius-m_min_y)/m_y_step));
	int tmp_max_z = int(ceil((current_concave->probe_center[2]+m_probe_radius-m_min_z)/m_z_step));

	//sanity check
	int min_x_id = (tmp_min_x>0 ? tmp_min_x:0);
	int min_y_id = (tmp_min_y>0 ? tmp_min_y:0);
	int min_z_id = (tmp_min_z>0 ? tmp_min_z:0);

	int max_x_id = (tmp_max_x<m_x_num-1 ? tmp_max_x : m_x_num-1);
	int max_y_id = (tmp_max_y<m_y_num-1 ? tmp_max_y : m_y_num-1);
	int max_z_id = (tmp_max_z<m_z_num-1 ? tmp_max_z : m_z_num-1);

	//first check the neighboring grid points and also reset the tmp mark for grid edges
	for(int x_i=min_x_id; x_i<=max_x_id; x_i++)
		for(int y_i=min_y_id; y_i<=max_y_id; y_i++)
			for(int z_i=min_z_id; z_i<=max_z_id; z_i++){
				grid_pointIter current_point = m_grid_points[grid_index(x_i,y_i,z_i)];
				current_point->possible_concavesphere_interior.push_back(current_concave);

				if(current_point->is_processed)
					continue;

				//check whether it is inside the concave face or not
				if(is_inside_concave_sphere(current_point->grid_pos,current_concave->probe_center)){
				//	current_point->is_inside_concave_sphere = true;
					current_point->is_processed = true;
				}
				

				if(is_inside_tetrahedron(current_concave->center_i,current_concave->center_j,current_concave->center_k,current_concave->probe_center,current_point->grid_pos))
					current_point->is_inside_VS_or_dualTet = true;
				
			}
}

void process_each_torus_newer(TorusIter current_torus){

	AtomIter atom_i = current_torus->atom_i;
	AtomIter atom_j = current_torus->atom_j;


	//decide which gird points are inside this atom sphere
	int tmp_min_x = int(floor((std::min((atom_i->center[0]-atom_i->r), (atom_j->center[0]-atom_j->r))- m_probe_radius - m_min_x)/m_x_step));
	int tmp_min_y = int(floor((std::min((atom_i->center[1]-atom_i->r), (atom_j->center[1]-atom_j->r))- m_probe_radius - m_min_y)/m_y_step));
	int tmp_min_z = int(floor((std::min((atom_i->center[2]-atom_i->r), (atom_j->center[2]-atom_j->r))- m_probe_radius - m_min_z)/m_z_step));

	int tmp_max_x = int(ceil((std::max((atom_i->center[0]+atom_i->r), (atom_j->center[0]+atom_j->r))+ m_probe_radius - m_min_x)/m_x_step));
	int tmp_max_y = int(ceil((std::max((atom_i->center[1]+atom_i->r), (atom_j->center[1]+atom_j->r))+ m_probe_radius - m_min_y)/m_y_step));
	int tmp_max_z = int(ceil((std::max((atom_i->center[2]+atom_i->r), (atom_j->center[2]+atom_j->r))+ m_probe_radius - m_min_z)/m_z_step));

	//sanity check
	int min_x_id = (tmp_min_x>0 ? tmp_min_x:0);
	int min_y_id = (tmp_min_y>0 ? tmp_min_y:0);
	int min_z_id = (tmp_min_z>0 ? tmp_min_z:0);

	int max_x_id = (tmp_max_x<m_x_num-1 ? tmp_max_x : m_x_num-1);
	int max_y_id = (tmp_max_y<m_y_num-1 ? tmp_max_y : m_y_num-1);
	int max_z_id = (tmp_max_z<m_z_num-1 ? tmp_max_z : m_z_num-1);

	//first check the neighboring grid points and also reset the tmp mark for grid edges
	for(int x_i=min_x_id; x_i<=max_x_id; x_i++)
		for(int y_i=min_y_id; y_i<=max_y_id; y_i++)
			for(int z_i=min_z_id; z_i<=max_z_id; z_i++){
				grid_pointIter current_point = m_grid_points[grid_index(x_i,y_i,z_i)];

				if(current_point->is_regular)
					continue;

				if(!current_torus->is_blocked && (current_torus->is_free || current_torus->saddle_face_list.size()>0))
					current_point->possible_torus_interior.push_back(current_torus);

				//if the point is inside the atoms or inside the concave spheres, leave the tag unchanged
				if(current_point->is_processed)
					continue;

				if(is_inside_accessible_torus(current_torus,current_point)){
					//current_point->is_inside_tori = true;
					current_point->is_processed = true;
					continue;
				}
				
				if(is_inside_Visibility_Sphere(current_torus,current_point->grid_pos)){
					current_point->is_inside_VS_or_dualTet = true;
				}
				
			}
}

//if return true, then this grid point is outside the molecular surface
bool is_inside_accessible_torus_old(TorusIter current_torus,grid_pointIter current_point){

	//if the torus does not exist
	if(current_torus->is_blocked)
		return false;

	double R = current_torus->torus_r;
	double r = m_probe_radius;
	CVector3d para_pos = (current_torus->torus_to_para).MultMatVec(current_point->grid_pos - current_torus->torus_center);
	double x = para_pos[0];
	double y = para_pos[1];
	double z = para_pos[2];

	//not inside the free torus
	if(is_f_InsideTorus_positive(x,y,z,R,r))
		return false;

	//for the self-intersecting part
	if(R<r && is_located_in_selfintersecting_part(x,y,z,R,r))
		return true;

	//current torus is the free torus
	if(current_torus->is_free)
		return true;

	if(is_covered_by_saddles(current_point->grid_pos,current_torus))
		return true;
	else
		return false;
}


//if return true, then this grid point is outside the molecular surface
bool is_inside_accessible_torus(TorusIter current_torus,grid_pointIter current_point){

	//if the torus does not exist
	if(current_torus->is_blocked)
		return false;

	double R = current_torus->torus_r;
	double r = m_probe_radius;
	CVector3d para_pos = (current_torus->torus_to_para).MultMatVec(current_point->grid_pos - current_torus->torus_center);
	double x = para_pos[0];
	double y = para_pos[1];
	double z = para_pos[2];

	if(is_f_quartic_equation_positive(x,y,z,R,r)){
		if(is_located_within_lemon_part(x,y,z,R,r))
			return true;
		else
			return false;
	}
	
	if(current_torus->is_free || is_covered_by_saddles(current_point->grid_pos,current_torus))
		return true;
	else
		return false;
}

bool is_inside_accessible_torus_with_output(TorusIter current_torus,grid_pointIter current_point){

	//if the torus does not exist
	if(current_torus->is_blocked)
		return false;

	double R = current_torus->torus_r;
	double r = m_probe_radius;
	CVector3d para_pos = (current_torus->torus_to_para).MultMatVec(current_point->grid_pos - current_torus->torus_center);
	double x = para_pos[0];
	double y = para_pos[1];
	double z = para_pos[2];

	if(is_f_quartic_equation_positive(x,y,z,R,r)){
		if(is_located_within_lemon_part(x,y,z,R,r))
			return true;
		else
			return false;
	}

	if(current_torus->is_free || is_covered_by_saddles(current_point->grid_pos,current_torus))
		return true;
	else
		return false;
}

void check_regular(){
	for(int i=0;i<m_grid_points.size();i++){
		grid_pointIter current_point = m_grid_points[i];

		if(current_point->is_outside_aug_atom)
			current_point->is_processed = true;
				
	}
	
	for(int i=0;i<m_grid_points.size();i++){
		grid_pointIter current_point = m_grid_points[i];
		if(!current_point->is_processed)
			continue;

		int x = current_point->grid_x;
		int y = current_point->grid_y;
		int z = current_point->grid_z;

		if(x==0 || y==0 || z==0)
			continue;
		if(x==m_x_num-1 || y==m_y_num-1 || z==m_z_num-1)
			continue;

		if(current_point->is_inside){
			grid_pointIter p1 = m_grid_points[grid_index(x-1,y,z)];
			if(!p1->is_processed || !p1->is_inside)
				continue;

			grid_pointIter p2 = m_grid_points[grid_index(x+1,y,z)];
			if(!p2->is_processed || !p2->is_inside)
				continue;

			grid_pointIter p3 = m_grid_points[grid_index(x,y-1,z)];
			if(!p3->is_processed || !p3->is_inside)
				continue;

			grid_pointIter p4 = m_grid_points[grid_index(x,y+1,z)];
			if(!p4->is_processed || !p4->is_inside)
				continue;

			grid_pointIter p5 = m_grid_points[grid_index(x,y,z-1)];
			if(!p5->is_processed || !p5->is_inside)
				continue;

			grid_pointIter p6 = m_grid_points[grid_index(x,y,z+1)];
			if(!p6->is_processed || !p6->is_inside)
				continue;

			current_point->is_regular = true;
			current_point->possible_atom_interior.clear();
			current_point->possible_concavesphere_interior.clear();

		}

		else{
			grid_pointIter p1 = m_grid_points[grid_index(x-1,y,z)];
			if(!p1->is_processed || p1->is_inside)
				continue;

			grid_pointIter p2 = m_grid_points[grid_index(x+1,y,z)];
			if(!p2->is_processed || p2->is_inside)
				continue;

			grid_pointIter p3 = m_grid_points[grid_index(x,y-1,z)];
			if(!p3->is_processed || p3->is_inside)
				continue;

			grid_pointIter p4 = m_grid_points[grid_index(x,y+1,z)];
			if(!p4->is_processed || p4->is_inside)
				continue;

			grid_pointIter p5 = m_grid_points[grid_index(x,y,z-1)];
			if(!p5->is_processed || p5->is_inside)
				continue;

			grid_pointIter p6 = m_grid_points[grid_index(x,y,z+1)];
			if(!p6->is_processed || p6->is_inside)
				continue;

			current_point->is_regular = true;
			current_point->possible_atom_interior.clear();
			current_point->possible_concavesphere_interior.clear();
		}
		

	}
	
}


void check_the_grid_points(){
	for(int i=0;i<m_grid_points.size();i++){
		grid_pointIter current_point = m_grid_points[i];
		if(current_point->is_processed)
			continue;
		if(current_point->is_inside_VS_or_dualTet){
			current_point->is_processed = true;
			current_point->is_inside = true;
		}
		
	}
	
}

//detect the other unprocessed grid points by counting the number of intersection points
void process_others_ray_counting(){
	for(int i=0;i<m_grid_points.size();i++){
		grid_pointIter current_point = m_grid_points[i];
		if(current_point->is_processed)
			continue;
		process_ray_racing_each_point(current_point);
	}
}

void process_ray_racing_each_point(grid_pointIter current_point){
	grid_pointIter processed_ngb;
	int found_k,mode;
	mode = -1;
	if(search_the_processed_neighbor(current_point,processed_ngb,found_k,mode)){
		int direction_x = 1;
		if(found_k > 0)
			direction_x = -1;
		int process_x;
		grid_pointIter current_process=NULL;

		for(int k=0;k<abs(found_k);k++){
			if(mode==0){
				process_x = processed_ngb->grid_x + direction_x;
				current_process = m_grid_points[grid_index(process_x,processed_ngb->grid_y,processed_ngb->grid_z)];
			}
			if(mode==1){
				process_x = processed_ngb->grid_y + direction_x;
				current_process = m_grid_points[grid_index(processed_ngb->grid_x,process_x,processed_ngb->grid_z)];
			}

			if(mode==2){
				process_x = processed_ngb->grid_z + direction_x;
				current_process = m_grid_points[grid_index(processed_ngb->grid_x,processed_ngb->grid_y,process_x)];
			}


			int num_intersection = compute_surface_intersection_num(current_process,processed_ngb);
			if((num_intersection%2 && !processed_ngb->is_inside) || (num_intersection%2==0 && processed_ngb->is_inside))
				current_process->is_inside = true;
			current_process->is_processed = true;
			processed_ngb = current_process;
		}
	}
	else
		current_point->is_processed = true; //default: outside the surface

}

void process_each_augmented_atom(AtomIter current_atom){
	//decide which gird points are inside this atom sphere
	int tmp_min_x = int(floor((current_atom->center[0] - current_atom->r - m_probe_radius - m_min_x)/m_x_step));
	int tmp_min_y = int(floor((current_atom->center[1] - current_atom->r - m_probe_radius - m_min_y)/m_y_step));
	int tmp_min_z = int(floor((current_atom->center[2] - current_atom->r - m_probe_radius - m_min_z)/m_z_step));

	int tmp_max_x = int(ceil((current_atom->center[0] + current_atom->r + m_probe_radius - m_min_x)/m_x_step));
	int tmp_max_y = int(ceil((current_atom->center[1] + current_atom->r + m_probe_radius - m_min_y)/m_y_step));
	int tmp_max_z = int(ceil((current_atom->center[2] + current_atom->r + m_probe_radius - m_min_z)/m_z_step));

	//sanity check
	int min_x_id = (tmp_min_x>0 ? tmp_min_x:0);
	int min_y_id = (tmp_min_y>0 ? tmp_min_y:0);
	int min_z_id = (tmp_min_z>0 ? tmp_min_z:0);

	int max_x_id = (tmp_max_x<m_x_num-1 ? tmp_max_x : m_x_num-1);
	int max_y_id = (tmp_max_y<m_y_num-1 ? tmp_max_y : m_y_num-1);
	int max_z_id = (tmp_max_z<m_z_num-1 ? tmp_max_z : m_z_num-1);

	for(int x_i=min_x_id; x_i<=max_x_id; x_i++)
		for(int y_i=min_y_id; y_i<=max_y_id; y_i++)
			for(int z_i=min_z_id; z_i<=max_z_id; z_i++){
				grid_pointIter current_point = m_grid_points[grid_index(x_i,y_i,z_i)];
				if(current_point->is_processed)
					continue;
				grid_pointIter processed_ngb;
				int found_k,mode;
				mode = -1;
				if(search_the_processed_neighbor(current_point,processed_ngb,found_k,mode)){
					int direction_x = 1;
					if(found_k > 0)
						direction_x = -1;
					int process_x;
					grid_pointIter current_process;

					for(int k=0;k<abs(found_k);k++){
						if(mode==0){
							process_x = processed_ngb->grid_x + direction_x;
							current_process = m_grid_points[grid_index(process_x,processed_ngb->grid_y,processed_ngb->grid_z)];
						}
						if(mode==1){
							process_x = processed_ngb->grid_y + direction_x;
							current_process = m_grid_points[grid_index(processed_ngb->grid_x,process_x,processed_ngb->grid_z)];
						}

						if(mode==2){
							process_x = processed_ngb->grid_z + direction_x;
							current_process = m_grid_points[grid_index(processed_ngb->grid_x,processed_ngb->grid_y,process_x)];
						}


						int num_intersection = compute_surface_intersection_num(current_process,processed_ngb);
						if((num_intersection%2 && !processed_ngb->is_inside) || (num_intersection%2==0 && processed_ngb->is_inside))
							current_process->is_inside = true;
						current_process->is_processed = true;
						processed_ngb = current_process;
					}
				}
				else
					current_point->is_processed = true; //default: outside the surface


			}
}


int compute_surface_intersection_num(grid_pointIter current_gridpoint,grid_pointIter processed_ngb){
	grid_pointIter inside_point, outside_point;
	if(processed_ngb->is_inside){
		inside_point = processed_ngb;
		outside_point = current_gridpoint;
	}
	else{
		inside_point = current_gridpoint;
		outside_point = processed_ngb;
	}

	return(compute_atom_intersection_num(inside_point,outside_point) +
		compute_torus_intersection_num(inside_point,outside_point) +
		compute_concave_intersection_num(inside_point,outside_point));

}
int compute_intersection_single_atom(grid_pointIter point_inside, grid_pointIter point_outside,AtomIter atom,std::vector<double> &t){
	/***
	compute the intersection point given the equation from wiki page:
	http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
	******/
	double r = atom->r;
	CVector3d o = point_inside->grid_pos;
	CVector3d c = atom->center;
	CVector3d l = point_outside->grid_pos - point_inside->grid_pos;
	l.Normalize();

	int num_intersection = 0;

	//for the special case: when the grid edge is tangential to the atom
	if((o-c).LengthSquared()- r*r ==0 && (point_outside->grid_pos-c).LengthSquared() - r*r >0){
		double d1 = 0;
		if(validate_atom_intersection(atom,o)){
			num_intersection++;
			t.push_back(d1);
			
		}

		return num_intersection;
	}
	

	double tmp1 = -l.Dot(o-c);
	double check_value = tmp1*tmp1 - (o-c).LengthSquared() + r*r;

	if(check_value < 0)
		return 0;

	double tmp2 = sqrt(check_value);

	double d1;
	
	for(int j=0;j<2;j++){
		if(j==0)
			d1 = tmp1 + tmp2; 
		else{
			/*if(fabs(tmp2)<ZERO_TOL)
				continue;*/
			d1 = tmp1 - tmp2;
		}


		if((o-c).LengthSquared()- r*r<=0 && (point_outside->grid_pos-c).LengthSquared() - r*r >0)
			if(fabs(d1)<ZERO_TOL)
				d1 = ZERO_TOL;

		// d1 has to be positive and smaller than the edge length
		if(d1 < 0 || d1>(point_inside->grid_pos - point_outside->grid_pos).Length())
			continue;

		//check if it is on the accessible patch
		CVector3d tmp_pos = o + d1 * l;
		
		if(validate_atom_intersection(atom,tmp_pos)){
			num_intersection++;
			t.push_back(d1);
		}
	}

	return num_intersection;


}
int compute_atom_intersection_num(grid_pointIter point_inside, grid_pointIter point_outside){
	int num_intersection = 0;
	CVector3d intersect_pos, intersect_normal;
			
	//check the intersection points with all the POSSIBLE atoms which include current grid point
	for(int i=0;i<point_inside->possible_atom_interior.size();i++){
		AtomIter atom = point_inside->possible_atom_interior[i];
		std::vector<double> t;

		num_intersection += compute_intersection_single_atom(point_inside,point_outside,atom,t);
	}
	
	return num_intersection;	
}



int compute_torus_intersection_num(grid_pointIter point_inside,grid_pointIter point_outside){
	int num_intersection = 0;

	for(int i=0;i<point_inside->possible_torus_interior.size();i++){
		TorusIter current_torus = point_inside->possible_torus_interior[i];

		if(current_torus->is_blocked)
			continue;

		//then there is an intersection with the torus
		std::vector<double> t;
		std::vector<CVector3d> inter_normal;

		//compute the parametrized pos in a standard torus
		CVector3d pos_interior = point_inside->grid_pos;
		CVector3d pos_outside = point_outside->grid_pos;
		CVector3d para_interior = (current_torus->torus_to_para).MultMatVec(pos_interior - current_torus->torus_center);
		CVector3d para_outside = (current_torus->torus_to_para).MultMatVec(pos_outside - current_torus->torus_center);

		num_intersection += compute_line_torus_intersection(current_torus,pos_interior,pos_outside,para_interior,para_outside,t,inter_normal);

	}

	return num_intersection;

}

//----------------------------------------------------------------------------
bool solveQuartic(double a, double b, double c, double d, double e, double *root)
{

	//switch to the Jenkins-Traub algorithm
	solve_quartic_equation(a,b,c,d,e,root);
	return true;
}

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
									CVector3d pos_outside,std::vector<double> &t,std::vector<CVector3d> &inter_normal){
	CVector3d p = pos_interior;
	CVector3d td = pos_outside - pos_interior;  //then t should be (0,1)
	double R = current_torus->torus_r;
	double r = m_probe_radius;

	double alpha = td.Dot(td);
	double beta = 2*p.Dot(td);
	double y = p.Dot(p) - R*R - r*r;

	//the coefficients for quartic equations
	double a = alpha*alpha;
	double b = 2*alpha*beta;
	double c = beta*beta + 2*alpha*y + 4*R*R*td[2]*td[2];
	double d = 2*beta*y + 8*R*R*p[2]*td[2];
	double e = y*y + 4*R*R*(p[2]*p[2] - r*r);

	int num_intersection = 0;

	double root[4];
	for(int i=0;i<4;i++)
		root[i] = -1;

	//if(check_value){
	//	std::cout<<std::setprecision(10)<<"a: "<<a<<" b: "<<b<<" c: "<<c<<" d: "<<d<<" e: "<<e<<std::endl;
	//	solveQuartic(a,b,c,d,e,root,true);
	//}

	solveQuartic(a,b,c,d,e,root);
	
	for(int i=0;i<4;i++){
		if(root[i] >1.0 || root[i]<0)
			continue;

		//compute the normal 
		CVector3d point = p + root[i]*td;
		double x = point[0]; double y=point[1]; double z = point[2];


		//the intersection point should be within the visibility sphere
		CVector3d pos = grid_pos_interior + root[i]*(grid_pos_outside-grid_pos_interior);
		

		if(!is_located_on_lemon_part(x,y,z,R,r) && is_inside_Visibility_Sphere(current_torus,pos)){
			//	std::cout<<"f value "<<f_value(a,b,c,d,e,t)<<std::endl;
			
			t.push_back(root[i]);
			double normal_x = 4*x*(x*x + y*y + z*z - r*r - R*R);
			double normal_y = 4*y*(x*x + y*y + z*z - r*r - R*R);
			double normal_z = 4*z*(x*x + y*y + z*z - r*r - R*R) + 8*R*R*z;
			CVector3d tmp_normal = CVector3d(normal_x,normal_y,normal_z);
			tmp_normal.Normalize();
			inter_normal.push_back(tmp_normal);

			
			num_intersection++;
		}
	}
	return num_intersection;
}

int compute_concave_intersection_num(grid_pointIter point_inside,grid_pointIter point_outside){

	int num_intersection = 0;

	std::vector<CVector3d> intersection_pos;
	std::vector<double> intersect_d1;
	std::vector<CVector3d> intersect_probe_center;
	std::vector<Concave_SphereIter> intersect_concave;

	for(int i=0;i<point_inside->possible_concavesphere_interior.size();i++){
		Concave_SphereIter current_concave = point_inside->possible_concavesphere_interior[i];
		CVector3d o, pos_outside;
		o = point_inside->grid_pos;
		pos_outside = point_outside->grid_pos;


		/****************************************************************
		compute the intersection point given the equation from wiki page:
		http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
		*******************************************************************/
		double r = m_probe_radius;
		CVector3d c = current_concave->probe_center;
		CVector3d l = pos_outside - o ;
		l.Normalize();

		double tmp1 = -l.Dot(o-c);
		double check_value = tmp1*tmp1 - (o-c).LengthSquared() + r*r;


		if(check_value<0)
			continue;

		double tmp2 = sqrt(check_value);
		double d1;
		for(int j=0;j<2;j++){
			if(j==0)
				d1 = tmp1 + tmp2;
			else
				d1 = tmp1 - tmp2;



				if(d1 > (pos_outside-o).Length() || d1<0)
					continue;

				CVector3d intersect_point = o + d1 * l;

				if(is_inside_spherical_triangle(current_concave,intersect_point)){
					bool is_exist = false;
					//in case we have the probe sphere touching more than three atoms at the same time
					for(int k=0;k<intersect_probe_center.size();k++){
						if((c-intersect_probe_center[k]).Length() < ZERO_TOL && fabs(d1 - intersect_d1[k])<ZERO_TOL){
							is_exist = true;
							break;
						}
					}
					if(!is_exist){
						intersect_probe_center.push_back(c);
						intersect_d1.push_back(d1);
						intersection_pos.push_back(intersect_point);
						intersect_concave.push_back(current_concave);

					}
				}
		}
	}

	//outside the concave sphere
	for(int i=0;i<intersection_pos.size();i++){
		bool is_inside_concave = false;
		for(int j=0;j<point_inside->possible_concavesphere_interior.size();j++){
			Concave_SphereIter current_concave = point_inside->possible_concavesphere_interior[j];
			bool is_duplex = false;
			if(current_concave==intersect_concave[i])
				is_duplex = true;
			if(current_concave->duplex_concave_faces.size()>0){
				for(int k=0;k<current_concave->duplex_concave_faces.size();k++){
					if(current_concave->duplex_concave_faces[k]==intersect_concave[i]){
							is_duplex = true;
							break;
					}
				}
			}

			if(is_duplex)
				continue;
		
			if(is_inside_concave_sphere(intersection_pos[i],current_concave->probe_center)){
				is_inside_concave = true;
				break;
			}
		}
		if(!is_inside_concave)
			num_intersection++;
	}


		
	return num_intersection;
}

/**********************labeling the grid points ends**********************************/
/**********************labeling the grid points ends**********************************/
/**********************labeling the grid points ends**********************************/







/********************compute intersection starts**************************************/
/********************compute intersection starts**************************************/
/********************compute intersection starts**************************************/


 int process_grid_intersection_new(grid_edgeIter current_edge){
	if(current_edge->point1->is_inside && current_edge->point2->is_inside)
		return 0;
	if(!current_edge->point1->is_inside && !current_edge->point2->is_inside)
		return 0;

	grid_pointIter point_inside, point_outside;
	if(current_edge->point1->is_inside){
		point_inside = current_edge->point1;
		point_outside = current_edge->point2;
	}
	else{
		point_inside = current_edge->point2;
		point_outside = current_edge->point1;
	}
	std::vector<CVector3d> intersect_points; //store the possible intersection points on this grid edge
	std::vector<CVector3d> intersect_normal; //store the normal direction of the possible intersection points
	std::vector<int> intersect_type; //store the type of intersection (0: atom; 1: torus; 2: concave);
	std::vector<cell_unit> intersect_cell;

	process_grid_atom_intersection_new(point_inside,point_outside,intersect_points,intersect_normal,intersect_type,intersect_cell);
	for(int i=0;i<point_inside->possible_torus_interior.size();i++)
		process_grid_torus_intersection_new(point_inside,point_outside,point_inside->possible_torus_interior[i],intersect_points,intersect_normal,intersect_type,intersect_cell);
	process_grid_edge_concave_new(point_inside,point_outside,intersect_points,intersect_normal,intersect_type,intersect_cell);

	return (process_each_grid_edge(point_inside,point_outside,intersect_points,intersect_normal,intersect_type,intersect_cell));
}

void process_grid_atom_intersection_new(grid_pointIter point_inside, grid_pointIter point_outside,std::vector<CVector3d> &intersect_points,
										std::vector<CVector3d> &intersect_normal,std::vector<int> &intersect_type,std::vector<cell_unit> &intersect_cell){

	double track_l= -1;
	CVector3d o = point_inside->grid_pos;
	CVector3d l = point_outside->grid_pos - point_inside->grid_pos;
	l.Normalize();

	bool is_found = false;
	CVector3d intersect_pos, intersect_normal_tmp;

	AtomIter intersect_atom=NULL;		

	//check the intersection points with all the POSSIBLE atoms which include current grid point
	for(int i=0;i<point_inside->possible_atom_interior.size();i++){
		AtomIter atom = point_inside->possible_atom_interior[i];
		std::vector<double> t;
		if(compute_intersection_single_atom(point_inside,point_outside,atom,t)){
			double d1 = t[0];
			if(t.size()>1)
				d1 = max(t[0],t[1]);

			if(track_l < d1){
				track_l = d1;
				is_found = true;

				intersect_atom = atom;
				//we know d1, d2 should be one positive and one negative, and the positive one (d1)is the intersection we are looking for
				intersect_pos = o + d1 * l;
				intersect_normal_tmp = intersect_pos - atom->center;
				intersect_normal_tmp.Normalize();
			}
		}
	}

	if(is_found){
		intersect_points.push_back(intersect_pos);
		intersect_normal.push_back(intersect_normal_tmp);
		intersect_type.push_back(0);

		cell_unit current_cell;
		current_cell.current_atom = intersect_atom;
		intersect_cell.push_back(current_cell);
	}
}

void process_grid_torus_intersection_new(grid_pointIter point_inside,grid_pointIter point_outside,TorusIter current_torus, std::vector<CVector3d> &intersect_points,
										 std::vector<CVector3d> &intersect_normal,std::vector<int> &intersect_type,std::vector<cell_unit> &intersect_cell){
	
	if(current_torus->is_blocked)
		return;

	//then there is an intersection with the torus
	std::vector<double> t;
	std::vector<CVector3d> inter_normal;
	
	//compute the parametrized pos in a standard torus
	CVector3d pos_interior = point_inside->grid_pos;
	CVector3d pos_outside = point_outside->grid_pos;
	CVector3d para_interior = (current_torus->torus_to_para).MultMatVec(pos_interior - current_torus->torus_center);
	CVector3d para_outside = (current_torus->torus_to_para).MultMatVec(pos_outside - current_torus->torus_center);



	if(compute_line_torus_intersection(current_torus,pos_interior,pos_outside,para_interior,para_outside,t,inter_normal)){
			
		for(int i=0;i<t.size();i++){
			CVector3d intersection_point = pos_interior + t[i]*(pos_outside - pos_interior);
			intersect_points.push_back(intersection_point);
			CVector3d tmp_normal = current_torus->para_to_torus.MultMatDir(inter_normal[i]);
			tmp_normal.Normalize();
			intersect_normal.push_back(-tmp_normal);
			intersect_type.push_back(1);

			cell_unit current_cell;
			current_cell.current_torus = current_torus;
			intersect_cell.push_back(current_cell);
		}
	}	
}

void process_grid_edge_concave_new(grid_pointIter point_inside,grid_pointIter point_outside,std::vector<CVector3d> &intersect_points,
								   std::vector<CVector3d> &intersect_normalvec,std::vector<int> &intersect_type,std::vector<cell_unit> &intersect_cell){

	//then there is an intersection with the torus
	CVector3d o, pos_outside;
	o = point_inside->grid_pos;
	pos_outside = point_outside->grid_pos;


	/****************************************************************
	compute the intersection point given the equation from wiki page:
	http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
	*******************************************************************/
	double r = m_probe_radius;
	std::vector<CVector3d> possible_intersection;
	std::vector<CVector3d> possible_normal;
	std::vector<cell_unit> possible_cell;

	for(int i=0;i<point_inside->possible_concavesphere_interior.size();i++){
		Concave_SphereIter current_concave = point_inside->possible_concavesphere_interior[i];
		CVector3d c = current_concave->probe_center;
		CVector3d l = pos_outside - o ;
		l.Normalize();

		double tmp1 = -l.Dot(o-c);
		double check_value = tmp1*tmp1 - (o-c).LengthSquared() + r*r;
		if(check_value<0)
			continue;
		double tmp2 = sqrt(check_value);
		double d1;
		for(int j=0;j<2;j++){
			if(j==0)
				d1 = tmp1 + tmp2;
			else
				d1 = tmp1 - tmp2;

			if(d1 > (pos_outside-o).Length() || d1<0)
				continue;

			

			CVector3d intersect_point = o + d1 * l;

			if(is_inside_spherical_triangle(current_concave,intersect_point)){

						CVector3d intersect_normal = current_concave->probe_center - intersect_point;
						intersect_normal.Normalize();

						possible_intersection.push_back(intersect_point);
						possible_normal.push_back(intersect_normal);

						cell_unit current_cell;
						current_cell.current_concave = current_concave;
						possible_cell.push_back(current_cell);
			}
		}
	}

	for(int i=0;i<possible_intersection.size();i++){
		bool is_inside_concave = false;
		
		for(int j=0;j<point_inside->possible_concavesphere_interior.size();j++){
			Concave_SphereIter current_concave = point_inside->possible_concavesphere_interior[j];
			bool is_duplex = false;

		
			if(current_concave==possible_cell[i].current_concave)
				is_duplex = true;

			
			if(current_concave->duplex_concave_faces.size()>0){
				
				for(int k=0;k<current_concave->duplex_concave_faces.size();k++){
					if(current_concave->duplex_concave_faces[k]==possible_cell[i].current_concave){
							is_duplex = true;
							break;
					}
				}
			}

			if(is_duplex)
				continue;

			if(is_inside_concave_sphere(possible_intersection[i],current_concave->probe_center)){
				
				is_inside_concave = true;
				break;
			}
		}

		if(!is_inside_concave){
			intersect_points.push_back(possible_intersection[i]);
			intersect_normalvec.push_back(possible_normal[i]);
			intersect_type.push_back(2);
			intersect_cell.push_back(possible_cell[i]);
		}
	}
		
}

/********************compute intersection ends**************************************/
/********************compute intersection ends**************************************/
/********************compute intersection ends**************************************/


void clean_memory(){
    //clean the grids
  //cout<<"@@@@@@@@@@@@"<<endl;
  for(int i=0;i<m_grid_edges.size();i++)
    delete m_grid_edges[i];
  //cout<<"@@@@@@@@@@@@"<<endl;
  for(int i=0;i<m_grid_points.size();i++)
    delete m_grid_points[i];
  //cout<<"@@@@@@@@@@@@"<<endl;
  m_grid_edges.clear();
  m_grid_points.clear();

  m_intersection_points.clear(); //store the output: location of intersection points
  m_intersection_normals.clear(); //store the normal direction of the intersection points
  m_intersection_interior.clear();
  m_intersection_outside.clear();
  m_intersection_cell.clear();
  m_intersection_types.clear(); // 0: atom; 1:torus; 2:concave
  
	//clean the surface
	for(int i=0;i<m_atoms_vec.size();i++)
		delete m_atoms_vec[i];
	for(int i=0;i<m_torus_vec.size();i++)
		delete m_torus_vec[i];
	for(int i=0;i<m_concave_spheres.size();i++)
		delete m_concave_spheres[i];
	for(int i=0;i<m_saddle_faces.size();i++)
		delete m_saddle_faces[i];
	for(int i=0;i<m_concave_edges.size();i++)
		delete m_concave_edges[i];

}


/**********************
	here we check whether it is inside the torus or outside 
	 the parametrized expression for torus is:
	 (R - sqrt(x*x + y*y) ) ^2 + z^2 = r^2

	 To be in the 'interior' of this torus patch, it should require
	 (1) [R - sqrt(x*x+y*y)]^2 + z^2 - r^2 >0
	 (2) within the visibility sphere (Fig 7. of paper "Interactive Visualization of Molecular Surface Dynamics")
************************************************************/
bool is_inside_torus(TorusIter current_torus,grid_pointIter current_point){
//check if it is within the visibility sphere
	double check2 = (current_point->grid_pos - current_torus->vs_sphere_center).Length() - current_torus->vs_sphere_r;
	
	CVector3d para_pos = (current_torus->torus_to_para).MultMatVec(current_point->grid_pos - current_torus->torus_center);
	

	if(check2>0)
		return false;

	double x = para_pos[0];
	double y = para_pos[1];
	double z = para_pos[2];
	double R = current_torus->torus_r;
	double r = m_probe_radius;

	//check if it is outside the torus
	double tmp = R - sqrt(x*x+y*y);
	double check1 = tmp*tmp + z*z - r*r;

	if(check1<0)
		return false;

	return true;



}

/**********************
	here we check whether it is inside the torus or outside 
	 the parametrized expression for torus is:
	 (R - sqrt(x*x + y*y) ) ^2 + z^2 = r^2

	 To be in the 'interior' of this torus patch, it should require
	 (1) [R - sqrt(x*x+y*y)]^2 + z^2 - r^2 >0
	 (2) within the visibility sphere (Fig 7. of paper "Interactive Visualization of Molecular Surface Dynamics")
************************************************************/
bool is_inside_torus_new(TorusIter current_torus,grid_pointIter current_point){
//check if it is within the visibility sphere
	double check2 = (current_point->grid_pos - current_torus->vs_sphere_center).Length() - current_torus->vs_sphere_r;
	if(check2>0)
		return false;

	CVector3d para_pos = (current_torus->torus_to_para).MultMatVec(current_point->grid_pos - current_torus->torus_center);
	double x = para_pos[0];
	double y = para_pos[1];
	double z = para_pos[2];
	double R = current_torus->torus_r;
	double r = m_probe_radius;

	//check if it is outside the torus
	double tmp = R - sqrt(x*x+y*y);
	double check1 = tmp*tmp + z*z - r*r;

	if(check1<0)
		return false;

	return true;

}







bool is_outside_concave_sphere(TorusIter current_torus,CVector3d current_point){
	for(int i=0;i<current_torus->concave_edge_list.size();i++){
		Concave_SphereIter current_concave = current_torus->concave_edge_list[i]->current_concave;
		if((current_point-current_concave->probe_center).LengthSquared() < m_probe_radius*m_probe_radius)
			return false;
		
	}

	return true;
}


//find the closest neighboring grid point which is already processed
bool search_the_processed_neighbor(grid_pointIter current_point, grid_pointIter &processed_neighbor,int &found_k,int &mode){
	std::vector<int> store_k;
	std::vector<int> store_mode;
	std::vector<grid_pointIter> store_point;

	for(int k=1;k<NEIGHBOR_LENGTH;k++){
		int neighbor_x = current_point->grid_x + k;
		if(neighbor_x<0 || neighbor_x>m_x_num-1)
			break;
		grid_pointIter tmp_point = m_grid_points[grid_index(neighbor_x,current_point->grid_y,current_point->grid_z)];
		if(tmp_point->is_processed){
			processed_neighbor = tmp_point;
			found_k = k;
			mode = 0;

			store_k.push_back(found_k);
			store_mode.push_back(mode);
			store_point.push_back(processed_neighbor);
			break;
			
			//return true;
		}
	}
	for(int k=-1;k>-NEIGHBOR_LENGTH;k--){
		int neighbor_x = current_point->grid_x + k;
		if(neighbor_x<0 || neighbor_x>m_x_num-1)
			break;
		grid_pointIter tmp_point = m_grid_points[grid_index(neighbor_x,current_point->grid_y,current_point->grid_z)];
		if(tmp_point->is_processed){
			processed_neighbor = tmp_point;
			found_k = k;
			mode = 0;
			
			store_k.push_back(found_k);
			store_mode.push_back(mode);
			store_point.push_back(processed_neighbor);
			break;
			//return true;
		}
	}
	
	for(int k=1;k<NEIGHBOR_LENGTH;k++){
		int neighbor_y = current_point->grid_y + k;
		if(neighbor_y<0 || neighbor_y>m_y_num-1)
			break;
		grid_pointIter tmp_point = m_grid_points[grid_index(current_point->grid_x,neighbor_y,current_point->grid_z)];
		if(tmp_point->is_processed){
			processed_neighbor = tmp_point;
			found_k = k;
			mode = 1;

			store_k.push_back(found_k);
			store_mode.push_back(mode);
			store_point.push_back(processed_neighbor);
			break;
			//return true;
		}
	}

	for(int k=-1;k>-NEIGHBOR_LENGTH;k--){
		int neighbor_y = current_point->grid_y + k;
		if(neighbor_y<0 || neighbor_y>m_y_num-1)
			break;
		grid_pointIter tmp_point = m_grid_points[grid_index(current_point->grid_x,neighbor_y,current_point->grid_z)];
		if(tmp_point->is_processed){ 
			processed_neighbor = tmp_point;
			found_k = k;
			mode = 1;

			store_k.push_back(found_k);
			store_mode.push_back(mode);
			store_point.push_back(processed_neighbor);
			break;
			//return true;
		}
	}

	for(int k=0;k<NEIGHBOR_LENGTH;k++){
		int neighbor_z = current_point->grid_z + k;
		if(neighbor_z<0 || neighbor_z>m_z_num-1)
			break;
		grid_pointIter tmp_point = m_grid_points[grid_index(current_point->grid_x,current_point->grid_y,neighbor_z)];
		if(tmp_point->is_processed){
			processed_neighbor = tmp_point;
			found_k = k;
			mode = 2;
			
			store_k.push_back(found_k);
			store_mode.push_back(mode);
			store_point.push_back(processed_neighbor);
			break;

			//return true;
		}
	}
	for(int k=-1;k>-NEIGHBOR_LENGTH;k--){
		int neighbor_z = current_point->grid_z + k;
		if(neighbor_z<0 || neighbor_z>m_z_num-1)
			break;
		grid_pointIter tmp_point = m_grid_points[grid_index(current_point->grid_x,current_point->grid_y,neighbor_z)];
		if(tmp_point->is_processed){
			processed_neighbor = tmp_point;
			found_k = k;
			mode = 2;
			
			store_k.push_back(found_k);
			store_mode.push_back(mode);
			store_point.push_back(processed_neighbor);
			break;

			//return true;
		}
	}




	if(store_k.size()==0){
		std::cout<<std::endl<<"may need to extend the neighbor length!!"<<std::endl;
		return false;
	}
	else{
		found_k = store_k[0];
		mode = store_mode[0];
		processed_neighbor = store_point[0];

		for(int l=1;l<store_k.size();l++){
			if(abs(found_k) > abs(store_k[l])){
				found_k = store_k[l];
				mode = store_mode[l];
				processed_neighbor = store_point[l];
			}
		}

		return true;
	}
	

}



bool solveQuadraticOther(double a, double b, double c, double &root)
{
	if(a == 0.0 || fabs(a/b) < 1.0e-6)
	{
		if(fabs(b) < 1.0e-6) 
			return false;
		else
		{
			root = -c/b;
			return true;
		}
	}

	double discriminant = b * b - 4.0 * a * c;
	if(discriminant >= 0.0)
	{
		discriminant = sqrt(discriminant);
		root = (b - discriminant) * -0.5 / a;
		return true;
	}

	return false;
}

//----------------------------------------------------------------------------
bool solveQuadratic(double a, double b, double c, double &root)
{
	if(a == 0.0 || fabs(a/b) < 1.0e-6)
	{
		if(fabs(b) < 1.0e-6) 
			return false;
		else
		{
			root = -c/b;
			return true;
		}
	}

	double discriminant = b * b - 4.0 * a * c;
	if(discriminant >= 0.0)
	{
		discriminant = sqrt(discriminant);
		root = (b + discriminant) * -0.5 / a;
		return true;
	}

	return false;
}

bool solveCubic(double a, double b, double c, double d, double &root)
{
	if(a == 0.0 || fabs(a/b) < 1.0e-6)
		return solveQuadratic(b, c, d, root);
	double B = b/a, C = c/a, D = d/a;

	double Q = (B*B - C*3.0)/9.0, QQQ = Q*Q*Q;
	double R = (2.0*B*B*B - 9.0*B*C + 27.0*D)/54.0, RR = R*R;

	// 3 real roots
	if(RR<QQQ)
	{
		/* This sqrt and division is safe, since RR >= 0, so QQQ > RR,    */
		/* so QQQ > 0.  The acos is also safe, since RR/QQQ < 1, and      */
		/* thus R/sqrt(QQQ) < 1.                                     */
		double theta = acos(R/sqrt(QQQ));
		/* This sqrt is safe, since QQQ >= 0, and thus Q >= 0             */
		double r1, r2, r3;
		r1 = r2 = r3 = -2.0*sqrt(Q);
		r1 *= cos(theta/3.0);
		r2 *= cos((theta+2*PI)/3.0);
		r3 *= cos((theta-2*PI)/3.0);

		r1 -= B/3.0;
		r2 -= B/3.0;
		r3 -= B/3.0; 

	
		//initialization of root!
		root = -1e10;

		if(r1 >= 0.0) root = r1;
		if(r2 >= 0.0 && r2 > root) root = r2;
		if(r3 >= 0.0 && r3 > root) root = r3;

		return true;
	}
	// 1 real root
	else
	{
		double tmp = R+sqrt(RR-QQQ);
		//for numerical issue
		if(fabs(QQQ/RR)<1e-8 && QQQ>0){
			tmp = 2*R;
		}
		double A2;
		if(tmp<0)
			A2 = -pow(-tmp,1/3.0);
		else
			A2 = pow(tmp,1/3.0);

		if (A2!=0.0) {		
			root = - A2 - Q/A2; 
		}
		root -= B/3.0;
		return true;
	}
}







double f_value(double a ,double b, double c, double d, double e, double x){
	return (a*pow(x,4) + b*pow(x,3) + c*pow(x,2) + d*x + e);
}

double f_derivative_value(double a ,double b, double c, double d, double e, double x){
	return (4*a*pow(x,3) + 3*b*pow(x,2) + 2*c*x + d);
}

bool solve_newton_method(double a ,double b, double c, double d, double e,double x_0, double &t){
	if(fabs(f_value(a,b,c,d,e,x_0)) < ZERO_TOL){
		t = x_0;
		return true;
	}

	double x_1 = x_0 - f_value(a,b,c,d,e,x_0)/f_derivative_value(a,b,c,d,e,x_0);
	
	int ntimes = 0;

	while(fabs(x_1 - x_0) > ZERO_TOL  && ntimes<100){
		x_0 = x_1;

		double test1 = f_value(a,b,c,d,e,x_0);
		double test2 = f_derivative_value(a,b,c,d,e,x_0);

		double check_0 = f_derivative_value(a,b,c,d,e,0.0);
		double check_1 = f_derivative_value(a,b,c,d,e,1.0);

		x_1 =  x_0 - f_value(a,b,c,d,e,x_0)/f_derivative_value(a,b,c,d,e,x_0);

		ntimes++;
	}

	if(ntimes>= 100)
		return false;

	if(x_1 >= 0 && x_1<= 1.0){
		t = x_1;
		std::cout<<"f value:"<<f_value(a,b,c,d,e,t)<<std::endl;
		return true;
	}
	else
		return false;
}

bool is_inside_spherical_triangle(Concave_SphereIter current_concave,CVector3d test_point){
	CVector3d p1 = current_concave->center_i;
	CVector3d p2 = current_concave->center_j;
	CVector3d p3 = current_concave->center_k;
	CVector3d p4 = current_concave->probe_center;

	CMatrix44 D0 = CMatrix44(p1[0],p1[1],p1[2],1.0,p2[0],p2[1],p2[2],1.0,
		p3[0],p3[1],p3[2],1.0,p4[0],p4[1],p4[2],1.0);

	CMatrix44 D1 = CMatrix44(test_point[0],test_point[1],test_point[2],1.0,p2[0],p2[1],p2[2],1.0,
		p3[0],p3[1],p3[2],1.0,p4[0],p4[1],p4[2],1.0);

	CMatrix44 D2 = CMatrix44(p1[0],p1[1],p1[2],1.0,test_point[0],test_point[1],test_point[2],1.0,
		p3[0],p3[1],p3[2],1.0,p4[0],p4[1],p4[2],1.0);

	CMatrix44 D3 = CMatrix44(p1[0],p1[1],p1[2],1.0,p2[0],p2[1],p2[2],1.0,
		test_point[0],test_point[1],test_point[2],1.0,p4[0],p4[1],p4[2],1.0);



	double check_D0 = D0.Determinant();
	double check_D1 = D1.Determinant();
	double check_D2 = D2.Determinant();
	double check_D3 = D3.Determinant();

	
	if((check_D0>=0 && check_D1>=0 && check_D2>=0 && check_D3>=0) ||
		(check_D0<=-0 && check_D1<=-0 && check_D2<=-0 && check_D3<=-0))
		return true;
	else
		return false;

}

void accum_surface_area(grid_pointIter point1,grid_pointIter point2,CVector3d intersect_normal){
	if(point1->grid_x != point2->grid_x){
		if ((point1->grid_y != m_y_num-1) && (point1->grid_z != m_z_num-1))
		{
			m_surface_area += fabs(intersect_normal[0])*m_area_unit;
			return;
		}
	}
	if(point1->grid_y != point2->grid_y){
		if ((point1->grid_x != m_x_num-1) && (point1->grid_z != m_z_num-1))
		{
			m_surface_area += fabs(intersect_normal[1])*m_area_unit;
			return;
		}
	}	
	if(point1->grid_z != point2->grid_z){
		if ((point1->grid_x != m_x_num-1) && (point1->grid_y != m_y_num-1))
		{
			m_surface_area += fabs(intersect_normal[2])*m_area_unit;
			return;
		}
	}
}

 int read_pqr_and_classify_inte_points(std::string filename)
 {
   //read
   std::vector<int> atom_amino_type; 
   std::fstream file;
   file.open(filename.c_str(), std::ios::in);
   std::string line;
   getline(file, line);

   while(getline(file, line))
     {
       std::stringstream ss(line);
       std::string trash;
       int type;
       for(int i=0; i<4; i++)
	 ss>>trash;
       ss>>type;
       atom_amino_type.push_back(type);
     }
   file.close();

   //classify
   std::vector<int> inte_amino_type;
   for(int i=0; i<m_intersection_points.size(); i++)
     {
       if(m_intersection_types[i]==0)
	 inte_amino_type.push_back(atom_amino_type[(m_intersection_cell[i].current_atom)->index]);
       else if(m_intersection_types[i]==1)
	 {
	   double d_0=pow((m_intersection_points[i][0]-m_atoms_vec[m_intersection_cell[i].current_torus->atom_i->index]->center[0]),2)+pow((m_intersection_points[i][1]-m_atoms_vec[m_intersection_cell[i].current_torus->atom_i->index]->center[1]),2)+pow((m_intersection_points[i][2]-m_atoms_vec[m_intersection_cell[i].current_torus->atom_i->index]->center[2]),2)-pow(m_atoms_vec[m_intersection_cell[i].current_torus->atom_i->index]->r,2);
	   double d_1=pow((m_intersection_points[i][0]-m_atoms_vec[m_intersection_cell[i].current_torus->atom_j->index]->center[0]),2)+pow((m_intersection_points[i][1]-m_atoms_vec[m_intersection_cell[i].current_torus->atom_j->index]->center[1]),2)+pow((m_intersection_points[i][2]-m_atoms_vec[m_intersection_cell[i].current_torus->atom_j->index]->center[2]),2)-pow(m_atoms_vec[m_intersection_cell[i].current_torus->atom_j->index]->r,2);
	   if(d_0<d_1)
	     inte_amino_type.push_back(atom_amino_type[m_intersection_cell[i].current_torus->atom_i->index]);
	   else
	     inte_amino_type.push_back(atom_amino_type[m_intersection_cell[i].current_torus->atom_j->index]);   	     
	 }
       else
	 {
	   double d_0=pow((m_intersection_points[i][0]-m_atoms_vec[m_intersection_cell[i].current_concave->index_i]->center[0]),2)+pow((m_intersection_points[i][1]-m_atoms_vec[m_intersection_cell[i].current_concave->index_i]->center[1]),2)+pow((m_intersection_points[i][2]-m_atoms_vec[m_intersection_cell[i].current_concave->index_i]->center[2]),2)-pow(m_atoms_vec[m_intersection_cell[i].current_concave->index_i]->r,2);
	   double d_1=pow((m_intersection_points[i][0]-m_atoms_vec[m_intersection_cell[i].current_concave->index_j]->center[0]),2)+pow((m_intersection_points[i][1]-m_atoms_vec[m_intersection_cell[i].current_concave->index_j]->center[1]),2)+pow((m_intersection_points[i][2]-m_atoms_vec[m_intersection_cell[i].current_concave->index_k]->center[2]),2)-pow(m_atoms_vec[m_intersection_cell[i].current_concave->index_j]->r,2);
	   double d_2=pow((m_intersection_points[i][0]-m_atoms_vec[m_intersection_cell[i].current_concave->index_k]->center[0]),2)+pow((m_intersection_points[i][1]-m_atoms_vec[m_intersection_cell[i].current_concave->index_k]->center[1]),2)+pow((m_intersection_points[i][2]-m_atoms_vec[m_intersection_cell[i].current_concave->index_k]->center[2]),2)-pow(m_atoms_vec[m_intersection_cell[i].current_concave->index_k]->r,2);

	   if(d_0<d_1 && d_0<d_2)
	     inte_amino_type.push_back(atom_amino_type[m_intersection_cell[i].current_concave->index_i]);
	   else if(d_1<d_0 && d_1<d_2)
	     inte_amino_type.push_back(atom_amino_type[m_intersection_cell[i].current_concave->index_j]);
	   else
	     inte_amino_type.push_back(atom_amino_type[m_intersection_cell[i].current_concave->index_k]);
				       
	 }
     }
  
   //output result
   file.open("intersection_amino.txt", std::ios::out);
   for(int i=0; i<inte_amino_type.size(); i++)
     file<<inte_amino_type[i]<<std::endl;
   file.close();

   return 1;
 }

int partition_area(){
   m_partition_area.assign(m_atoms_vec.size(),0);
   for(int i=0; i<m_intersection_points.size(); i++)
     {
       if(m_intersection_types[i]==0)
	 {
	   /*	   
	   if(m_intersection_interior[i]->grid_x!=m_intersection_outside[i]->grid_x)
	     m_partition_area[(m_intersection_cell[i].current_atom)->index]+=m_area_unit*fabs(m_intersection_normals[i][0]);
	   else if(m_intersection_interior[i]->grid_y!=m_intersection_outside[i]->grid_y)
	     m_partition_area[(m_intersection_cell[i].current_atom)->index]+=m_area_unit*fabs(m_intersection_normals[i][1]);
	   else
	     m_partition_area[(m_intersection_cell[i].current_atom)->index]+=m_area_unit*fabs(m_intersection_normals[i][2]);
	   */
	 }
       else if(m_intersection_types[i]==1)
	 {
	   
	   double d_0=pow((m_intersection_points[i][0]-m_atoms_vec[m_intersection_cell[i].current_torus->atom_i->index]->center[0]),2)+pow((m_intersection_points[i][1]-m_atoms_vec[m_intersection_cell[i].current_torus->atom_i->index]->center[1]),2)+pow((m_intersection_points[i][2]-m_atoms_vec[m_intersection_cell[i].current_torus->atom_i->index]->center[2]),2)-pow(m_atoms_vec[m_intersection_cell[i].current_torus->atom_i->index]->r,2);
	   double d_1=pow((m_intersection_points[i][0]-m_atoms_vec[m_intersection_cell[i].current_torus->atom_j->index]->center[0]),2)+pow((m_intersection_points[i][1]-m_atoms_vec[m_intersection_cell[i].current_torus->atom_j->index]->center[1]),2)+pow((m_intersection_points[i][2]-m_atoms_vec[m_intersection_cell[i].current_torus->atom_j->index]->center[2]),2)-pow(m_atoms_vec[m_intersection_cell[i].current_torus->atom_j->index]->r,2);
	   
	   if(m_intersection_interior[i]->grid_x!=m_intersection_outside[i]->grid_x)
	     {      
	       m_partition_area[m_intersection_cell[i].current_torus->atom_i->index]+=m_area_unit*fabs(m_intersection_normals[i][0])*(d_1/(d_0+d_1));
	       m_partition_area[m_intersection_cell[i].current_torus->atom_j->index]+=m_area_unit*fabs(m_intersection_normals[i][0])*(d_0/(d_0+d_1));
	       }
	   else if(m_intersection_interior[i]->grid_y!=m_intersection_outside[i]->grid_y)
	     {
	       m_partition_area[m_intersection_cell[i].current_torus->atom_i->index]+=m_area_unit*fabs(m_intersection_normals[i][1])*(d_1/(d_0+d_1));
	       m_partition_area[m_intersection_cell[i].current_torus->atom_j->index]+=m_area_unit*fabs(m_intersection_normals[i][1])*(d_0/(d_0+d_1));
	       }
	   else
	     {
	       m_partition_area[m_intersection_cell[i].current_torus->atom_i->index]+=m_area_unit*fabs(m_intersection_normals[i][2])*(d_1/(d_0+d_1));
	       m_partition_area[m_intersection_cell[i].current_torus->atom_j->index]+=m_area_unit*fabs(m_intersection_normals[i][2])*(d_0/(d_0+d_1));
	       }
	   
	 }
       else
	 {
	   /*
	   double d_0=pow((m_intersection_points[i][0]-m_atoms_vec[m_intersection_cell[i].current_concave->index_i]->center[0]),2)+pow((m_intersection_points[i][1]-m_atoms_vec[m_intersection_cell[i].current_concave->index_i]->center[1]),2)+pow((m_intersection_points[i][2]-m_atoms_vec[m_intersection_cell[i].current_concave->index_i]->center[2]),2)-pow(m_atoms_vec[m_intersection_cell[i].current_concave->index_i]->r,2);
	   double d_1=pow((m_intersection_points[i][0]-m_atoms_vec[m_intersection_cell[i].current_concave->index_j]->center[0]),2)+pow((m_intersection_points[i][1]-m_atoms_vec[m_intersection_cell[i].current_concave->index_j]->center[1]),2)+pow((m_intersection_points[i][2]-m_atoms_vec[m_intersection_cell[i].current_concave->index_k]->center[2]),2)-pow(m_atoms_vec[m_intersection_cell[i].current_concave->index_j]->r,2);
	   double d_2=pow((m_intersection_points[i][0]-m_atoms_vec[m_intersection_cell[i].current_concave->index_k]->center[0]),2)+pow((m_intersection_points[i][1]-m_atoms_vec[m_intersection_cell[i].current_concave->index_k]->center[1]),2)+pow((m_intersection_points[i][2]-m_atoms_vec[m_intersection_cell[i].current_concave->index_k]->center[2]),2)-pow(m_atoms_vec[m_intersection_cell[i].current_concave->index_k]->r,2);

	   if(m_intersection_interior[i]->grid_x!=m_intersection_outside[i]->grid_x)
	     {
	       m_partition_area[m_intersection_cell[i].current_concave->index_i]+=m_area_unit*fabs(m_intersection_normals[i][0])*(1/d_0/(1/d_0+1/d_1+1/d_2));
	       m_partition_area[m_intersection_cell[i].current_concave->index_j]+=m_area_unit*fabs(m_intersection_normals[i][0])*(1/d_1/(1/d_0+1/d_1+1/d_2));
	       m_partition_area[m_intersection_cell[i].current_concave->index_k]+=m_area_unit*fabs(m_intersection_normals[i][0])*(1/d_2/(1/d_0+1/d_1+1/d_2));
	     }
	   else if(m_intersection_interior[i]->grid_y!=m_intersection_outside[i]->grid_y)
	     {
	       m_partition_area[m_intersection_cell[i].current_concave->index_i]+=m_area_unit*fabs(m_intersection_normals[i][1])*(1/d_0/(1/d_0+1/d_1+1/d_2));
	       m_partition_area[m_intersection_cell[i].current_concave->index_j]+=m_area_unit*fabs(m_intersection_normals[i][1])*(1/d_1/(1/d_0+1/d_1+1/d_2));
	       m_partition_area[m_intersection_cell[i].current_concave->index_k]+=m_area_unit*fabs(m_intersection_normals[i][1])*(1/d_2/(1/d_0+1/d_1+1/d_2));
	     }
	   else
	     {
	       m_partition_area[m_intersection_cell[i].current_concave->index_i]+=m_area_unit*fabs(m_intersection_normals[i][2])*(1/d_0/(1/d_0+1/d_1+1/d_2));
	       m_partition_area[m_intersection_cell[i].current_concave->index_j]+=m_area_unit*fabs(m_intersection_normals[i][2])*(1/d_1/(1/d_0+1/d_1+1/d_2));
	       m_partition_area[m_intersection_cell[i].current_concave->index_k]+=m_area_unit*fabs(m_intersection_normals[i][2])*(1/d_2/(1/d_0+1/d_1+1/d_2));
	     }
	   */
	 }
     }
   return 1;
}

int process_each_grid_edge(grid_pointIter point_inside, grid_pointIter point_outside,std::vector<CVector3d> &intersect_points,
			   std::vector<CVector3d> &intersect_normal,std::vector<int> &intersect_type,std::vector<cell_unit> &intersect_cell){

	if (!point_outside->is_irregular) {
	  //int x_pos = point_outside->grid_x;
	  //int y_pos = point_outside->grid_y;
	  //int z_pos = point_outside->grid_z;
	  //if (m_grid_regularity[x_pos + y_pos*m_x_num + z_pos*m_x_num*m_y_num] == true)
	  //{
	  //m_grid_regularity[x_pos + y_pos*m_x_num + z_pos*m_x_num*m_y_num] = false;
	  point_outside->is_irregular = true;
	  //m_surface_volume += m_volume_unit / 2;
	  //}
	}
	if (!point_inside->is_irregular) {
	  //int x_pos = point_inside->grid_x;
	  //int y_pos = point_inside->grid_y;
	  //int z_pos = point_inside->grid_z;
	  //if (m_grid_regularity[x_pos + y_pos*m_x_num + z_pos*m_x_num*m_y_num] == true)
	  //{
	  //m_grid_regularity[x_pos + y_pos*m_x_num + z_pos*m_x_num*m_y_num] = false;
	  point_inside->is_irregular = true;
	  //m_surface_volume += m_volume_unit / 2;
	  //}
	}
	double dist = -1; int f_id = -1;
	CVector3d farthest_point, farthest_normal;
	int farthest_type;
	cell_unit farthest_cell;
	//the intersection point with the surface is the one farthest to the interior point
	for(int i=0;i<intersect_points.size(); i++){
		double current_dist = (intersect_points[i] - point_inside->grid_pos).Length();
		if(current_dist > dist){
			dist = current_dist;
			f_id = i;
		}
	}

	m_intersection_interior.push_back(point_inside);
	m_intersection_outside.push_back(point_outside);


	if(intersect_points.size()>0){
		farthest_point = intersect_points[f_id];
		farthest_normal = intersect_normal[f_id];
		farthest_type = intersect_type[f_id];
		farthest_cell = intersect_cell[f_id];

		m_intersection_points.push_back(farthest_point);
		m_intersection_normals.push_back(farthest_normal);
		m_intersection_types.push_back(farthest_type);
		m_intersection_cell.push_back(farthest_cell);

		accum_surface_area(point_inside,point_outside,farthest_normal);

		return 0;
	}
	else{
		std::cout<<"missing case: "<<std::endl<<point_inside->grid_x<<' '<<point_inside->grid_y<<
			' '<<point_inside->grid_z<<std::endl;

		std::cout<<point_outside->grid_x<<' '<<point_outside->grid_y<<
			' '<<point_outside->grid_z<<std::endl;

		//std::cout<<"#####################################################################"<<std::endl;

		farthest_point = (point_inside->grid_pos+point_outside->grid_pos)/2;
		farthest_normal = CVector3d(1.0,0.0,0.0);
		farthest_type = -1;
		
		m_intersection_points.push_back(farthest_point);
		m_intersection_normals.push_back(farthest_normal);
		m_intersection_types.push_back(farthest_type);
		m_intersection_cell.push_back(farthest_cell);
		return 1;
	}
}


 void output_intersection(bool& first_write, int a, int b, int c, int block_size, int idx)
{
	//output the intersection info (first the location of interior grid, location of outside grid, then location and normal of the intersection)
  std::ofstream file_intersection;
  //std::ofstream file_poisson;
  //std::ofstream file_block;
  if(first_write==true)
    {
      file_intersection.open("intersection_info.txt", ios::out);
      //file_poisson.open("poisson.xyzn", ios::out);
	  //file_block.open("block.txt", ios::out);
      first_write=false;
    }
  else
    {
      file_intersection.open("intersection_info.txt", ios::app);
      //file_poisson.open("poisson.xyzn", ios::app);
	  //file_block.open("block.txt", ios::app);
    }
  for (int i = 0; i<m_intersection_points.size(); i++) {
    //intersection file----------------------------------------
    if((m_intersection_interior[i]->grid_x==a*block_size+m_x_num-1 && m_intersection_outside[i]->grid_x==a*block_size+m_x_num-1) 
       || (m_intersection_interior[i]->grid_y==b*block_size+m_y_num-1 && m_intersection_outside[i]->grid_y==b*block_size+m_y_num-1)
       || (m_intersection_interior[i]->grid_z==c*block_size+m_z_num-1 && m_intersection_outside[i]->grid_z==c*block_size+m_z_num-1))
      continue;
    else
    {
      file_intersection << m_intersection_interior[i]->grid_x << ' ' << m_intersection_interior[i]->grid_y << ' ' << m_intersection_interior[i]->grid_z << ' ';
      file_intersection << m_intersection_outside[i]->grid_x << ' ' << m_intersection_outside[i]->grid_y << ' ' << m_intersection_outside[i]->grid_z << ' ';
      file_intersection << std::scientific;
      file_intersection << std::setprecision(12) << m_intersection_points[i][0] << ' ' << m_intersection_points[i][1] << ' ' << m_intersection_points[i][2] << ' ';
      file_intersection << std::setprecision(12) << m_intersection_normals[i][0] << ' ' << m_intersection_normals[i][1] << ' ' << m_intersection_normals[i][2] << ' ';
      
	  if (m_intersection_types[i] == 0)
		  file_intersection << m_local_to_global_atom_idx[(m_intersection_cell[i].current_atom)->index];
      if (m_intersection_types[i] == 1) {
		TorusIter current_torus = m_intersection_cell[i].current_torus;
		file_intersection << m_local_to_global_atom_idx[current_torus->atom_i->index] << ' ' 
			<< m_local_to_global_atom_idx[current_torus->atom_j->index];
      }
      if (m_intersection_types[i] == 2) {
		Concave_SphereIter current_concave = m_intersection_cell[i].current_concave;
		file_intersection << m_local_to_global_atom_idx[current_concave->index_i] << ' '
			<< m_local_to_global_atom_idx[current_concave->index_j] << ' '
			<< m_local_to_global_atom_idx[current_concave->index_k];
      }
    }
    file_intersection << std::endl;

	//// for block info
	//if (m_intersection_interior[i]->grid_x%block_size == 0 
	//	&& m_intersection_outside[i]->grid_x%block_size == 0)
	//	file_block << "-1" << std::endl;
	//else if (m_intersection_interior[i]->grid_y%block_size == 0
	//	&& m_intersection_outside[i]->grid_y%block_size == 0)
	//	file_block << "-1" << std::endl;
	//else if (m_intersection_interior[i]->grid_z%block_size == 0
	//	&& m_intersection_outside[i]->grid_z%block_size == 0)
	//	file_block << "-1" << std::endl;
	//else
	//	file_block << idx << std::endl;
		
    // for poisson reconstruction ----------------------------------------------------------
    //file_poisson << std::setprecision(12) << m_intersection_points[i][0] << ' ' << m_intersection_points[i][1] << ' ' << m_intersection_points[i][2] << ' ';
    //file_poisson << std::setprecision(12) << m_intersection_normals[i][0] << ' ' << m_intersection_normals[i][1] << ' ' << m_intersection_normals[i][2] << ' ';
    
    //file_poisson << std::endl;
  }
  
  file_intersection.close();
  //file_poisson.close();
}
 


//given the status of the each grid point (inside and outside), check the missing intersections
void test_intersection(){
	std::cout<<m_atoms_vec.size()<<" atoms..."<<std::endl;
	std::cout<<m_torus_vec.size()<<" torus..."<<std::endl;
	std::cout<<m_concave_spheres.size()<<" concave faces..."<<std::endl;
	std::cout<<"bounding box..."<<std::endl;
	std::cout<<m_min_x<<","<<m_min_y<<","<<m_min_z<<std::endl;
	std::cout<<m_max_x<<","<<m_max_y<<","<<m_max_z<<std::endl;
	//first find the possible intersections
	std::vector<grid_pointIter> possible_interior;
	std::vector<grid_pointIter> possible_outside;

	for(int i=0;i<m_grid_points.size();i++){
		grid_pointIter grid_point = m_grid_points[i];
		if(grid_point->is_inside){
			if(grid_point->grid_x >0){
				grid_pointIter  backward_x_point = m_grid_points[grid_index(grid_point->grid_x-1,grid_point->grid_y,grid_point->grid_z)];

				if(!backward_x_point->is_inside){
					possible_interior.push_back(grid_point);
					possible_outside.push_back(backward_x_point);
				}
			}
			if(grid_point->grid_y >0){
				grid_pointIter  backward_y_point = m_grid_points[grid_index(grid_point->grid_x,grid_point->grid_y-1,grid_point->grid_z)];

				if(!backward_y_point->is_inside){
					possible_interior.push_back(grid_point);
					possible_outside.push_back(backward_y_point);
				}
			}
			if(grid_point->grid_z >0){
				grid_pointIter  backward_z_point = m_grid_points[grid_index(grid_point->grid_x,grid_point->grid_y,grid_point->grid_z-1)];

				if(!backward_z_point->is_inside){
					possible_interior.push_back(grid_point);
					possible_outside.push_back(backward_z_point);
				}
			}

			if(grid_point->grid_x <m_x_num-1){
				grid_pointIter  forward_x_point = m_grid_points[grid_index(grid_point->grid_x+1,grid_point->grid_y,grid_point->grid_z)];

				if(!forward_x_point->is_inside){
					possible_interior.push_back(grid_point);
					possible_outside.push_back(forward_x_point);
				}
			}

			if(grid_point->grid_y <m_y_num-1){
				grid_pointIter  forward_y_point = m_grid_points[grid_index(grid_point->grid_x,grid_point->grid_y+1,grid_point->grid_z)];

				if(!forward_y_point->is_inside){
					possible_interior.push_back(grid_point);
					possible_outside.push_back(forward_y_point);
				}
			}

			if(grid_point->grid_z <m_z_num-1){
				grid_pointIter  forward_z_point = m_grid_points[grid_index(grid_point->grid_x,grid_point->grid_y,grid_point->grid_z+1)];

				if(!forward_z_point->is_inside){
					possible_interior.push_back(grid_point);
					possible_outside.push_back(forward_z_point);
				}
			}
		}
	}

	int num_found = 0;
	//now scan all the intersection to make sure there are no 'missing' here
	for(int i=0;i<m_intersection_points.size();i++){
		for(int j=0;j<possible_interior.size();j++){
			if((possible_interior[j]->grid_x==m_intersection_interior[i]->grid_x)
				&& (possible_interior[j]->grid_y==m_intersection_interior[i]->grid_y)
				&& (possible_interior[j]->grid_z==m_intersection_interior[i]->grid_z)
				&& (possible_outside[j]->grid_x==m_intersection_outside[i]->grid_x)
				&& (possible_outside[j]->grid_y==m_intersection_outside[i]->grid_y)
				&& (possible_outside[j]->grid_z==m_intersection_outside[i]->grid_z)){
					possible_interior.erase(possible_interior.begin()+j);
					possible_outside.erase(possible_outside.begin()+j);
					num_found++;
					break;
			}
		}
	}

	for(int i=0;i<(std::min((int)(possible_interior.size()),5));i++){
		std::cout<<std::endl<<possible_interior[i]->grid_x<<' '<<possible_interior[i]->grid_y<<' '<<possible_interior[i]->grid_z<<std::endl;
		std::cout<<possible_outside[i]->grid_x<<' '<<possible_outside[i]->grid_y<<' '<<possible_outside[i]->grid_z<<std::endl;
	}
	std::cout<<num_found<<" number of points found..."<<std::endl;
	std::cout<<possible_interior.size()<<" number of points missing..."<<std::endl;

}


};

#endif
