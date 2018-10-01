/********************************
March/07/2014

This project is for computing the intersection of the regular grid
with the union of spheres/atoms

Input: (stored in a file)

number of atoms; (N)
center of each atom;
radius of each atom; (r)
bounding box; (min_x,min_y,min_z;max_x,max_y,max_z)
number of grid points in each dimension; (num_x,num_y,num_z)

Output:
number of intersection points;
location of each intersection point
normal direction of each intersection point


ATTENTION: for now we assume that each grid edge intesects with the surface at most ONCE

*********************************/
#ifndef _VANDERWAALS_SURFACE_
#define _VANDERWAALS_SURFACE_

#include <iostream>
#include <fstream> 
#include <string>
#include <vector>
#include <cmath>
#include "Vector3d.h"
#include "Atom_sphere.h"


class VanderWaals_surface{
public:
	VanderWaals_surface() {

	}
	int m_N_atoms;  //# of atoms
	Atom *m_atoms_list;
	double m_min_x,m_min_y,m_min_z;  //bounding box
	double m_max_x,m_max_y,m_max_z;
	int m_x_num,m_y_num,m_z_num; //number of grid points in each dimension
	double m_x_step,m_y_step,m_z_step;
	std::vector<std::vector<int>> m_grid_point; //stores the id of atom sphere which includes corresponding the grid point

	std::vector<CVector3d> m_intersection_point; //store the output: location of intersection points
	std::vector<CVector3d> m_intersection_normal; //store the normal direction of the inteserction points



	//read the input variables from file
bool read_info(char *filename){
	std::ifstream file(filename);
	if(!file.is_open()){
		std::cout<<"cannot find input file!"<<std::endl;
		return false;
	}

	std::string line;
	std::getline(file,line);
	file>>m_N_atoms;  //read the number of atoms
	m_atoms_list = new Atom[m_N_atoms];

	std::getline(file,line);
	std::getline(file,line);
	//read the atoms info
	for(int i=0;i<m_N_atoms;i++){
		double x,y,z,r;
		file>>x>>y>>z>>r;
	/*	m_atoms_list[i].center = CVector3d(x,y,z);
		m_atoms_list[i].r = r;*/

		m_atoms_list[i] = Atom(CVector3d(x,y,z),r,i);
	}

	//read the bounding box info
	std::getline(file,line);
	std::getline(file,line);
	file>>m_min_x>>m_min_y>>m_min_z;
	std::getline(file,line);
	std::getline(file,line);
	file>>m_max_x>>m_max_y>>m_max_z;
	std::getline(file,line);
	std::getline(file,line);
	file>>m_x_num>>m_y_num>>m_z_num;
	file.close();

	m_x_step = (m_max_x - m_min_x)/(m_x_num-1);
	m_y_step = (m_max_y - m_min_y)/(m_y_num-1);
	m_z_step = (m_max_z - m_min_z)/(m_z_num-1);

	m_grid_point.resize(m_x_num*m_y_num*m_z_num);
	
	return true;

}

/*******
check which grid point in included in this atom sphere

say current the gird point represented as (X,Y,Z) [index starts from 0]

then its position in 3D is (X*step_x+m_min_x, Y*step_y+m_min_y, Z*step_z+m_min_z)
where 

step_x = (m_max_x - m_min_x)/(m_x_num-1)
step_y = (m_max_y - m_min_y)/(m_y_num-1)
step_z = (m_max_z - m_min_z)/(m_z_num-1)

also the corresponding index in m_grid_point is X + m_x_num*Y + (m_x_num*m_y_num)*Z
******/

//return the global index of the grid point
int grid_index(int X,int Y, int Z){
	return (X + m_x_num*Y + m_x_num*m_y_num*Z);
}

//return the location of grid point in 3D
CVector3d grid_pos(int X,int Y,int Z){
	return CVector3d(X * m_x_step + m_min_x,Y * m_y_step + m_min_y,Z * m_z_step + m_min_z);

}

//check whether current gird point is the atom sphere or not
bool is_inside_atom(int atom_id,int grid_x,int grid_y,int grid_z){

	CVector3d gird_loc = grid_pos(grid_x,grid_y,grid_z);
	
	double radius = m_atoms_list[atom_id].r;

	double tmp = (gird_loc - m_atoms_list[atom_id].center).LengthSquared() - radius*radius;
	
	if(tmp>0)
		return false;
	else
		return true;

}



void process_each_atom(int atom_id){
	//decide which gird points are inside this atom sphere
	int tmp_min_x = int(ceil((m_atoms_list[atom_id].center[0] - m_atoms_list[atom_id].r - m_min_x)/m_x_step)+0.5);
	int tmp_min_y = int(ceil((m_atoms_list[atom_id].center[1] - m_atoms_list[atom_id].r - m_min_y)/m_y_step)+0.5);
	int tmp_min_z = int(ceil((m_atoms_list[atom_id].center[2] - m_atoms_list[atom_id].r - m_min_z)/m_z_step)+0.5);

	int tmp_max_x = int(floor((m_atoms_list[atom_id].center[0] + m_atoms_list[atom_id].r - m_min_x)/m_x_step)+0.5);
	int tmp_max_y = int(floor((m_atoms_list[atom_id].center[1] + m_atoms_list[atom_id].r - m_min_y)/m_y_step)+0.5);
	int tmp_max_z = int(floor((m_atoms_list[atom_id].center[2] + m_atoms_list[atom_id].r - m_min_z)/m_z_step)+0.5);

	//sanity check
	int min_x_id = (tmp_min_x>0 ? tmp_min_x:0);
	int min_y_id = (tmp_min_y>0 ? tmp_min_y:0);
	int min_z_id = (tmp_min_z>0 ? tmp_min_z:0);

	int max_x_id = (tmp_max_x<m_x_num ? tmp_max_x : m_x_num);
	int max_y_id = (tmp_max_y<m_y_num ? tmp_max_y : m_y_num);
	int max_z_id = (tmp_max_z<m_z_num ? tmp_max_z : m_z_num);

	for(int x_i=min_x_id; x_i<=max_x_id; x_i++)
		for(int y_i=min_y_id; y_i<=max_y_id; y_i++)
			for(int z_i=min_z_id; z_i<=max_z_id; z_i++)
				if(is_inside_atom(atom_id,x_i,y_i,z_i))
					m_grid_point[grid_index(x_i,y_i,z_i)].push_back(atom_id);

}

//compute the distance of the grid point to the atom center
double distance_atom(int atom_id,int grid_x,int grid_y,int grid_z){
	CVector3d grid_loc = grid_pos(grid_x,grid_y,grid_z);

	return (grid_loc-m_atoms_list[atom_id].center).Length();

}
//compute the intersection
//input: the index of points of the grid edge
void compute_intesection(int outside_x,int outside_y,int outside_z, int inside_x,int inside_y,int inside_z){

	double track_l= -1;
	CVector3d intersect_pos, intersect_normal;
	
	for(int i=0;i<m_grid_point[grid_index(inside_x,inside_y,inside_z)].size();i++){
		int track_atom = track_atom = m_grid_point[grid_index(inside_x,inside_y,inside_z)][i];
		/***
	compute the intersection point given the equation from wiki page:
	http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
		******/
		double r = m_atoms_list[track_atom].r;
		CVector3d o = grid_pos(inside_x,inside_y,inside_z);
		CVector3d c = m_atoms_list[track_atom].center;
		CVector3d l = (-grid_pos(inside_x,inside_y,inside_z) + grid_pos(outside_x,outside_y,outside_z));
		l.Normalize();

		double tmp1 = -l.Dot(o-c);
		double tmp2 = sqrt(tmp1*tmp1 - (o-c).LengthSquared() + r*r);
		double d1 = tmp1 + tmp2;
		double d2 = tmp1 - tmp2;

		if(track_l < d1){
			track_l = d1;
			//we know d1, d2 should be one positive and one negative, and the positive one (d1)is the intersection we are looking for
			intersect_pos = o + d1 * l;
			intersect_normal = intersect_pos - c;
			intersect_normal.Normalize();
		}

	}
	m_intersection_point.push_back(intersect_pos);
	m_intersection_normal.push_back(intersect_normal);
}
void traverse_each_grid_point(int x_i,int y_i,int z_i){
	
	if(m_grid_point[grid_index(x_i,y_i,z_i)].size()==0)
		return;

	//now (x_i, y_i, z_i) is inside
	//check whether (x_i-1, y_i, z_i) is outside
	if(x_i>0 && m_grid_point[grid_index(x_i-1,y_i,z_i)].size()==0)
		compute_intesection(x_i-1,y_i,z_i,x_i,y_i,z_i);
	
	//check whether (x_i, y_i-1, z_i) is outside
	if(y_i>0 && m_grid_point[grid_index(x_i,y_i-1,z_i)].size()==0)
		compute_intesection(x_i,y_i-1,z_i,x_i,y_i,z_i);


	//check whether (x_i, y_i, z_i-1) is outside
	if(z_i>0 && m_grid_point[grid_index(x_i,y_i,z_i-1)].size()==0)
		compute_intesection(x_i,y_i,z_i-1,x_i,y_i,z_i);
	

	//check whether (x_i+1, y_i, z_i) is outside
	if(x_i<m_x_num-1 && m_grid_point[grid_index(x_i+1,y_i,z_i)].size()==0)
		compute_intesection(x_i+1,y_i,z_i,x_i,y_i,z_i);
		

	//check whether (x_i, y_i+1, z_i) is outside
	if(y_i<m_y_num-1 && m_grid_point[grid_index(x_i,y_i+1,z_i)].size()==0)
		compute_intesection(x_i,y_i+1,z_i,x_i,y_i,z_i);

	//check whether (x_i, y_i, z_i+1) is outside
	if(z_i<m_z_num-1 && m_grid_point[grid_index(x_i,y_i,z_i+1)].size()==0)
		compute_intesection(x_i,y_i,z_i+1,x_i,y_i,z_i);
}


int process_vanderwaals_surface(char *filename){
		//first read the input from file 
	if(!read_info(filename))
		return 0;

	/******
	process each atom and store the info for each grid point:
	which atom sphere this sphere is in
	*********/
	for(int i=0;i<m_N_atoms;i++)
		process_each_atom(i);

	//process each grid edge with one grid point inside and the other outside
	for(int x_i=0;x_i<m_x_num;x_i++)
		for(int y_i=0;y_i<m_y_num;y_i++)
			for(int z_i=0;z_i<m_z_num;z_i++)
				traverse_each_grid_point(x_i,y_i,z_i);

	std::ofstream file("intersection.txt");
	//output the intersection points and normal
	for(int i=0;i<m_intersection_point.size();i++){
		std::cout<<i<<"-th intersection pos:"<<m_intersection_point[i].x()<<" "<<m_intersection_point[i].y()<<" "<<m_intersection_point[i].z()<<std::endl;
		std::cout<<i<<"-th intersection normal:"<<m_intersection_normal[i].x()<<" "<<m_intersection_normal[i].y()<<" "<<m_intersection_normal[i].z()<<std::endl<<std::endl;

		file<<i<<"-th intersection pos:"<<m_intersection_point[i].x()<<" "<<m_intersection_point[i].y()<<" "<<m_intersection_point[i].z()<<std::endl;
		file<<i<<"-th intersection normal:"<<m_intersection_normal[i].x()<<" "<<m_intersection_normal[i].y()<<" "<<m_intersection_normal[i].z()<<std::endl<<std::endl;

	}

	file.close();

	//clear space
	delete []m_atoms_list;
	m_grid_point.clear();
	m_intersection_normal.clear();
	m_intersection_point.clear();

	return 1;

}

};

#endif