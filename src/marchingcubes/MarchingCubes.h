#ifndef _MarchingCubes_H_
#define _MarchingCubes_H_

//#define FLT_EPSILON 1.19209290E-07F
#define DEBUG

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include "LookUpTable.h"
#include "../intersection/parallel_wrapper.h"

//Vertex structure
typedef struct
{
  double x, y, z; //Vertex coordinates
  double nx, ny, nz; //Vertex normal
} Vertex;

//Triangle structure
typedef struct
{
  int v1, v2, v3; //Triangle vertex indices
} Triangle;


//Marching Cube Algorithm
class MarchingCubes
{
 public:
  //Constructor
  MarchingCubes()
    {
      _cube.resize(8);
      _nvertices=0;
      _ntriangles=0;
    };
  MarchingCubes(parallel_wrapper* pw)
    {
      _pw = pw;
      _cube.resize(8);
      _nvertices=0;
      _ntriangles=0;
    };
  /*
  //Destructor
  ~MarchingCubes()
    {
      std::cout<<"in destructor"<<std::endl;
      
      std::cout<<_vertices.size()<<std::endl;
      std::cout<<_triangles.size()<<std::endl;
      _vertices.clear();
      _triangles.clear();
    };
*/
  //Get # of vertices in the generated mesh
  inline const int nvertices()
  {
    return _nvertices;
  }
  
  //Get # of triangles in the generated mesh
  inline const int ntriangles()
  {
    return _ntriangles;
  }

  //Get a specific vertex in the generated mesh
  inline Vertex vertex(const int i) 
  {
    if(i<0 || i>=_nvertices)
      {
	Vertex v;
	v.x=0;
	v.y=0;
	v.z=0;
	v.nx=0;
	v.ny=0;
	v.nz=0;
	std::cout<<"Error: Vertex index out of range!"<<std::endl;
	return v;
      } 
    else
      return _vertices[i];
  }

  //Get a specific triangle in the generated mesh
  inline Triangle triangle(const int i)
  {
    if(i<0 || i>=_ntriangles)
      {
	Triangle t;
	t.v1=0;
	t.v2=0;
	t.v3=0;
	std::cout<<"Error: Triangle index out of range!"<<std::endl;
	return t;
      } 
    else
      return _triangles[i];
  }

  //Get the vertex buffer of the generated mesh
  inline std::vector<Vertex>& vertex_buffer()
  {
    return _vertices;
  }

  //Get the triangle buffer of the generated mesh
  inline std::vector<Triangle>& triangle_buffer()
  {
    return _triangles;
  }

    int process();
    int set_bounding_box();
    int set_grid();
    int set_intersection();

  //Read bounding box and the grid size
  int read_bounding_box(std::string filename);

  //Read grid info
  int read_grid_info(std::string filename);

  //Read intersection info
  int read_intersection_info(std::string filename);

  //Save grid info into data
  int save_grid_info();

  //Get grid info from data
  inline double get_data(const int i, const int j, const int k)
  {
    return _data[i+j*_x_size+k*_x_size*_y_size];
  };
  
  //initialize _x_verts, _y_verts and _z_verts
  int initialize_verts();

  //Add a vertex on the horizontal edge 
  int add_x_vert(const int i);

  //Add a vertex on the longitudinal edge 
  int add_y_vert(const int i);
  
  //Add a vertex on the vertical edge 
  int add_z_vert(const int i);
  
  //Add a vertex inside the current cube
  int add_c_vert();

  //Set the vertex index on the lower horizontal edge of a specific cube
  int set_x_vert( int val, const int i, const int j, const int k);
  
  //Set the vertex index on the lower longitudinal edge of a specific cube
  int set_y_vert(const int val, const int i, const int j, const int k);
  
  //Set the vertex index on the lower vertical edge of a specific cube
  int set_z_vert(const int val, const int i, const int j, const int k);

  //Get the vertex index on the lower horizontal edge of a specific cube
  inline int get_x_vert(const int i, const int j, const int k)
  {
    return _x_verts[i+j*_x_size+k*_x_size*_y_size];
  };

  //Get the vertex index on the lower longitudinal edge of a specific cube
  inline int get_y_vert(const int i, const int j, const int k)
  {
    return _y_verts[i+j*_x_size+k*_x_size*_y_size];
  };

  //Get the vertex index on the lower vertical edge of a specific cube
  inline int get_z_vert(const int i, const int j, const int k)
  {
    return  _z_verts[i+j*_x_size+k*_x_size*_y_size];
  };

  //Save intersection vertices into _nvertices
  int save_intersection_info();

  //Test if the components of the tesselation of the cube should be connected by the interior of an ambiguous face
  int test_face(signed char face);

  //Test if the components of the tesselation of the cube should be connected through the interior of the cube
  bool test_interior(signed char s)    ;

  //add a triangle to the mesh
  int add_triangle(const char* trig, int n, int v12=-1);
  
  //Tessellate one cube
  int process_cube();

  //run
  int run();

  //calculate triangle normal
  int triangle_normal();

  //output the wavefront OBJ file
  int write_OBJ_file();
  int write_OFF_file();
  int write_VTK_file();

#ifdef DEBUG
 public:
#else
 private:
#endif
    parallel_wrapper* _pw;

  std::vector<std::vector<double> > _boun_box;
  int _x_size, _y_size, _z_size;
  std::vector<std::vector<double> > _grid_info;
  std::vector<std::vector<int> > _inte_dual;
  std::vector<std::vector<double> > _inte_verts;
  std::vector<std::vector<double> > _inte_normal;
  
  int _nvertices;
  int _ntriangles;
  std::vector<Vertex> _vertices;
  std::vector<Triangle> _triangles;
  std::vector<std::vector<double> > _triangle_normal;
  std::vector<double> _data;

  std::vector<int> _x_verts;
  std::vector<int> _y_verts;
  std::vector<int> _z_verts;

  int _i;
  int _j;
  int _k;
  
  std::vector<double> _cube;
  unsigned char _lut_entry;
  unsigned char _case;
  unsigned char _config;
  unsigned char _subconfig;
};


#endif 
