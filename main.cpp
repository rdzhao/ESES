//#define TEST_VANDERWAALS

//grid intersection for molecular surface
/*********
input:
*.xyzr file
radius of the probe sphere
resolution of grid
extension of the bounding box

output:
1.grid_info.txt: 
index of the grid point in each dimension (X,Y,Z) and
whether it is inside (1) or outside (-1) the molecular surface

2.bounding_box.txt:
the bounding box (min_x,min_y,min_z \n max_x,max_y,max_z)
and the number of grid points in each dimension (X,Y,Z)

3.intersection_info.txt
output the info for intersection points:
the index of the inside grid point in each dimension;
the index of the outside grid point in each diemsnion;
location of the intersection points;
normal of the intersection;
the index/indice of the corresponding atoms
***********************/


#include <time.h>
#include <chrono>
#include <numeric>
#include <omp.h>
#include "molecular_surface.h"

//#include "viewer\Viewer.h"

////grid intersection for VanderWaals surface
///*********
//input: file which contains:)
//
//number of atoms; (N)
//center of each atom;
//radius of each atom; (r)
//bounding box; (min_x,min_y,min_z;max_x,max_y,max_z)
//number of grid points in each dimension; (num_x,num_y,num_z)
//
//******/

#define BLOCK_SIZE 63 //grid block-1 

class basic_info
{
public:
  double m_probe_radius;
  double m_extend_bounding;
  double m_grid_resolution;
  double m_area_unit;
  double m_volume_unit;
  int m_N_atoms;

  double m_x_step, m_y_step, m_z_step;

  std::vector<Atom *> m_atoms_vec;
  CVector3d m_mesh_min,m_mesh_max;
  
  double t_min_x, t_min_y, t_min_z;  //total bounding box
  double t_max_x, t_max_y, t_max_z;
  int t_x_num, t_y_num, t_z_num; //total number of grid points in each dimension
};
  
//read the input variables from file
  bool read_info(char *filename_xyzr, basic_info &bi, double probe_radius,double grid_size, double bounding_residual){
    bi.m_probe_radius = probe_radius;
    bi.m_extend_bounding = bounding_residual;
    bi.m_grid_resolution = grid_size;
    bi.m_area_unit = pow(bi.m_grid_resolution,2);
    bi.m_volume_unit = pow(bi.m_grid_resolution,3);
    
    bi.m_x_step = bi.m_y_step = bi.m_z_step = bi.m_grid_resolution;
    
    std::ifstream file(filename_xyzr);
    if(!file.is_open()){
      std::cout<<"cannot find input file for atom info!"<<std::endl;
      return false;
    }
    
    int i=0;
    while(!file.eof()){
      double x,y,z,r;
      file>>x>>y>>z>>r;
      if(file.fail())
	break;
      Atom *new_atom = new Atom(CVector3d(x,y,z),r,i);
      bi.m_atoms_vec.push_back(new_atom);
      
      if(i==0){
	bi.m_mesh_min = CVector3d(x-r, y-r, z-r);
	bi.m_mesh_max = CVector3d(x+r, y+r, z+r);
		}
      else{
	if( (x-r) < bi.m_mesh_min[0])
	  bi.m_mesh_min[0] = x-r;
	if( (y-r) < bi.m_mesh_min[1])
	  bi.m_mesh_min[1] = y-r;
	if( (z-r) < bi.m_mesh_min[2])
	  bi.m_mesh_min[2] = z-r;

	if( (x+r) > bi.m_mesh_max[0])
	  bi.m_mesh_max[0] = x+r;
	if( (y+r) > bi.m_mesh_max[1])
	  bi.m_mesh_max[1] = y+r;
	if( (z+r) > bi.m_mesh_max[2])
	  bi.m_mesh_max[2] = z+r;
      }
      i++;
    }
    file.close();
    bi.m_N_atoms = bi.m_atoms_vec.size();
    
    bi.t_min_x = bi.m_mesh_min[0] - bi.m_extend_bounding; 
    bi.t_min_y = bi.m_mesh_min[1] - bi.m_extend_bounding; 
    bi.t_min_z = bi.m_mesh_min[2] - bi.m_extend_bounding;
    bi.t_max_x = bi.m_mesh_max[0] + bi.m_extend_bounding; 
    bi.t_max_y = bi.m_mesh_max[1] + bi.m_extend_bounding; 
    bi.t_max_z = bi.m_mesh_max[2] + bi.m_extend_bounding;
    
    bi.t_x_num = int((bi.t_max_x - bi.t_min_x)/(bi.m_x_step)) + 1;
    bi.t_y_num = int((bi.t_max_y - bi.t_min_y)/(bi.m_y_step)) + 1;
    bi.t_z_num = int((bi.t_max_z - bi.t_min_z)/(bi.m_z_step)) + 1;
    
    bi.t_max_x = bi.t_min_x + (bi.t_x_num - 1)*bi.m_x_step;
    bi.t_max_y = bi.t_min_y + (bi.t_y_num - 1)*bi.m_y_step;
    bi.t_max_z = bi.t_min_z + (bi.t_z_num - 1)*bi.m_z_step;

    //m_min_x += ((m_x_num + 1) / 2)*m_x_step;
    //m_min_y += ((m_y_num + 1) / 2)*m_y_step;
    //m_min_z += ((m_z_num + 1) / 2)*m_z_step;
    
    //m_max_x = m_min_x + (((m_x_num - 1) / 2))*m_x_step;
    //m_max_y = m_min_y + (((m_y_num - 1) / 2))*m_y_step;
    //m_max_z = m_min_z + (((m_z_num - 1) / 2))*m_z_step;

    //m_x_num = (m_x_num+1) / 2;
    //m_y_num = (m_y_num+1) / 2;
    //m_z_num = (m_z_num+1) / 2;
    
    std::cout<<"grid dimension: "<<bi.t_x_num<<' '<<bi.t_y_num<<' '<<bi.t_z_num<<std::endl;
    
    return true;
  }


void initialize_molecular_surface(molecular_surface &ms, basic_info bi, int a, int b, int c)
{
  //critical: find the nearest cluster
  double extension=4+2*bi.m_probe_radius;
  
  double x_max, y_max, z_max;
  double x_min, y_min, z_min;

  ms.m_x_step=ms.m_y_step=ms.m_z_step=bi.m_grid_resolution;  

  ms.m_min_x=bi.t_min_x+a*BLOCK_SIZE*bi.m_grid_resolution;
  ms.m_min_y=bi.t_min_y+b*BLOCK_SIZE*bi.m_grid_resolution;
  ms.m_min_z=bi.t_min_z+c*BLOCK_SIZE*bi.m_grid_resolution;
  //ms.m_max_x=ms.m_min_x+(BLOCK_SIZE)*bi.m_grid_resolution;
  //ms.m_max_y=ms.m_min_y+(BLOCK_SIZE)*bi.m_grid_resolution;
  //ms.m_max_z=ms.m_min_z+(BLOCK_SIZE)*bi.m_grid_resolution;
  
if (ms.m_min_x + BLOCK_SIZE*ms.m_x_step > bi.t_max_x)
    {
      ms.m_x_num = int((bi.t_max_x - ms.m_min_x) / ms.m_x_step) + 1;
      ms.m_max_x = ms.m_min_x + (ms.m_x_num - 1)*ms.m_x_step;
    }
  else
    {
      ms.m_max_x = ms.m_min_x + BLOCK_SIZE*ms.m_x_step;
      ms.m_x_num = BLOCK_SIZE + 1;
    }
  
  if (ms.m_min_y + BLOCK_SIZE*ms.m_y_step > bi.t_max_y)
    {
      ms.m_y_num = int((bi.t_max_y - ms.m_min_y) / ms.m_y_step) + 1;
      ms.m_max_y = ms.m_min_y + (ms.m_y_num - 1)*ms.m_y_step;
    }
  else
    {
      ms.m_max_y = ms.m_min_y + BLOCK_SIZE*ms.m_y_step;
      ms.m_y_num = BLOCK_SIZE + 1;
    }
  
  if (ms.m_min_z + BLOCK_SIZE*ms.m_z_step > bi.t_max_z)
    {
      ms.m_z_num = int((bi.t_max_z - ms.m_min_z) / ms.m_z_step) + 1;
      ms.m_max_z = ms.m_min_z + (ms.m_z_num - 1)*ms.m_z_step;
    }
  else
    {
      ms.m_max_z = ms.m_min_z + BLOCK_SIZE*ms.m_z_step;
      ms.m_z_num = BLOCK_SIZE + 1;
    }
  cout << "Box info: "<< endl;
  cout << "Min: " << ms.m_min_x << " " << ms.m_min_y << " " << ms.m_min_z << endl;
  cout << "Max: " << ms.m_max_x << " " << ms.m_max_y << " " << ms.m_max_z << endl;
  cout << "Num: " << ms.m_x_num << " " << ms.m_y_num << " " << ms.m_z_num << endl;

  x_min=ms.m_min_x-extension;
  y_min=ms.m_min_y-extension;
  z_min=ms.m_min_z-extension;
  x_max=ms.m_max_x+extension;
  y_max=ms.m_max_y+extension;
  z_max=ms.m_max_z+extension;
  
  int idx=0;
  for(int i=0; i<bi.m_N_atoms; i++)
    {
      if((bi.m_atoms_vec[i]->center.x()>x_min && bi.m_atoms_vec[i]->center.x()<x_max) &&
	 (bi.m_atoms_vec[i]->center.y()>y_min && bi.m_atoms_vec[i]->center.y()<y_max) &&
	 (bi.m_atoms_vec[i]->center.z()>z_min && bi.m_atoms_vec[i]->center.z()<z_max))
	{
	  Atom *atom_ptr=new Atom(bi.m_atoms_vec[i]->center,bi.m_atoms_vec[i]->r,idx);;
	  ms.m_atoms_vec.push_back(atom_ptr);
	  ms.m_local_to_global_atom_idx[idx] = i;
	  idx++;
	}
    }
  
  ms.m_probe_radius=bi.m_probe_radius;
  ms.m_N_atoms=ms.m_atoms_vec.size();
  ms.m_grid_resolution=bi.m_grid_resolution;
  ms.m_surface_area=ms.m_surface_volume=0;
  ms.m_area_unit=bi.m_area_unit;
  ms.m_volume_unit=bi.m_volume_unit;
}

void compute_volume(double &total_volume, std::vector<bool> t_grid_status, std::vector<bool> t_grid_regularity, double m_volume_unit){
  //int total_inside=0;
  //int total_irregular=0;
  for(int i=0; i<t_grid_status.size(); i++)
    {
      if(t_grid_status[i]==true && t_grid_regularity[i]==true)
	total_volume+=m_volume_unit;
      else if(t_grid_regularity[i]==false)
	total_volume+=m_volume_unit/2;
    }
    //cout<<"------------------------total irregular:  "<<total_irregular<<endl;
  //cout<<"------------------------total inside:  "<<total_inside<<endl;
  /*for(int i=0;i<m_grid_status.size();i++){
    grid_pointIter grid_point = m_grid_points[i];
      if (grid_point->is_inside && !grid_point->is_irregular)
      {
      m_surface_volume += m_volume_unit;
      }
    int z_pos = int(i / (m_x_num*m_y_num));
    int y_pos = int((i - m_x_num*m_y_num*z_pos) / m_x_num);
    int x_pos = i - m_x_num*m_y_num*z_pos - m_x_num*y_pos;
    if(m_grid_status[i]==true && m_grid_fix_volume[i]==false)
      m_surface_volume += m_volume_unit;
  }*/
}

void output_info(std::vector<bool>& t_grid_status, vector<double>& t_partition_area, basic_info bi){
	//output each grid point location and whether it is inside (1) or outside (-1)
	std::ofstream file_grid("grid_info.txt");
	/*for(int i=0;i<m_grid_points.size();i++){
		grid_pointIter grid_point = m_grid_points[i];
		file_grid<<grid_point->grid_x<<' '<<grid_point->grid_y<<' '<<grid_point->grid_z<<' ';
		if(grid_point->is_inside)
			file_grid<<1<<std::endl;
		
		else
			file_grid<<-1<<std::endl;
	}*/
	for (int i = 0; i < t_grid_status.size(); i++)
	{
		int z_pos = int(i / (bi.t_x_num*bi.t_y_num));
		int y_pos = int((i - bi.t_x_num*bi.t_y_num*z_pos) / bi.t_x_num);
		int x_pos = i - bi.t_x_num*bi.t_y_num*z_pos - bi.t_x_num*y_pos;
		file_grid << x_pos<<' '<< y_pos << ' ' << z_pos << ' ';
		if(t_grid_status[i]==true)
		  file_grid << 1 << std::endl;
		else
		  file_grid << -1 << std::endl;
	}

	file_grid.close();

	//output the bounding box info for Cartesian mesh
	std::ofstream file_bounding("bounding_box.txt");
	file_bounding<<std::scientific;
	file_bounding<<std::setprecision(12)<<bi.t_min_x<<' '<<bi.t_min_y<<' '<<bi.t_min_z<<std::endl;
	file_bounding<<std::setprecision(12)<<bi.t_max_x<<' '<<bi.t_max_y<<' '<<bi.t_max_z<<std::endl;
	file_bounding<<bi.t_x_num<<' '<<bi.t_y_num<<' '<<bi.t_z_num<<std::endl;
	file_bounding.close();

	//output the partition area. atom index and corresponding area
	std::ofstream file_partition("partition_area.txt");
	for(int i=0; i<bi.m_N_atoms; i++)
	{
	  file_partition<<i<<" ";
	  //file_partition<<std::scientific;
	  file_partition<<t_partition_area[i]<<endl;
	}
	file_partition.close();

  int sum = accumulate(t_partition_area.begin(), t_partition_area.end(), 0.0);
  cout<<"Total area from partition: "<<sum<<endl;

}

void process_molecular_surface_newer(double probe_radius,double grid_size,double bounding_residual,char *filename_xyzr/*,int& numAtoms, double& area, double& volume, int& gx, int& gy, int& gz, int& bx, int& by, int& bz*/){
  //basic info
  basic_info bi;
  if(!read_info(filename_xyzr, bi, probe_radius, grid_size, bounding_residual))
    return;

  //read the input info
  //if(!ms.read_info(probe_radius,grid_size,bounding_residual,filename_xyzr))
  //return;
  
  int block_dim_x, block_dim_y, block_dim_z;
  block_dim_x = (bi.t_x_num - 2) / BLOCK_SIZE+1;
  block_dim_y = (bi.t_y_num - 2) / BLOCK_SIZE+1;
  block_dim_z = (bi.t_z_num - 2) / BLOCK_SIZE+1;
  cout << "Block info: " << block_dim_x << " " << block_dim_y << " " << block_dim_z << endl;
  
  // inside -- true
  std::vector<bool> t_grid_status;
  std::vector<bool> t_grid_regularity;
  t_grid_status.resize(bi.t_x_num*bi.t_y_num*bi.t_z_num);
  t_grid_regularity.resize(bi.t_x_num*bi.t_y_num*bi.t_z_num);
  for(int i=0; i<t_grid_regularity.size(); i++)
    t_grid_regularity[i]=true;
  
  vector<double> t_partition_area;
  t_partition_area.resize(bi.m_N_atoms, 0);
  double t_area_from_partition = 0;

  double total_surface_area, total_volume;
  total_surface_area=0;
  total_volume=0;

  omp_set_num_threads(8);
  bool first_write=true;
#pragma omp parallel for schedule(dynamic)

  //for (int c = 0; c < block_dim_z; c++)
  //for (int b = 0; b < block_dim_y; b++)
  //for (int a = 0; a < block_dim_x; a++)
  for(int index=0; index<block_dim_x*block_dim_y*block_dim_z; index++)	
  {
    int a,b,c;
    c=index/(block_dim_x*block_dim_y);
    b=(index-c*block_dim_x*block_dim_y)/block_dim_x;
    a=index-c*block_dim_x*block_dim_y-b*block_dim_x;


	  //cout << "--------------------------------------------------------" << endl;
	  cout << "Start block ("<<a<<", " << b << ", " << c << ")" << endl;

	  molecular_surface ms;
	  initialize_molecular_surface(ms,bi,a,b,c);
	  
	  std::cout << "initializing the surface..." << std::endl;
	  //initialize the molecular surface
	  ms.initialize_molecular_surface();
	  
	  int missing_intersection = 0;
	  //ms.m_surface_area = ms.m_surface_volume = 0;
		
	  //ms.m_grid_status.resize(ms.m_x_num*ms.m_y_num*ms.m_z_num);
	  //ms.m_grid_fix_volume.resize(ms.m_x_num*ms.m_y_num*ms.m_z_num);
	  //for (int i = 0; i < ms.m_grid_fix_volume.size(); i++)
	  //ms.m_grid_fix_volume[i] = false;
	  //ms.m_grid_regularity.resize(ms.m_x_num*ms.m_y_num*ms.m_z_num);
	  //for (int i = 0; i < ms.m_grid_regularity.size(); i++)
	  //ms.m_grid_regularity[i] = true;
	  
		
	  //std::cout << "initializing the grid..." << std::endl;
	  //initialize the grid point and edge info
	  ms.initialize_grid();
		
	  //std::cout << "classifying the grid points..." << std::endl;
	  //check whether it is inside the atoms
	  for (int i = 0; i < ms.m_atoms_vec.size(); i++)
	    ms.process_each_atom_newer(ms.m_atoms_vec[i]);
	  //cout<<"processing convex sphere complete ... "<<endl;
	  //check whether it is inside the concave sphere and whether it is insie the dual Tet
	  for (int i = 0; i < ms.m_concave_spheres.size(); i++)
	    ms.process_each_concave_sphere_newer(ms.m_concave_spheres[i]);
	  //cout << "processing concave sphere complete ... " << endl;
	  ms.check_regular();
		
	  //check whether it is inside the visibility sphere and inside the tori/saddle faces
	  //if the point is already inside the atoms or inside the concave spheres,leave the tag unchanged
	  for (int i = 0; i < ms.m_torus_vec.size(); i++)
	    ms.process_each_torus_newer(ms.m_torus_vec[i]);
	  //cout << "processing saddle surface complete ... " << endl;

	  ms.check_the_grid_points();
	  ms.process_others_ray_counting();
		
	  //std::cout << "compute intersections..." << std::endl;
		
	  for (int i = 0; i < ms.m_grid_edges.size(); i++)
	    missing_intersection += ms.process_grid_intersection_new(ms.m_grid_edges[i]);

#pragma omp critical
	  {		
	  if (missing_intersection == 0)
	    std::cout << std::endl << "No missing cases!" << std::endl;
	  else
	    {
	      std::cout << std::endl << missing_intersection << " missing cases!!!!!" << std::endl;
	      std::cout<<"Press Enter to continue..."<<std::endl;
	      cin.ignore();
	    }
  }	
	  cout << "mapping grid index..." << endl;
	  for (int l = 0; l < ms.m_z_num; l++)
	    for (int m = 0; m < ms.m_y_num; m++)
	      for (int n = 0; n < ms.m_x_num; n++)
		{
		  //map index of grid_info
		  int x_pos = n + a*BLOCK_SIZE;
		  int y_pos = m + b*BLOCK_SIZE;
		  int z_pos = l + c*BLOCK_SIZE;
			
		  grid_pointIter grid_point = ms.m_grid_points[n + ms.m_x_num*m + ms.m_x_num*ms.m_y_num*l];
		
#pragma omp critical
		  {
		  //ms.m_grid_status[n + m*ms.m_x_num + l*ms.m_x_num*ms.m_y_num] = grid_point->is_inside;
		  t_grid_status[x_pos + y_pos*bi.t_x_num + z_pos*bi.t_x_num*bi.t_y_num] = grid_point->is_inside;	  
		  if(grid_point->is_irregular==true)
		      t_grid_regularity[x_pos + y_pos*bi.t_x_num + z_pos*bi.t_x_num*bi.t_y_num]=false;	
		  }
		  grid_point->grid_x += a*BLOCK_SIZE;
		  grid_point->grid_y += b*BLOCK_SIZE;
		  grid_point->grid_z += c*BLOCK_SIZE;
		}
	  
#pragma omp critical
	  {
    //cout<<"output intersection..."<<endl;
    //int block_size=BLOCK_SIZE;
	  ms.output_intersection(first_write, a, b, c, BLOCK_SIZE, index);
  }	
          /*int for_amino = 0;
	    if (for_amino == 1)
	    read_pqr_and_classify_inte_points("1ajj.pqr");*/
	  //cout<<"clean memory..."<<endl;

    // compute partition area
    ms.partition_area(a, b, c, BLOCK_SIZE);

#pragma omp critical
	  {
      //cout<<"Exporting partition area!!!"<<endl;
     // ms.export_partition_area(t_partition_area);
      
      for(int i=0; i<ms.m_partition_area.size(); ++i){
		    //cout<<i<<" "<<ms.m_local_to_global_atom_idx[i]<<" "<<ms.m_partition_area[i]<<endl;
		    t_partition_area[ms.m_local_to_global_atom_idx[i]] += ms.m_partition_area[i];
      }

      //cout<<"Exporting done!!!"<<endl;
    }
	  ms.clean_memory();
	  
	  
	  cout << "End block (" << a << ", " << b << ", " << c << ")" << endl;
	  //cout << "--------------------------------------------------------" << endl;
	  cout << endl;
	  cout << endl;
	  
	  //ms.compute_volume();
#pragma omp critical
	  {
	    total_surface_area+=ms.m_surface_area;
      cout<<"Area: "<<ms.m_surface_area <<" "<<ms.m_area_from_partition<<" "<<ms.m_another_area_from_partiion<<endl;
	  }
	  //total_volume+=ms.m_surface_volume;
	}
	
  
  compute_volume(total_volume, t_grid_status, t_grid_regularity, bi.m_volume_unit);

  std::cout << std::endl << "surface area: " << total_surface_area << std::endl << "surface volume: " << total_volume << std::endl;
  //numAtoms = bi.m_N_atoms;
  //area = total_surface_area;
  //volume = total_volume;
  //gx = bi.t_x_num;
  //gy = bi.t_y_num;
  //gz = bi.t_z_num;
  //bx = block_dim_x;
  //by = block_dim_y;
  //bz = block_dim_z;

  std::cout << "output info into the files..." << std::endl;
  output_info(t_grid_status, t_partition_area, bi);
  //std::cout<<"clean the memory..."<<std::endl;
}



int main(int argc, char *argv[]){
	auto start = std::chrono::system_clock::now();
 
  
  //molecular_surface mesh;
	int atomNum;
	int gx, gy, gz, bx, by, bz;
	double area, volume;

	process_molecular_surface_newer(atof(argv[2]), atof(argv[3]), atof(argv[4]), argv[1]/*, atomNum, area, volume, gx, gy, gz, bx, by, bz*/);
	//process_molecular_surface_newer(1.4, 0.9, 2.0, "5vkq.xyzr");
   auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  cout << "Running time: " << elapsed_seconds.count() << endl;

  //std::string fn = argv[1];
  //std::stringstream ss(fn);
  //std::string token;
  //std::getline(ss, token, '.');
  /*std::fstream stat("stat.txt", std::ios::app);
  std::string str(argv[5]);
  stat << str << " "
	  << atomNum << " "
	  << gx << "x" << gy << "x" << gz << " "
	  << bx << "x" << by << "x" << bz << " "
	  << area << " "
	  << volume << " "
	  << elapsed_seconds.count() << std::endl;*/

  return 0;
  
}

/*************
output file format11

pos for interior grid, outside grid point
pos for intersection point
normal 
*********************/
