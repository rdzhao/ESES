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

#include "parser/pqr_parser.h"
#include "intersection/molecular_surface.h"
#include "intersection/parallel_wrapper.h"
#include "marchingcubes/MarchingCubes.h"
#include "marchingcubes/LevelSet.h"
#include "visualization/space_fill.h"
#include "visualization/surface_patch.h"
#include "visualization/intersection_visual.h"

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

int main(int argc, char *argv[]){
	auto start = std::chrono::system_clock::now();

    // process pqr file
    std::cout<<CONSOLE_YELLOW<<"Parsing PQR file ... "<<CONSOLE_WHITE<<std::endl;
    pqr_parser pp(argv[1]);
    pp.parse();

    // generate catesian mesh
	//process_molecular_surface_newer(atof(argv[2]), atof(argv[3]), atof(argv[4]), pp.getAtoms(), pp.getRadii());
    std::cout<<CONSOLE_YELLOW<<"Parallel computing Eulerian mesh ... "<<CONSOLE_WHITE<<std::endl;
    parallel_wrapper pw(pp.getAtoms(), pp.getRadii(), atof(argv[2]), atof(argv[3]), atof(argv[4]));
    pw.process_molecular_surface();

    std::cout<<CONSOLE_YELLOW<<"Marching cubes: Lagrangian mesh ... "<<CONSOLE_WHITE<<std::endl;
    MarchingCubes mc(&pw);
    mc.process();

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout <<CONSOLE_GREEN<< "Running time: " << elapsed_seconds.count() <<CONSOLE_WHITE<< std::endl;

    std::cout<<CONSOLE_YELLOW<<"Visualization ... "<<CONSOLE_WHITE<<std::endl;
    // visualization: space fill
    space_fill sf(pp.getAtoms(), pp.getRadii(), pp.getType());
    sf.write();

    // visualization: surface patch
    surface_patch sp(&pw, &mc);
    //sp.check_mc_verts_and_pw_verts();
    sp.compute();
    sp.write();

    // visualization: intersection visual
    intersection_visual iv(&pw);
    iv.compute();
    iv.write();

    return 0;
}

/*************
output file format11

pos for interior grid, outside grid point
pos for intersection point
normal 
*********************/
