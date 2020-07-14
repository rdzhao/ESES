#ifndef _PARALLEL_WRAPPER_H_
#define _PARALLEL_WRAPPER_H_

#include "molecular_surface.h"

const int BLOCK_SIZE = 63;

class basic_information{
public:
    double m_probe_radius;
    double m_extend_bounding;
    double m_grid_spacing;
    double m_area_unit;
    double m_volume_unit;
    int m_N_atoms;

    std::vector<Atom *> m_atoms_vec;
    CVector3d m_mesh_min,m_mesh_max;
  
    double t_min_x, t_min_y, t_min_z;  //total bounding box
    double t_max_x, t_max_y, t_max_z;
    int t_x_num, t_y_num, t_z_num; //total number of grid points in each dimension
};

class parallel_wrapper{
private:
    basic_information m_basic_info;

    int m_block_dim_x;
    int m_block_dim_y;
    int m_block_dim_z;

    std::vector<bool> m_grid_status;
    std::vector<bool> m_grid_regularity;
    
    intersection_info m_intersection;

    std::vector<double> m_partition_area;
    double m_surface_area;
    double m_volume;

public:
    parallel_wrapper();
    parallel_wrapper(std::vector<CVector3d> atoms, std::vector<double> radii,
                        double probe_size, double grid_spacing, double bounding_residual);
    ~parallel_wrapper();

    basic_information& get_basic_info();
    intersection_info& get_intersection_info();
    std::vector<bool>& get_grid_status();

    void process_molecular_surface();

    void initialize();
    void parallel_compute();
    void write();

    void prepare_molecular_surface(molecular_surface& ms, basic_information& bi, int block_x, int block_y, int block_z);
    void process_grid_status_regularity(molecular_surface& ms, int block_x, int block_y, int block_z);
    void compute_volume();

    void write_bounding_box();
    void write_grid_info();
    void write_intersection_info();
    void write_partition_area();
};

#endif

