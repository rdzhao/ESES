#include "parallel_wrapper.h"

parallel_wrapper::parallel_wrapper(){

}

parallel_wrapper::parallel_wrapper(std::vector<CVector3d> atoms, std::vector<double> radii,
    double probe_radius, double grid_spacing, double bounding_residual){
        m_basic_info.m_probe_radius = probe_radius;
        m_basic_info.m_extend_bounding = bounding_residual;
        m_basic_info.m_grid_spacing = grid_spacing;
        m_basic_info.m_area_unit = pow(grid_spacing, 2);
        m_basic_info.m_volume_unit = pow(grid_spacing, 3);

        for(int i=0; i<atoms.size(); ++i){
            CVector3d& p = atoms[i];
            double& r = radii[i];
            Atom *new_atom = new Atom(p,r,i);
            m_basic_info.m_atoms_vec.push_back(new_atom);

            if(i==0){
	            m_basic_info.m_mesh_min = CVector3d(p.x()-r, p.y()-r, p.z()-r);
	            m_basic_info.m_mesh_max = CVector3d(p.x()+r, p.y()+r, p.z()+r);
		    }
            else{
                m_basic_info.m_mesh_min[0] = min(m_basic_info.m_mesh_min[0], p.x()-r);
                m_basic_info.m_mesh_min[1] = min(m_basic_info.m_mesh_min[1], p.y()-r);
                m_basic_info.m_mesh_min[2] = min(m_basic_info.m_mesh_min[2], p.z()-r);
            
                m_basic_info.m_mesh_max[0] = max(m_basic_info.m_mesh_max[0], p.x()+r);
                m_basic_info.m_mesh_max[1] = max(m_basic_info.m_mesh_max[1], p.y()+r);
                m_basic_info.m_mesh_max[2] = max(m_basic_info.m_mesh_max[2], p.z()+r);
            }
        }

        m_basic_info.m_N_atoms = m_basic_info.m_atoms_vec.size();
        
        m_basic_info.t_min_x = m_basic_info.m_mesh_min[0] - m_basic_info.m_extend_bounding; 
        m_basic_info.t_min_y = m_basic_info.m_mesh_min[1] - m_basic_info.m_extend_bounding; 
        m_basic_info.t_min_z = m_basic_info.m_mesh_min[2] - m_basic_info.m_extend_bounding;
        m_basic_info.t_max_x = m_basic_info.m_mesh_max[0] + m_basic_info.m_extend_bounding; 
        m_basic_info.t_max_y = m_basic_info.m_mesh_max[1] + m_basic_info.m_extend_bounding; 
        m_basic_info.t_max_z = m_basic_info.m_mesh_max[2] + m_basic_info.m_extend_bounding;

        m_basic_info.t_x_num = int((m_basic_info.t_max_x - m_basic_info.t_min_x)/(m_basic_info.m_grid_spacing)) + 1;
        m_basic_info.t_y_num = int((m_basic_info.t_max_y - m_basic_info.t_min_y)/(m_basic_info.m_grid_spacing)) + 1;
        m_basic_info.t_z_num = int((m_basic_info.t_max_z - m_basic_info.t_min_z)/(m_basic_info.m_grid_spacing)) + 1;
    
        m_basic_info.t_max_x = m_basic_info.t_min_x + (m_basic_info.t_x_num - 1)*m_basic_info.m_grid_spacing;
        m_basic_info.t_max_y = m_basic_info.t_min_y + (m_basic_info.t_y_num - 1)*m_basic_info.m_grid_spacing;
        m_basic_info.t_max_z = m_basic_info.t_min_z + (m_basic_info.t_z_num - 1)*m_basic_info.m_grid_spacing;

        std::cout<<"grid dimension: "
                    <<m_basic_info.t_x_num<<' '
                    <<m_basic_info.t_y_num<<' '
                    <<m_basic_info.t_z_num<<std::endl;
    }

parallel_wrapper::~parallel_wrapper(){

}

basic_information& parallel_wrapper::get_basic_info(){
    return m_basic_info;
}

intersection_info& parallel_wrapper::get_intersection_info(){
    return m_intersection;
}

std::vector<bool>& parallel_wrapper::get_grid_status(){
    return m_grid_status;
}

void parallel_wrapper::process_molecular_surface(){
    std::cout<<"Initialize ..."<<std::endl;
    initialize();
    std::cout<<"Parallel compute ..."<<std::endl;
    parallel_compute();
    std::cout<<"Write ..."<<std::endl;
    write();
}

void parallel_wrapper::initialize(){
    m_block_dim_x = (m_basic_info.t_x_num - 2) / BLOCK_SIZE+1;
    m_block_dim_y = (m_basic_info.t_y_num - 2) / BLOCK_SIZE+1;
    m_block_dim_z = (m_basic_info.t_z_num - 2) / BLOCK_SIZE+1;
    cout << "Block info: " << m_block_dim_x << " " << m_block_dim_y << " " << m_block_dim_z << endl;

    int size = m_basic_info.t_x_num*m_basic_info.t_y_num*m_basic_info.t_z_num;
    m_grid_status.resize(size);
    m_grid_regularity.resize(size, true);
    m_partition_area.resize(m_basic_info.m_N_atoms, 0);
    m_surface_area = 0;
    m_volume = 0;
}

void parallel_wrapper::parallel_compute(){
    int block_size = m_block_dim_x*m_block_dim_y*m_block_dim_z;
    omp_set_num_threads(12);

    #pragma omp parallel for schedule(dynamic)
    for(int index = 0; index<block_size; ++index){
        int block_x, block_y, block_z;
        block_z=index/(m_block_dim_x*m_block_dim_y);
        block_y=(index-block_z*m_block_dim_x*m_block_dim_y)/m_block_dim_x;
        block_x=index-block_z*m_block_dim_x*m_block_dim_y-block_y*m_block_dim_x;

        #pragma omp critical
        {
            cout << "Start block ("<<block_x<<", " << block_y << ", " << block_z << ")" << endl;
        }

        molecular_surface ms;
        prepare_molecular_surface(ms, m_basic_info, block_x, block_y, block_z);

        ms.initialize_molecular_surface();
        int missing_intersection = 0;
        ms.initialize_grid();

        for (int i = 0; i < ms.m_atoms_vec.size(); i++){
	        ms.process_each_atom_newer(ms.m_atoms_vec[i]);
        }
        for (int i = 0; i < ms.m_concave_spheres.size(); i++){
	        ms.process_each_concave_sphere_newer(ms.m_concave_spheres[i]);
        }
        ms.check_regular();
        for (int i = 0; i < ms.m_torus_vec.size(); i++){
	        ms.process_each_torus_newer(ms.m_torus_vec[i]);
        }

        ms.check_the_grid_points();
	    ms.process_others_ray_counting();

        for (int i = 0; i < ms.m_grid_edges.size(); i++){
	        missing_intersection += ms.process_grid_intersection_new(ms.m_grid_edges[i]);
        }

        #pragma omp critical 
        {
	        if (missing_intersection == 0)
	            std::cout << std::endl << "No missing cases!" << std::endl;
	        else{
	            std::cout << std::endl << missing_intersection << " missing cases!!!!!" << std::endl;
	            std::cout<<"Press Enter to continue..."<<std::endl;
	            cin.ignore();
	        }
        }

        #pragma omp critical
		{
            process_grid_status_regularity(ms, block_x, block_y, block_z);
            ms.write_intersections(m_intersection, block_x, block_y, block_z, BLOCK_SIZE);
        }

        ms.partition_area(block_x, block_y, block_z, BLOCK_SIZE);
        #pragma omp critical
	    {
            for(int i=0; i<ms.m_partition_area.size(); ++i){
		        m_partition_area[ms.m_local_to_global_atom_idx[i]] += ms.m_partition_area[i];
            }
        }

        ms.clean_memory();

        #pragma omp critical
	    {
        	cout << "End block (" << block_x << ", " << block_y << ", " << block_z << ")" << endl;
	        //cout << "--------------------------------------------------------" << endl;
        }
    }

    compute_volume();
}

void parallel_wrapper::write(){
    std::cout<<"Write bounding box ..."<<std::endl;
    write_bounding_box();
    std::cout<<"Write grid info ..."<<std::endl;
    write_grid_info();
    std::cout<<"Write intersection info ..."<<std::endl;
    write_intersection_info();
    std::cout<<"Write partition area ..."<<std::endl;
    write_partition_area();
}

void parallel_wrapper::prepare_molecular_surface(molecular_surface& ms, basic_information& bi, int block_x, int block_y, int block_z){
    //critical: find the nearest cluster
    double extension = 2+2*bi.m_probe_radius;
    
    double x_max, y_max, z_max;
    double x_min, y_min, z_min;

    ms.m_x_step=ms.m_y_step=ms.m_z_step=bi.m_grid_spacing;
    
    ms.m_min_x=bi.t_min_x+block_x*BLOCK_SIZE*bi.m_grid_spacing;
    ms.m_min_y=bi.t_min_y+block_y*BLOCK_SIZE*bi.m_grid_spacing;
    ms.m_min_z=bi.t_min_z+block_z*BLOCK_SIZE*bi.m_grid_spacing;

    if (ms.m_min_x + BLOCK_SIZE*ms.m_x_step > bi.t_max_x){
        ms.m_x_num = int((bi.t_max_x - ms.m_min_x) / ms.m_x_step) + 1;
        ms.m_max_x = ms.m_min_x + (ms.m_x_num - 1)*ms.m_x_step;
    }
    else{
        ms.m_max_x = ms.m_min_x + BLOCK_SIZE*ms.m_x_step;
        ms.m_x_num = BLOCK_SIZE + 1;
    }

    if (ms.m_min_y + BLOCK_SIZE*ms.m_y_step > bi.t_max_y){
        ms.m_y_num = int((bi.t_max_y - ms.m_min_y) / ms.m_y_step) + 1;
        ms.m_max_y = ms.m_min_y + (ms.m_y_num - 1)*ms.m_y_step;
    }
    else{
        ms.m_max_y = ms.m_min_y + BLOCK_SIZE*ms.m_y_step;
        ms.m_y_num = BLOCK_SIZE + 1;
    }
  
    if(ms.m_min_z + BLOCK_SIZE*ms.m_z_step > bi.t_max_z){
        ms.m_z_num = int((bi.t_max_z - ms.m_min_z) / ms.m_z_step) + 1;
        ms.m_max_z = ms.m_min_z + (ms.m_z_num - 1)*ms.m_z_step;
    }
    else{
        ms.m_max_z = ms.m_min_z + BLOCK_SIZE*ms.m_z_step;
        ms.m_z_num = BLOCK_SIZE + 1;
    }

    x_min=ms.m_min_x-extension;
    y_min=ms.m_min_y-extension;
    z_min=ms.m_min_z-extension;
    x_max=ms.m_max_x+extension;
    y_max=ms.m_max_y+extension;
    z_max=ms.m_max_z+extension;

    cout<<x_min<<" "<<y_min<<" "<<z_min<<endl;
    cout<<x_max<<" "<<y_max<<" "<<z_max<<endl;

    int idx=0;
    for(int i=0; i<bi.m_N_atoms; i++){
        if((bi.m_atoms_vec[i]->center.x()>x_min && bi.m_atoms_vec[i]->center.x()<x_max) 
            && (bi.m_atoms_vec[i]->center.y()>y_min && bi.m_atoms_vec[i]->center.y()<y_max) 
            && (bi.m_atoms_vec[i]->center.z()>z_min && bi.m_atoms_vec[i]->center.z()<z_max)){
	            Atom *atom_ptr=new Atom(bi.m_atoms_vec[i]->center,bi.m_atoms_vec[i]->r,idx);
	            ms.m_atoms_vec.push_back(atom_ptr);
	            ms.m_local_to_global_atom_idx[idx] = i;
	            idx++;
	        }
    }
    cout<<"Num: "<<idx<<endl;
    ms.m_probe_radius=bi.m_probe_radius;
    ms.m_N_atoms=ms.m_atoms_vec.size();
    ms.m_grid_resolution=bi.m_grid_spacing;
    ms.m_surface_area=ms.m_surface_volume=0;
    ms.m_area_unit=bi.m_area_unit;
    ms.m_volume_unit=bi.m_volume_unit;
}

void parallel_wrapper::process_grid_status_regularity(molecular_surface& ms, int block_x, int block_y, int block_z){
    for (int l = 0; l < ms.m_z_num; l++){
	    for (int m = 0; m < ms.m_y_num; m++){
	        for (int n = 0; n < ms.m_x_num; n++){
		        //map index of grid_info
		        int x_pos = n + block_x*BLOCK_SIZE;
		        int y_pos = m + block_y*BLOCK_SIZE;
		        int z_pos = l + block_z*BLOCK_SIZE;
			
		        grid_pointIter grid_point = ms.m_grid_points[n + ms.m_x_num*m + ms.m_x_num*ms.m_y_num*l];
		
		                //ms.m_grid_status[n + m*ms.m_x_num + l*ms.m_x_num*ms.m_y_num] = grid_point->is_inside;
		        m_grid_status[x_pos + y_pos*m_basic_info.t_x_num + z_pos*m_basic_info.t_x_num*m_basic_info.t_y_num] = grid_point->is_inside;	  
		        if(grid_point->is_irregular==true)
		            m_grid_regularity[x_pos + y_pos*m_basic_info.t_x_num + z_pos*m_basic_info.t_x_num*m_basic_info.t_y_num]=false;	

		        grid_point->grid_x += block_x*BLOCK_SIZE;
		        grid_point->grid_y += block_y*BLOCK_SIZE;
		        grid_point->grid_z += block_z*BLOCK_SIZE;
		    }
        }
    }
}

void parallel_wrapper::compute_volume(){
    for(int i=0; i<m_grid_status.size(); i++){
        if(m_grid_status[i]==true && m_grid_regularity[i]==true)
	        m_volume+=m_basic_info.m_volume_unit;
        else if(m_grid_regularity[i]==false)
	        m_volume+=m_basic_info.m_volume_unit/2;
    }
}

void parallel_wrapper::write_bounding_box(){
    std::ofstream file_bounding("bounding_box.txt");
	file_bounding<<std::scientific;
	file_bounding<<std::setprecision(12)<<m_basic_info.t_min_x<<' '<<m_basic_info.t_min_y<<' '<<m_basic_info.t_min_z<<std::endl;
	file_bounding<<std::setprecision(12)<<m_basic_info.t_max_x<<' '<<m_basic_info.t_max_y<<' '<<m_basic_info.t_max_z<<std::endl;
	file_bounding<<m_basic_info.t_x_num<<' '<<m_basic_info.t_y_num<<' '<<m_basic_info.t_z_num<<std::endl;
	file_bounding.close();
}

void parallel_wrapper::write_grid_info(){
    std::ofstream file_grid("grid_info.txt");
	for (int i = 0; i < m_grid_status.size(); i++){
		int z_pos = int(i / (m_basic_info.t_x_num*m_basic_info.t_y_num));
		int y_pos = int((i - m_basic_info.t_x_num*m_basic_info.t_y_num*z_pos) / m_basic_info.t_x_num);
		int x_pos = i - m_basic_info.t_x_num*m_basic_info.t_y_num*z_pos - m_basic_info.t_x_num*y_pos;
		file_grid << x_pos<<' '<< y_pos << ' ' << z_pos << ' ';
		if(m_grid_status[i]==true)
		    file_grid << 1 << std::endl;
		else
		    file_grid << -1 << std::endl;
	}
	file_grid.close();
}

void parallel_wrapper::write_intersection_info(){
    std::ofstream file_intersection("intersection_info.txt");
    int size = m_intersection.m_intersection.size();
    for(int i=0; i<size; ++i){
        CVector3d& gv_in = m_intersection.m_grid_inside[i];
        CVector3d& gv_out = m_intersection.m_grid_outside[i];
        file_intersection
            <<static_cast<int>(gv_in.x())<<" "
            <<static_cast<int>(gv_in.y())<<" "
            <<static_cast<int>(gv_in.z())<<" ";
        file_intersection
            <<static_cast<int>(gv_out.x())<<" "
            <<static_cast<int>(gv_out.y())<<" "
            <<static_cast<int>(gv_out.z())<<" ";

        CVector3d& pos = m_intersection.m_intersection[i];
        CVector3d& normal = m_intersection.m_normal[i];
        file_intersection<<std::scientific;
        file_intersection<<std::setprecision(12)<<pos.x()<<" "<<pos.y()<<" "<<pos.z()<<" ";
        file_intersection<<std::setprecision(12)<<normal.x()<<" "<<normal.y()<<" "<<normal.z()<<" ";
        
        for(int k=0; k<m_intersection.m_adj_atoms[i].size(); ++k){
            file_intersection<<m_intersection.m_adj_atoms[i][k]<<" ";
        }
        file_intersection<<std::endl;
    }
}

void parallel_wrapper::write_partition_area(){
    std::ofstream file_partition("partition_area.txt");
	for(int i=0; i<m_basic_info.m_N_atoms; i++){
	    file_partition<<i<<" ";
	    file_partition<<std::scientific;
	    file_partition<<std::setprecision(12)<<m_partition_area[i]<<endl;
	}
	file_partition.close();
}