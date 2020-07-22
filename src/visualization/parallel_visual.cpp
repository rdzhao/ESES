#include "parallel_visual.h"

parallel_visual::parallel_visual(){

}

parallel_visual::parallel_visual(parallel_wrapper* pw, MarchingCubes* mc)
    : m_pw(pw), m_mc(mc){

}

parallel_visual::~parallel_visual(){

}

void parallel_visual::compute(){
    determine_facet_block_id();
    build_surface_block_id();
}

void parallel_visual::write(){
    write_surface();
}

void parallel_visual::determine_facet_block_id(){
    double block_length = m_pw->get_basic_info().m_grid_spacing*BLOCK_SIZE;
    int block_dim_x = m_pw->get_block_dim_x();
    int block_dim_y = m_pw->get_block_dim_y();
    int block_dim_z = m_pw->get_block_dim_z();

    auto& vb = m_mc->vertex_buffer();
    auto& tb = m_mc->triangle_buffer();
    m_facet_block_id.reserve(tb.size());
    for(int i=0; i<tb.size(); ++i){
        CVector3d v1(vb[tb[i].v1].x, vb[tb[i].v1].y, vb[tb[i].v1].z);
        CVector3d v2(vb[tb[i].v2].x, vb[tb[i].v2].y, vb[tb[i].v2].z);
        CVector3d v3(vb[tb[i].v3].x, vb[tb[i].v3].y, vb[tb[i].v3].z);
        CVector3d ave = (v1+v2+v3)/3;

        int a = (ave.x() - m_pw->get_basic_info().t_min_x) / block_length;
        int b = (ave.y() - m_pw->get_basic_info().t_min_y) / block_length;
        int c = (ave.z() - m_pw->get_basic_info().t_min_z) / block_length;

        int block_id = a*block_dim_y*block_dim_z + b*block_dim_z + c;
        m_facet_block_id.push_back(block_id);
    }
}

void parallel_visual::build_surface_block_id(){
    std::set<int> block_id_set;
    for(int i=0; i<m_facet_block_id.size(); ++i){
        block_id_set.insert(m_facet_block_id[i]);
    }

    cout<<"Block size: "<<block_id_set.size()<<endl;

    auto& tb = m_mc->triangle_buffer();
    for(auto e:block_id_set){
        std::unordered_map<int, int> new_idx;
        std::vector<int> vertices;
        std::vector<Triangle> facets;

        int idx = 0;
        for(int i=0; i<tb.size(); ++i){
            if(m_facet_block_id[i] != e) continue;
            
            if(new_idx.count(tb[i].v1) == 0){
                new_idx[tb[i].v1] = idx++;
                vertices.push_back(tb[i].v1);
            }
            if(new_idx.count(tb[i].v2) == 0){
                new_idx[tb[i].v2] = idx++;
                vertices.push_back(tb[i].v2);
            } 
            if(new_idx.count(tb[i].v3) == 0){
                new_idx[tb[i].v3] = idx++;
                vertices.push_back(tb[i].v3);
            }

            Triangle t;
            t.v1 = new_idx[tb[i].v1];
            t.v2 = new_idx[tb[i].v3];
            t.v3 = new_idx[tb[i].v2];
            facets.push_back(t);
        }

        m_block_vertices.push_back(vertices);
        m_block_facets.push_back(facets);
    }
}

void parallel_visual::write_surface(){
    for(int i=0; i<m_block_vertices.size(); ++i){
        string filename = "parallel_visual/block_"+std::to_string(i)+".obj";
        ofstream out(filename);

        for(int j=0; j<m_block_vertices[i].size(); ++j){
            int& idx = m_block_vertices[i][j];
            auto& vertex = m_mc->vertex_buffer()[idx];
            out<<"v "<<vertex.x<<" "<<vertex.y<<" "<<vertex.z<<std::endl;
        }

        for(int j=0; j<m_block_facets[i].size(); ++j){
            auto& t = m_block_facets[i][j];
            out<<"f "<<t.v1+1<<" "<<t.v2+1<<" "<<t.v3+1<<std::endl;
        }
    }
}