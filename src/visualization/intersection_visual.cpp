#include "intersection_visual.h"

intersection_visual::intersection_visual(){

}

intersection_visual::intersection_visual(parallel_wrapper* pw)
    : m_pw(pw){

}

intersection_visual::~intersection_visual(){

}

void intersection_visual::compute(){
    //std::cout<<"Computing ..."<<std::endl;
    assemble_edges();
    //std::cout<<"End ..."<<std::endl;
}

void intersection_visual::write(){
    write_intersections();
    write_edges();
}

void intersection_visual::assemble_edges(){
    auto& info = m_pw->get_intersection_info();
    std::cout<<"Intersection Size: "<<info.m_intersection.size()<<std::endl;
    for(int i=0; i<info.m_intersection.size(); ++i){
        //cout<<"Intersection: "<<i<<std::endl;
        auto& tg_in = info.m_grid_inside[i];
        auto& tg_out = info.m_grid_outside[i];
        SVector3I g_in(static_cast<int>(tg_in.x()), 
                        static_cast<int>(tg_in.y()), 
                        static_cast<int>(tg_in.z()));
        SVector3I g_out(static_cast<int>(tg_out.x()), 
                        static_cast<int>(tg_out.y()), 
                        static_cast<int>(tg_out.z()));


        if(g_in.x != g_out.x){ // along x
            for(int k=-1; k<=1; ++k){
                put_edge(SVector3I(g_in.x, g_in.y, g_in.z+k), SVector3I(g_in.x, g_in.y+1, g_in.z+k));
                put_edge(SVector3I(g_in.x, g_in.y, g_in.z+k), SVector3I(g_in.x, g_in.y-1, g_in.z+k));
                put_edge(SVector3I(g_in.x, g_in.y+k, g_in.z), SVector3I(g_in.x, g_in.y+k, g_in.z+1));
                put_edge(SVector3I(g_in.x, g_in.y+k, g_in.z), SVector3I(g_in.x, g_in.y+k, g_in.z-1));

                put_edge(SVector3I(g_out.x, g_out.y, g_out.z+k), SVector3I(g_out.x, g_out.y+1, g_out.z+k));
                put_edge(SVector3I(g_out.x, g_out.y, g_out.z+k), SVector3I(g_out.x, g_out.y-1, g_out.z+k));
                put_edge(SVector3I(g_out.x, g_out.y+k, g_out.z), SVector3I(g_out.x, g_out.y+k, g_out.z+1));
                put_edge(SVector3I(g_out.x, g_out.y+k, g_out.z), SVector3I(g_out.x, g_out.y+k, g_out.z-1));
            }

            
            for(int k=-1; k<=1; ++k){
                for(int l=-1; l<=1; ++l){
                    put_edge(SVector3I(g_in.x, g_in.y+k, g_in.z+l), SVector3I(g_out.x, g_out.y+k, g_out.z+l));
                }
            }
            
        }
        else if(g_in.y != g_out.y){ // along y
            for(int k=-1; k<=1; ++k){
                put_edge(SVector3I(g_in.x, g_in.y, g_in.z+k), SVector3I(g_in.x+1, g_in.y, g_in.z+k));
                put_edge(SVector3I(g_in.x, g_in.y, g_in.z+k), SVector3I(g_in.x-1, g_in.y, g_in.z+k));
                put_edge(SVector3I(g_in.x+k, g_in.y, g_in.z), SVector3I(g_in.x+k, g_in.y, g_in.z+1));
                put_edge(SVector3I(g_in.x+k, g_in.y, g_in.z), SVector3I(g_in.x+k, g_in.y, g_in.z-1));

                put_edge(SVector3I(g_out.x, g_out.y, g_out.z+k), SVector3I(g_out.x+1, g_out.y, g_out.z+k));
                put_edge(SVector3I(g_out.x, g_out.y, g_out.z+k), SVector3I(g_out.x-1, g_out.y, g_out.z+k));
                put_edge(SVector3I(g_out.x+k, g_out.y, g_out.z), SVector3I(g_out.x+k, g_out.y, g_out.z+1));
                put_edge(SVector3I(g_out.x+k, g_out.y, g_out.z), SVector3I(g_out.x+k, g_out.y, g_out.z-1));
            }

            
            for(int k=-1; k<=1; ++k){
                for(int l=-1; l<=1; ++l){
                    put_edge(SVector3I(g_in.x+k, g_in.y, g_in.z+l), SVector3I(g_out.x+k, g_out.y, g_out.z+l));
                }
            }
            
        }
        else{ // g_in.y != g_out.y, along z
            for(int k=-1; k<=1; ++k){
                put_edge(SVector3I(g_in.x+k, g_in.y, g_in.z), SVector3I(g_in.x+k, g_in.y+1, g_in.z));
                put_edge(SVector3I(g_in.x+k, g_in.y, g_in.z), SVector3I(g_in.x+k, g_in.y-1, g_in.z));
                put_edge(SVector3I(g_in.x, g_in.y+k, g_in.z), SVector3I(g_in.x+1, g_in.y+k, g_in.z));
                put_edge(SVector3I(g_in.x, g_in.y+k, g_in.z), SVector3I(g_in.x-1, g_in.y+k, g_in.z));

                put_edge(SVector3I(g_out.x+k, g_out.y, g_out.z), SVector3I(g_out.x+k, g_out.y+1, g_out.z));
                put_edge(SVector3I(g_out.x+k, g_out.y, g_out.z), SVector3I(g_out.x+k, g_out.y-1, g_out.z));
                put_edge(SVector3I(g_out.x, g_out.y+k, g_out.z), SVector3I(g_out.x+1, g_out.y+k, g_out.z));
                put_edge(SVector3I(g_out.x, g_out.y+k, g_out.z), SVector3I(g_out.x-1, g_out.y+k, g_out.z));
            }
            
            for(int k=-1; k<=1; ++k){
                for(int l=-1; l<=1; ++l){
                    put_edge(SVector3I(g_in.x+l, g_in.y+k, g_in.z), SVector3I(g_out.x+l, g_out.y+k, g_out.z));
                }
            }
            
        }
    }
}

void intersection_visual::write_intersections(){
    std::ofstream out("intersection_visal_points.obj");

    auto& info = m_pw->get_intersection_info();
    for(int i=0; i<info.m_intersection.size(); ++i){
        out<<"v "
            <<info.m_intersection[i].x()<<" "
            <<info.m_intersection[i].y()<<" "
            <<info.m_intersection[i].z()<<std::endl;
    }

    out.close();
}

void intersection_visual::write_edges(){
    std::ofstream out("intersection_visual_edges.obj");

    for(int i=0; i<m_points.size(); ++i){
        CVector3d pos = get_position(m_points[i]);
        out<<"v "<<pos.x()<<" "<<pos.y()<<" "<<pos.z()<<std::endl;
    }

    for(int i=0; i<m_edges.size(); ++i){
        out<<"l "<<m_edges[i].x+1<<" "<<m_edges[i].y+1<<std::endl;
    }
}

void intersection_visual::put_edge(SVector3I p1, SVector3I p2){
    if(m_point_set.find(p1) == m_point_set.end()){
        m_index[p1] = m_points.size();
        m_points.push_back(p1);
        m_point_set.insert(p1);
    }
    if(m_point_set.find(p2) == m_point_set.end()){
        m_index[p2] = m_points.size();
        m_points.push_back(p2);
        m_point_set.insert(p2);
    }

    auto i1 = m_index[p1];
    auto i2 = m_index[p2];
    if(i1>i2) swap(i1, i2);

    SVector2I e(i1, i2);
    if(m_edge_set.find(e) == m_edge_set.end()){
        m_edges.push_back(e);
        m_edge_set.insert(e);
    }
}

CVector3d intersection_visual::get_position(SVector3I p){
    basic_information& info = m_pw->get_basic_info();
    return CVector3d(info.t_min_x + p.x*info.m_grid_spacing, 
                        info.t_min_y + p.y*info.m_grid_spacing,
                        info.t_min_z + p.z*info.m_grid_spacing);
}