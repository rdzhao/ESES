#include "surface_patch.h"

surface_patch::surface_patch(){

}

surface_patch::surface_patch(parallel_wrapper* pw, MarchingCubes* mc)
    : m_pw(pw), m_mc(mc){
    
}

surface_patch::~surface_patch(){

}

void surface_patch::compute(){
    populate_vertex_patch_type();
    determine_face_patch_type();

    build_patch(PATCH_CONVEX);
    build_patch(PATCH_SADDLE);
    build_patch(PATCH_CONCAVE);
}

void surface_patch::write(){
    write_patch_convex();
    write_patch_saddle();
    write_patch_concave();
}

void surface_patch::populate_vertex_patch_type(){
    int size = m_pw->get_intersection_info().m_intersection.size();
    m_vertex_patch_type.resize(size);
    for(int i=0; i<m_pw->get_intersection_info().m_intersection.size(); ++i){
        if(m_pw->get_intersection_info().m_adj_atoms[i].size() == 1){
            m_vertex_patch_type[i] = PATCH_CONVEX;
        }
        else if(m_pw->get_intersection_info().m_adj_atoms[i].size() == 2){
            m_vertex_patch_type[i] = PATCH_SADDLE;
        }
        else if(m_pw->get_intersection_info().m_adj_atoms[i].size() == 3){
            m_vertex_patch_type[i] = PATCH_CONCAVE;
        }
    }
}

void surface_patch::determine_face_patch_type(){
    int size = m_mc->triangle_buffer().size();
    m_face_patch_type.resize(size);
    auto& tb = m_mc->triangle_buffer();
    for(int i=0; i<tb.size(); ++i){
        int ncv, nsd, ncc;
        ncv = nsd = ncc = 0;
        
        if(m_vertex_patch_type[tb[i].v1] == PATCH_CONVEX) ncv++;
        else if(m_vertex_patch_type[tb[i].v1] == PATCH_SADDLE) nsd++;
        else if(m_vertex_patch_type[tb[i].v1] == PATCH_CONCAVE) ncc++;

        if(m_vertex_patch_type[tb[i].v2] == PATCH_CONVEX) ncv++;
        else if(m_vertex_patch_type[tb[i].v2] == PATCH_SADDLE) nsd++;
        else if(m_vertex_patch_type[tb[i].v2] == PATCH_CONCAVE) ncc++;

        if(m_vertex_patch_type[tb[i].v3] == PATCH_CONVEX) ncv++;
        else if(m_vertex_patch_type[tb[i].v3] == PATCH_SADDLE) nsd++;
        else if(m_vertex_patch_type[tb[i].v3] == PATCH_CONCAVE) ncc++;

        if(ncv > 1) m_face_patch_type[i] = PATCH_CONVEX;
        else if(nsd > 1) m_face_patch_type[i] = PATCH_SADDLE;
        else if(ncc > 1) m_face_patch_type[i] = PATCH_CONCAVE;
        else m_face_patch_type[i] = PATCH_CONVEX;
    }

    // check
    int a,b,c;
    a=b=c=0;
    for(int i=0; i<m_face_patch_type.size(); ++i){
        if(m_face_patch_type[i] == PATCH_CONVEX) a++;
        else if(m_face_patch_type[i] == PATCH_SADDLE) b++;
        else if(m_face_patch_type[i] == PATCH_CONCAVE) c++;
    }
    std::cout<<"Patch triangle num: "<<a<<" "<<b<<" "<<c<<std::endl;
}

void surface_patch::build_patch(PatchType type){
    std::unordered_map<int, int> new_vertex_idx;
    int idx = 0;
    auto& tb = m_mc->triangle_buffer();
    auto& vb = m_mc->vertex_buffer();
    for(int i=0; i<tb.size(); ++i){
        if(m_face_patch_type[i] != type) continue;

        if(new_vertex_idx.count(tb[i].v1) == 0){
            new_vertex_idx[tb[i].v1] = idx;
            ++idx;
            if(type == PATCH_CONVEX) m_convex_vertices.push_back(vb[tb[i].v1]);
            else if(type == PATCH_SADDLE) m_saddle_vertices.push_back(vb[tb[i].v1]);
            else if(type == PATCH_CONCAVE) m_concave_vertices.push_back(vb[tb[i].v1]);
        }
        if(new_vertex_idx.count(tb[i].v2) == 0){
            new_vertex_idx[tb[i].v2] = idx;
            ++idx;
            if(type == PATCH_CONVEX) m_convex_vertices.push_back(vb[tb[i].v2]);
            else if(type == PATCH_SADDLE) m_saddle_vertices.push_back(vb[tb[i].v2]);
            else if(type == PATCH_CONCAVE) m_concave_vertices.push_back(vb[tb[i].v2]);
        }
        if(new_vertex_idx.count(tb[i].v3) == 0){
            new_vertex_idx[tb[i].v3] = idx;
            ++idx;
            if(type == PATCH_CONVEX) m_convex_vertices.push_back(vb[tb[i].v3]);
            else if(type == PATCH_SADDLE) m_saddle_vertices.push_back(vb[tb[i].v3]);
            else if(type == PATCH_CONCAVE) m_concave_vertices.push_back(vb[tb[i].v3]);
        }

        Triangle t;
        t.v1 = new_vertex_idx[tb[i].v1];
        t.v2 = new_vertex_idx[tb[i].v3];
        t.v3 = new_vertex_idx[tb[i].v2];
        if(type == PATCH_CONVEX) m_convex_faces.push_back(t);
        else if(type == PATCH_SADDLE) m_saddle_faces.push_back(t);
        else if(type == PATCH_CONCAVE) m_concave_faces.push_back(t);
    }
}

void surface_patch::write_patch_convex(){
    std::ofstream out("patch_convex.obj");
    for(int i=0; i<m_convex_vertices.size(); ++i){
        out<<"v "
            <<m_convex_vertices[i].x<<" "
            <<m_convex_vertices[i].y<<" "
            <<m_convex_vertices[i].z<<std::endl;
    }
    for(int i=0; i<m_convex_faces.size(); ++i){
        out<<"f "
            <<m_convex_faces[i].v1+1<<" "
            <<m_convex_faces[i].v2+1<<" "
            <<m_convex_faces[i].v3+1<<std::endl;
    }
    out.close();
}

void surface_patch::write_patch_saddle(){
    std::ofstream out("patch_saddle.obj");
    for(int i=0; i<m_saddle_vertices.size(); ++i){
        out<<"v "
            <<m_saddle_vertices[i].x<<" "
            <<m_saddle_vertices[i].y<<" "
            <<m_saddle_vertices[i].z<<std::endl;
    }
    for(int i=0; i<m_saddle_faces.size(); ++i){
        out<<"f "
            <<m_saddle_faces[i].v1+1<<" "
            <<m_saddle_faces[i].v2+1<<" "
            <<m_saddle_faces[i].v3+1<<std::endl;
    }
    out.close();
}

void surface_patch::write_patch_concave(){
    std::ofstream out("patch_concave.obj");
    for(int i=0; i<m_concave_vertices.size(); ++i){
        out<<"v "
            <<m_concave_vertices[i].x<<" "
            <<m_concave_vertices[i].y<<" "
            <<m_concave_vertices[i].z<<std::endl;
    }
    for(int i=0; i<m_concave_faces.size(); ++i){
        out<<"f "
            <<m_concave_faces[i].v1+1<<" "
            <<m_concave_faces[i].v2+1<<" "
            <<m_concave_faces[i].v3+1<<std::endl;
    }
    out.close();
}

void surface_patch::check_mc_verts_and_pw_verts(){
    if(m_mc->vertex_buffer().size() != m_pw->get_intersection_info().m_intersection.size()){
        std::cout<<"Error: Different number of vertices ..."<<std::endl;
        exit(0);
    }

    for(int i=0; i<m_mc->vertex_buffer().size(); ++i){
        auto& v_mc = m_mc->vertex_buffer()[i];
        auto& v_pw = m_pw->get_intersection_info().m_intersection[i];

        double diff = pow((v_mc.x-v_pw.x()), 2)
                        + pow((v_mc.y-v_pw.y()), 2)
                        + pow((v_mc.z-v_pw.z()), 2);

        if(diff > 0.001) {
            std::cout<<diff<<std::endl;
        }
    }
}