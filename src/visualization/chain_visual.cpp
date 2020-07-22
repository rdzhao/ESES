#include "chain_visual.h"

chain_visual::chain_visual(){

}

chain_visual::chain_visual(pqr_parser* pp, parallel_wrapper* pw, MarchingCubes* mc) 
: m_pp(pp), m_pw(pw), m_mc(mc){

}

chain_visual::~chain_visual(){

}

void chain_visual::compute(){
    populate_vertex_chain_id();
    determine_triangle_chain_id();
    build_surfaces_chain_id();
}

void chain_visual::write(){
    write_surfaces();
}
    
void chain_visual::populate_vertex_chain_id(){
    std::vector<CVector3d>& atoms = m_pp->getAtoms();
    std::vector<int>& chain_id = m_pp->getChain();
    intersection_info& ii = m_pw->get_intersection_info();
    m_vertex_chain_id.reserve(ii.m_intersection.size());
    for(int i=0; i<ii.m_intersection.size(); ++i){
        int closest_atom = 0;
        double dist = std::numeric_limits<double>::max();
        for(int j=0; j<ii.m_adj_atoms[i].size(); ++j){
            double d = (ii.m_intersection[i] - atoms[ii.m_adj_atoms[i][j]]).Length();
            if(d<dist){
                closest_atom = ii.m_adj_atoms[i][j];
                dist = d;
            }
        }
        m_vertex_chain_id.push_back(chain_id[closest_atom]);
    }
}

void chain_visual::determine_triangle_chain_id(){
    auto& tb = m_mc->triangle_buffer();
    m_facet_chain_id.reserve(tb.size());
    for(int i=0; i<tb.size(); ++i){
        unordered_map<int, int> cnt;
        cnt[m_vertex_chain_id[tb[i].v1]]++;
        cnt[m_vertex_chain_id[tb[i].v2]]++;
        cnt[m_vertex_chain_id[tb[i].v3]]++;

        bool determined = false;
        for(auto e:cnt){
            if(e.second<2) continue;
            determined = true;
            m_facet_chain_id.push_back(e.first);
        }
        if(!determined) m_facet_chain_id.push_back(m_vertex_chain_id[tb[i].v1]);
    }
}

void chain_visual::build_surfaces_chain_id(){
    std::set<int> chain_id_set;
    for(int i=0; i<m_facet_chain_id.size(); ++i){
        chain_id_set.insert(m_facet_chain_id[i]);
    }

    auto& tb = m_mc->triangle_buffer();
    for(auto e:chain_id_set){
        std::unordered_map<int, int> new_idx;
        std::vector<int> vertices;
        std::vector<Triangle> facets;
        
        int idx = 0;
        for(int i=0; i<tb.size(); ++i){
            if(m_facet_chain_id[i] != e) continue;
            
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

        m_chain_vertices.push_back(vertices);
        m_chain_facets.push_back(facets);
    }
}

void chain_visual::write_surfaces(){
    for(int i=0; i<m_chain_vertices.size(); ++i){
        string filename = "chain_visual/chain_"+std::to_string(i)+".obj";
        ofstream out(filename);

        for(int j=0; j<m_chain_vertices[i].size(); ++j){
            int& idx = m_chain_vertices[i][j];
            auto& vertex = m_mc->vertex_buffer()[idx];
            out<<"v "<<vertex.x<<" "<<vertex.y<<" "<<vertex.z<<std::endl;
        }

        for(int j=0; j<m_chain_facets[i].size(); ++j){
            auto& t = m_chain_facets[i][j];
            out<<"f "<<t.v1+1<<" "<<t.v2+1<<" "<<t.v3+1<<std::endl;
        }
    }
}