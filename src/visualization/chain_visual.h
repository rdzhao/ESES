#ifndef _CHAIN_VISUAL_H_
#define _CHAIN_VISUAL_H_

#include "../utility/Types.h"
#include "../parser/pqr_parser.h"
#include "../intersection/parallel_wrapper.h"
#include "../marchingcubes/MarchingCubes.h"

class chain_visual{
private:
    pqr_parser* m_pp;
    parallel_wrapper* m_pw;
    MarchingCubes* m_mc;

    std::vector<int> m_vertex_chain_id;
    std::vector<int> m_facet_chain_id;

    std::vector<std::vector<int>> m_chain_vertices;
    std::vector<std::vector<Triangle>> m_chain_facets;

public:
    chain_visual();
    chain_visual(pqr_parser* pp, parallel_wrapper* pw, MarchingCubes* mc);
    ~chain_visual();

    void compute();
    void write();
    
    void populate_vertex_chain_id();
    void determine_triangle_chain_id();
    void build_surfaces_chain_id();
    void write_surfaces();
};

#endif