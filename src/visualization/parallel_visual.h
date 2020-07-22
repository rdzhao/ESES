#ifndef _PARALLEL_VISUAL_H_
#define _PARALLEL_VISUAL_H_

#include "../utility/Types.h"
#include "../intersection/parallel_wrapper.h"
#include "../marchingcubes/MarchingCubes.h"

class parallel_visual{
private:
    parallel_wrapper* m_pw;
    MarchingCubes* m_mc;

    vector<int> m_facet_block_id;

    std::vector<std::vector<int>> m_block_vertices;
    std::vector<std::vector<Triangle>> m_block_facets;

public:
    parallel_visual();
    parallel_visual(parallel_wrapper* pw, MarchingCubes* mc);
    ~parallel_visual();

    void compute();
    void write();

    void determine_facet_block_id();
    void build_surface_block_id();
    void write_surface();

};

#endif