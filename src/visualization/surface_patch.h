// Visialization of convex, saddle and concave patches

#ifndef _SURFACE_PATCH_H_
#define _SURFACE_PATCH_H_

#include "../utility/Types.h"
#include "../intersection/parallel_wrapper.h"
#include "../marchingcubes/MarchingCubes.h"

class surface_patch{
private:
    parallel_wrapper* m_pw;
    MarchingCubes* m_mc;

    std::vector<PatchType> m_vertex_patch_type;
    std::vector<PatchType> m_face_patch_type;

    std::vector<Vertex> m_convex_vertices;
    std::vector<Triangle> m_convex_faces;
    std::vector<Vertex> m_saddle_vertices;
    std::vector<Triangle> m_saddle_faces;
    std::vector<Vertex> m_concave_vertices;
    std::vector<Triangle> m_concave_faces;

public:
    surface_patch();
    surface_patch(parallel_wrapper* pw, MarchingCubes* mc);
    ~surface_patch();

    void compute();
    void write();

    void populate_vertex_patch_type();
    void determine_face_patch_type();
    
    void build_patch(PatchType type);
    void write_patch_convex();
    void write_patch_saddle();
    void write_patch_concave();

    void check_mc_verts_and_pw_verts();

};


#endif