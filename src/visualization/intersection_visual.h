#ifndef _INTERSECTION_VISUAL_H_
#define _INTERSECTION_VISUAL_H_

#include "../utility/Types.h"
#include "../intersection/parallel_wrapper.h"

class intersection_visual{
private:
    parallel_wrapper* m_pw;

    std::set<SVector3I, SVector3ICompare> m_point_set;
    std::set<SVector2I, SVector2ICompare> m_edge_set;
    std::vector<SVector3I> m_points;
    std::map<SVector3I, int, SVector3ICompare> m_index;
    std::vector<SVector2I> m_edges;

public:
    intersection_visual();
    intersection_visual(parallel_wrapper* pw);
    ~intersection_visual();

    void compute();
    void write();

    void assemble_edges();
    void write_intersections();
    void write_edges();

    void put_edge(SVector3I p1, SVector3I p2);

    CVector3d get_position(SVector3I p);
};

#endif