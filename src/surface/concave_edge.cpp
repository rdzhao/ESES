#include "concave_edge.h"


concave_edge::concave_edge(){

}
concave_edge::concave_edge(CVector3d p1,CVector3d p2,int n1,int n2){
	pos1 = p1;
	pos2 = p2;
	atom1 = n1;
	atom2 = n2;
}
void concave_edge::initial_concave_edge(CVector3d probe_center,CVector3d torus_center,CVector3d torus_axis){
	n_ijk = (probe_center - torus_center).Cross(torus_axis);
	n_ijk.Normalize();
}