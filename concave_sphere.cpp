#include <algorithm>
#include "concave_sphere.h"

Concave_Sphere::Concave_Sphere(){

}

Concave_Sphere::Concave_Sphere(concave_edgeIter e1, concave_edgeIter e2,concave_edgeIter e3, int i,int j, int k,CVector3d a_i,CVector3d a_j, CVector3d a_k,double radius_i,
							   double radius_j,double radius_k,CVector3d probe){
			   edge_ij= e1;
			   edge_jk = e2;
			   edge_ki = e3;
			   index_i = i;
			   index_j = j;
			   index_k= k;
			   probe_center = probe;

			   center_i = a_i;
			   center_j = a_j;
			   center_k = a_k;

			   r_i = radius_i;
			   r_j = radius_j;
			   r_k = radius_k;		
			   has_check_duplex = false;
}