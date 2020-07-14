#include "space_fill.h"

space_fill::space_fill(){

}
space_fill::space_fill(std::vector<CVector3d> atoms, std::vector<double> radii, std::vector<AtomType> types)
    : m_atoms(atoms), m_radii(radii), m_types(types){

}

space_fill::~space_fill(){

}

void space_fill::write(){
    write_atoms(ATOM_C, "C");
    //write_atoms(ATOM_H, "H");
    write_atoms(ATOM_O, "O");
    write_atoms(ATOM_N, "N");
    write_atoms(ATOM_P, "P");
    write_atoms(ATOM_P, "S");
}

void space_fill::write_all_atoms(){
    
}

void space_fill::write_atoms(AtomType at,std::string sfx){
    std::string fn = "space_fill_"+sfx+".obj";
    std::ofstream out(fn.c_str());

    for(int i=0; i<m_types.size(); ++i){
        if(m_types[i] != at) continue;

        for(int k=0; k<points.size(); k+=3){
            out<<"v "
                <<m_radii[i]*points[k]+m_atoms[i].x()<<" "
                <<m_radii[i]*points[k+1]+m_atoms[i].y()<<" "
                <<m_radii[i]*points[k+2]+m_atoms[i].z()<<std::endl;
        }
    }

    int cnt = 0;
    for(int i=0; i<m_types.size(); ++i){
        if(m_types[i] != at) continue;

        for(int k=0; k<faces.size(); k+=3){
            out<<"f "
                <<faces[k]+cnt+1<<" "
                <<faces[k+1]+cnt+1<<" "
                <<faces[k+2]+cnt+1<<std::endl;
        }
        cnt+=points.size()/3;
    }
}

