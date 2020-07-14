#include "pqr_parser.h"

pqr_parser::pqr_parser(){

}

pqr_parser::pqr_parser(std::string fn) : m_filename(fn){

}

pqr_parser::~pqr_parser(){
    
}

std::vector<CVector3d>& pqr_parser::getAtoms(){
    return m_atoms;
}

std::vector<double>& pqr_parser::getRadii(){
    return m_radii;
}

std::vector<AtomType>& pqr_parser::getType(){
    return m_types;
}

void pqr_parser::parse(){
    std::ifstream file(m_filename.c_str());

    std::string line;
    while(getline(file, line)){
        std::string marker = line.substr(0, 6);
        marker.erase(std::remove(marker.begin(), marker.end(), ' '), marker.end());
        if(marker != "ATOM") continue;

        std::string xs = line.substr(31, 8);
        std::string ys = line.substr(39, 8);
        std::string zs = line.substr(47, 8);
        xs.erase(std::remove(xs.begin(), xs.end(), ' '), xs.end());
        ys.erase(std::remove(ys.begin(), ys.end(), ' '), ys.end());
        zs.erase(std::remove(zs.begin(), zs.end(), ' '), zs.end());
        m_atoms.push_back(CVector3d(std::stod(xs), std::stod(ys), std::stod(zs)));

        std::string ts = line.substr(12, 4);
        ts.erase(std::remove(ts.begin(), ts.end(), ' '), ts.end());
        m_types.push_back(atmap.at(ts[0]));

        std::string tr = line.substr(63, 6);
        tr.erase(std::remove(tr.begin(), tr.end(), ' '), tr.end());
        m_radii.push_back(std::stod(tr));

        //std::cout<<xs<<" "<<ys<<" "<<zs<<" "<<ts[0]<<" "<<tr<<std::endl;
    }

    std::cout<<"Atom size: "<<m_atoms.size()<<std::endl;
}