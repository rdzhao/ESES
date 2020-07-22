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

std::vector<int>& pqr_parser::getChain(){
    return m_chain;
}

void pqr_parser::parse(){
    std::ifstream file(m_filename.c_str());

    std::string line;
    while(getline(file, line)){
        std::stringstream ss(line);
        
        //std::string marker;
        //ss>>marker;
        //if(marker != "ATOM") continue;

        //std::string atom_id, atom_type, res_name, chain_id, res_id;
        //std::string x, y, z, charge, radius;

        /*
        ss>>atom_id>>atom_type>>res_name>>chain_id>>res_id;
        ss>>x>>y>>z>>charge>>radius;
        std::cout<<atom_id<<" "<<x<<" "<<y<<" "<<z<<" "<<atom_type[0]<<" "<<chain_id[0]<<" "<<radius<<std::endl;

        m_atoms.push_back(CVector3d(std::stod(x), std::stod(y), std::stod(z)));
        m_types.push_back(atmap.at(atom_type[0]));
        m_chain.push_back(static_cast<int>(chain_id[0]));
        m_radii.push_back(std::stod(radius));
        */

        
        std::string marker = line.substr(0, 6);
        marker.erase(std::remove(marker.begin(), marker.end(), ' '), marker.end());
        if(marker != "ATOM") continue;

        std::string atom_id = line.substr(6,5);

        std::string xs = line.substr(30, 8);
        std::string ys = line.substr(38, 8);
        std::string zs = line.substr(46, 8);
        xs.erase(std::remove(xs.begin(), xs.end(), ' '), xs.end());
        ys.erase(std::remove(ys.begin(), ys.end(), ' '), ys.end());
        zs.erase(std::remove(zs.begin(), zs.end(), ' '), zs.end());
        m_atoms.push_back(CVector3d(std::stod(xs), std::stod(ys), std::stod(zs)));

        std::string ts = line.substr(11, 5);
        ts.erase(std::remove(ts.begin(), ts.end(), ' '), ts.end());
        m_types.push_back(atmap.at(ts[0]));

        std::string cs = line.substr(21,1);
        cs.erase(std::remove(cs.begin(), cs.end(), ' '), cs.end());
        m_chain.push_back(static_cast<int>(cs[0]));

        std::string tr = line.substr(62, 7);
        tr.erase(std::remove(tr.begin(), tr.end(), ' '), tr.end());
        m_radii.push_back(std::stod(tr));

        //std::cout<<atom_id<<" "<<xs<<" "<<ys<<" "<<zs<<" "<<ts[0]<<" "<<cs<<" "<<tr<<std::endl;
    }

    std::cout<<"Atom size: "<<m_atoms.size()<<std::endl;
}