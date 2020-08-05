#ifndef _PQR_PARSER_H_
#define _PQR_PARSER_H_

#include "../utility/Types.h"
#include "../utility/Vector3d.h"

class pqr_parser{
private:
    std::string m_filename;
    std::vector<CVector3d> m_atoms;
    std::vector<double> m_radii;
    std::vector<AtomType> m_types;
    std::vector<int> m_chain; 

public:
    pqr_parser();
    pqr_parser(std::string fn);
    ~pqr_parser();

    std::vector<CVector3d>& getAtoms();
    std::vector<double>& getRadii();
    std::vector<AtomType>& getType();
    std::vector<int>& getChain();

    void parse();

    void write_point_cloud();
};

#endif