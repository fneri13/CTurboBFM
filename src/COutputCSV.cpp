#include "COutputCSV.hpp"

void COutputCSV::writeSolution(){
    std::string filename = _config.getSolutionName() + ".csv";
    std::ofstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    size_t ni = _mesh.getNumberPointsI();
    size_t nj = _mesh.getNumberPointsJ();
    size_t nk = _mesh.getNumberPointsK();

    file << "NI=" << ni << "\n";
    file << "NJ=" << nj << "\n";
    file << "NK=" << nk << "\n";
    file << "x,y,z,rho,rhoU,rhoV,rhoW,rhoE\n";

    for (size_t i=0; i<ni; ++i){
        for (size_t j=0; j<nj; ++j){
            for (size_t k=0; k<nk; ++k){
                file << _mesh.getVertex(i,j,k).x() << "," << _mesh.getVertex(i,j,k).y() << "," << _mesh.getVertex(i,j,k).z() << ",";
                file << _solution._rho(i,j,k) << "," << _solution._rhoU(i,j,k) << "," << _solution._rhoV(i,j,k) << "," << _solution._rhoW(i,j,k) << "," << _solution._rhoE(i,j,k) << "\n";
            }
        }
    }

    file.close();
}