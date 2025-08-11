#include "COutputCSV.hpp"

void COutputCSV::writeSolution(size_t iterationCounter){
    std::string filename = getOutputFilename(iterationCounter);

    std::ofstream file(_outputDirectory + "/" + filename + ".csv");

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    std::map<std::string, Matrix3D<FloatType>> scalarFieldsMap;
    getScalarFieldsMap(scalarFieldsMap);

    size_t ni = _mesh.getNumberPointsI();
    size_t nj = _mesh.getNumberPointsJ();
    size_t nk = _mesh.getNumberPointsK();

    file << "NI=" << ni << "\n";
    file << "NJ=" << nj << "\n";
    file << "NK=" << nk << "\n";

    // write coordinates header
    file << "x,y,z";

    // write scalar fields header
    for (auto& field : scalarFieldsMap){
        file << "," << field.first;
    }
    file << "\n";

    for (size_t i=0; i<ni; ++i){
        for (size_t j=0; j<nj; ++j){
            for (size_t k=0; k<nk; ++k){
                
                //write coords
                file << _mesh.getVertex(i,j,k).x() << "," << _mesh.getVertex(i,j,k).y() << "," << _mesh.getVertex(i,j,k).z() ;
                
                // write scalar fields in the map
                for (auto& field : scalarFieldsMap){
                    file << "," << field.second(i,j,k);
                }
                file << "\n";

            }
        }
    }

    file.close();

    std::cout << std::endl;
    std::cout << "Solution written to file: " << filename << std::endl;
    std::cout << std::endl;}