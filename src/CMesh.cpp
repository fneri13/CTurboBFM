#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "../include/CMesh.hpp"

CMesh::CMesh(std::string filename) {
    _filename = filename;
    readPoints();
    allocateMemory();
    printMeshInfo();
}

void CMesh::readPoints() {
    std::ifstream file(_filename);
    if (!file) {
        std::cerr << "Error opening coordinates CSV file!\n";
        std::exit(EXIT_FAILURE);  // EXIT_FAILURE = standard failure code
    }

    // Read header lines
    std::string line;
    while (std::getline(file, line)) {
        if (line.rfind("NDIMENSIONS=", 0) == 0)
            _nDimensions = std::stoi(line.substr(12));
        else if (line.rfind("NI=", 0) == 0)
            _nPointsI = std::stoi(line.substr(3));
        else if (line.rfind("NJ=", 0) == 0)
            _nPointsJ = std::stoi(line.substr(3));
        else if (line.rfind("NK=", 0) == 0)
            _nPointsK = std::stoi(line.substr(3));
        else if (line == "x,y,z") // Column header line
            break;
    }
    _nPointsTotal = _nPointsI * _nPointsJ * _nPointsK;
    _vertices.resize(_nPointsI, std::vector<std::vector<Point3D>>(_nPointsJ, std::vector<Point3D>(_nPointsK)));

    // Read coordinate data
    unsigned long int nPoint = 0;
    unsigned long int i,j,k;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        Point3D point;
        char comma;
        if (ss >> point.x >> comma >> point.y >> comma >> point.z){
            i = nPoint / (_nPointsJ * _nPointsK);
            j = (nPoint % (_nPointsJ * _nPointsK)) / _nPointsK;
            k = nPoint % _nPointsK;
            _vertices[i][j][k] = point;
            nPoint ++;
        }
           
    }

    file.close();
}

void CMesh::allocateMemory() {
    _surfacesI.resize(_nPointsI, _nPointsJ, _nPointsK, 3);
    _surfacesJ.resize(_nPointsI, _nPointsJ, _nPointsK, 3);
    _surfacesK.resize(_nPointsI, _nPointsJ, _nPointsK, 3);
    _volumes.resize(_nPointsI, _nPointsJ, _nPointsK);
    _centers.resize(_nPointsI, _nPointsJ, _nPointsK, 3);
}



void CMesh::printMeshInfo() {
    std::cout << "=========================================\n";
    std::cout << "        INFORMATION OF GRID FILE        \n";
    std::cout << "=========================================\n";
    std::cout << std::endl;

    std::cout << "Number of Dimensions    : " << _nDimensions << "\n";
    std::cout << "Points along Axis I     : " << _nPointsI << "\n";
    std::cout << "Points along Axis J     : " << _nPointsJ << "\n";
    std::cout << "Points along Axis K     : " << _nPointsK << "\n";
    std::cout << std::endl;

    std::cout << "The grid file '" << _filename << "' contains a total of " << _nPointsTotal << " points.\n";
    std::cout << "=========================================\n";
}