#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "../include/CMesh.hpp"
#include "../include/Config.hpp"

CMesh::CMesh(Config &config) {
    _config = config;
    _filename = config.gridFilePath();
    readPoints();
    allocateMemory();

    auto topology = config.getTopology();
    if (topology == Topology::TWO_DIMENSIONAL || topology == Topology::AXISYMMETRIC) {
        computeDualGrid2D();
    } 
    else {
        computeDualGrid3D();
    }

    computeMeshInterfaces();
    computeMeshVolumes();
    computeMeshQuality();
    computeBoundaryAreas();
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
    _nDualPointsI = _nPointsI + 1;
    _nDualPointsJ = _nPointsJ + 1;
    _nDualPointsK = _nPointsK + 1;
    _nDualPointsTotal = _nDualPointsI * _nDualPointsJ * _nDualPointsK;

    _vertices.resize(_nPointsI, _nPointsJ, _nPointsK);

    // Read coordinate data
    unsigned long int nPoint = 0;
    unsigned long int i,j,k;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        Vector3D point;
        char comma;
        if (ss >> point.x() >> comma >> point.y() >> comma >> point.z()){
            i = nPoint / (_nPointsJ * _nPointsK);
            j = (nPoint % (_nPointsJ * _nPointsK)) / _nPointsK;
            k = nPoint % _nPointsK;
            _vertices(i,j,k) = point;
            nPoint ++;
        }
    }
    file.close();
}

void CMesh::allocateMemory() {
    _surfacesI.resize(_nDualPointsI, _nDualPointsJ-1, _nDualPointsK-1);
    _surfacesJ.resize(_nDualPointsI-1, _nDualPointsJ, _nDualPointsK-1);
    _surfacesK.resize(_nDualPointsI-1, _nDualPointsJ-1, _nDualPointsK);
    
    _centersI.resize(_nDualPointsI, _nDualPointsJ-1, _nDualPointsK-1);
    _centersJ.resize(_nDualPointsI-1, _nDualPointsJ, _nDualPointsK-1);
    _centersK.resize(_nDualPointsI-1, _nDualPointsJ-1, _nDualPointsK);

    _volumes.resize(_nPointsI, _nPointsJ, _nPointsK);
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


void CMesh::computeMeshInterfaces() {
    std::cout << "Computing cell elements interfaces\n";

    // Computing the i-interfaces
    for (int i = 0; i < _nDualPointsI; i++) {
        for (int j = 0; j < _nDualPointsJ-1; j++) {
            for (int k = 0; k < _nDualPointsK-1; k++) {
                computeSurfaceVectorAndCG(_vertices(i,j,k), _vertices(i,j+1,k), _vertices(i,j+1,k+1), _vertices(i,j,k+1), _surfacesI(i,j,k), _centersI(i,j,k));
            }
        }
    }


}

void CMesh::computeSurfaceVectorAndCG(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3, const Vector3D &p4, Vector3D &normal, Vector3D &center){
    // sides used to compute areas
    Vector3D a1, b1, a2, b2;
    a1 = p2 - p1;
    b1 = p3 - p2;
    a2 = p3 - p1;
    b2 = p4 - p3;

    // compute areas
    Vector3D surf1, surf2;
    surf1 = a1.cross(b1);
    surf2 = a2.cross(b2);
    normal = surf1 + surf2;
    normal *= 0.5;

    // compute centers
    Vector3D center1, center2;
    center1 = (p1 + p2 + p3) / 3.0;
    center2 = (p1 + p3 + p4) / 3.0;
    center = (center1 * surf1.magnitude() + center2 * surf2.magnitude()) / (surf1.magnitude() + surf2.magnitude());
}

void CMesh::computeMeshVolumes() {
    std::cout << "Copy from TurboBFM\n";
}

void CMesh::computeMeshQuality() {
    std::cout << "Copy from TurboBFM\n";
}

void CMesh::computeBoundaryAreas() {
    std::cout << "Copy from TurboBFM\n";
}