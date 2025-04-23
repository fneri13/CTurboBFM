#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>
#include "../include/CMesh.hpp"
#include "../include/Config.hpp"
#include "../include/commonFunctions.hpp"

CMesh::CMesh(Config& config) : _config(config) {
    _filename = config.gridFilePath();
    readPoints();
    allocateMemory();
    auto topology = config.getTopology();
    if (topology == Topology::TWO_DIMENSIONAL || topology == Topology::AXISYMMETRIC) {
        computeDualGrid2D();
    } 
    else {
        // computeDualGrid3D();
    }
    computeMeshInterfaces();
    computeMeshVolumes();
    // computeMeshQuality();
    computeBoundaryAreas();
    // printMeshInfo();
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
    _dualNodes.resize(_nDualPointsI, _nDualPointsJ, _nDualPointsK);

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
                computeSurfaceVectorAndCG(_dualNodes(i,j,k), _dualNodes(i,j+1,k), _dualNodes(i,j+1,k+1), _dualNodes(i,j,k+1), _surfacesI(i,j,k), _centersI(i,j,k));
            }
        }
    }

    // Computing the j-interfaces
    for (int i = 0; i < _nDualPointsI-1; i++) {
        for (int j = 0; j < _nDualPointsJ; j++) {
            for (int k = 0; k < _nDualPointsK-1; k++) {
                computeSurfaceVectorAndCG(_dualNodes(i,j,k), _dualNodes(i,j,k+1), _dualNodes(i+1,j,k+1), _dualNodes(i+1,j,k), _surfacesJ(i,j,k), _centersJ(i,j,k));
            }
        }
    }

    // Computing the k-interfaces
    for (int i = 0; i < _nDualPointsI-1; i++) {
        for (int j = 0; j < _nDualPointsJ-1; j++) {
            for (int k = 0; k < _nDualPointsK; k++) {
                computeSurfaceVectorAndCG(_dualNodes(i,j,k), _dualNodes(i+1,j,k), _dualNodes(i+1,j+1,k), _dualNodes(i,j+1,k), _surfacesK(i,j,k), _centersK(i,j,k));
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
    surf1 = a1.cross(b1) / 2.0;
    surf2 = a2.cross(b2) / 2.0;
    normal = surf1 + surf2;

    // compute centers
    Vector3D center1, center2;
    center1 = (p1 + p2 + p3) / 3.0;
    center2 = (p1 + p3 + p4) / 3.0;
    center = (center1 * surf1.magnitude() + center2 * surf2.magnitude()) / (surf1.magnitude() + surf2.magnitude());
}

FloatType CMesh::computeElementVolume(const std::vector<Vector3D> &boundSurfaces, const std::vector<Vector3D> &boundCenters) {
    // Compute the volume using the Green-Gauss theorem
    FloatType volume {0.0};
    for (size_t i = 0; i < boundSurfaces.size(); ++i) {
        volume += boundSurfaces[i].x() * boundCenters[i].x();
    }
    return volume;
}

void CMesh::computeMeshVolumes() {
    std::cout << "Compute element volumes\n";

    for (int i = 0; i< _nPointsI; i++) {
        for (int j = 0; j<_nPointsJ; j++) {
            for (int k = 0; k<_nPointsK; k++) {
                std::vector<Vector3D> boundSurfaces (6);
                std::vector<Vector3D> boundCenters (6);
        
                boundSurfaces[0] = - _surfacesI(i,j,k);
                boundSurfaces[1] = - _surfacesJ(i,j,k);
                boundSurfaces[2] = - _surfacesK(i,j,k);
                boundSurfaces[3] = _surfacesI(i+1,j,k);
                boundSurfaces[4] = _surfacesJ(i,j+1,k);
                boundSurfaces[5] = _surfacesK(i,j,k+1);
                boundCenters[0] = _centersI(i,j,k);
                boundCenters[1] = _centersJ(i,j,k);
                boundCenters[2] = _centersK(i,j,k); 
                boundCenters[3] = _centersI(i+1,j,k);
                boundCenters[4] = _centersJ(i,j+1,k);
                boundCenters[5] = _centersK(i,j,k+1);

                // Compute the volume using the Green-Gauss theorem
                _volumes(i,j,k) = computeElementVolume(boundSurfaces, boundCenters);
            }
        }
    }
}

void CMesh::computeMeshQuality() {
    std::cout << "Copy from TurboBFM\n";
}


void CMesh::computeBoundaryAreas() {
    std::array<BoundaryIndices, 6> BoundaryIndicesArray = {
        BoundaryIndices::I_START,
        BoundaryIndices::I_END,
        BoundaryIndices::J_START,
        BoundaryIndices::J_END,
        BoundaryIndices::K_START,
        BoundaryIndices::K_END
    };

    // get the slices of boundary areas
    for (const auto& index : BoundaryIndicesArray) {
        _boundarySurfaces[index] = getBoundarySurface(index);
    }

    // compute the areas for each boundary
    for (const auto& index : BoundaryIndicesArray) {
        _boundaryAreas[index] = computeSurfaceIntegral(_boundarySurfaces[index]);
    }

}

const Matrix2D<Vector3D> CMesh::getBoundarySurface(BoundaryIndices index) const {
    Matrix2D<Vector3D> boundary;

    // i slices
    if (index==BoundaryIndices::I_START || index==BoundaryIndices::I_END){
        auto boundaryNormals = getSurfacesI();
        int ni = boundaryNormals.sizeI();
        int nj = boundaryNormals.sizeJ();
        int nk = boundaryNormals.sizeK();
        boundary.resize(nj, nk);

        if (index==BoundaryIndices::I_START){
            for (int j=0; j<nj; j++){
                for (int k=0; k<nk; k++){
                    boundary(j,k) = boundaryNormals(0, j, k);
                }
            }
        } else {
            for (int j=0; j<nj; j++){
                for (int k=0; k<nk; k++){
                    boundary(j,k) = boundaryNormals(ni-1, j, k);
                }
            }
        }
    }

    // j slices
    else if (index==BoundaryIndices::J_START || index==BoundaryIndices::J_END){
        auto boundaryNormals = getSurfacesJ();
        int ni = boundaryNormals.sizeI();
        int nj = boundaryNormals.sizeJ();
        int nk = boundaryNormals.sizeK();
        boundary.resize(ni, nk);

        if (index==BoundaryIndices::J_START){
            for (int i=0; i<ni; i++){
                for (int k=0; k<nk; k++){
                    boundary(i,k) = boundaryNormals(i, 0, k);
                }
            }
        } else {
            for (int i=0; i<ni; i++){
                for (int k=0; k<nk; k++){
                    boundary(i,k) = boundaryNormals(i, nj-1, k);
                }
            }
        }
    }

    // k slices
    else{
        auto boundaryNormals = getSurfacesK();
        int ni = boundaryNormals.sizeI();
        int nj = boundaryNormals.sizeJ();
        int nk = boundaryNormals.sizeK();
        boundary.resize(ni, nj);

        if (index==BoundaryIndices::K_START){
            for (int i=0; i<ni; i++){
                for (int j=0; j<nj; j++){
                    boundary(i,j) = boundaryNormals(i, j, 0);
                }
            }
        } else {
            for (int i=0; i<ni; i++){
                for (int j=0; j<nj; j++){
                    boundary(i,j) = boundaryNormals(i, j, nk-1);
                }
            }
        }
    }
    return boundary;
}

void CMesh::computeDualGrid2D() {
    std::cout << "Computing dual grid coordinates for 2D mesh\n";

    Matrix2D<Vector3D> nodes(_nDualPointsI, _nDualPointsJ);

    // find internal dual nodes
    for (int i = 1; i < _nDualPointsI-1; i++) {
        for (int j = 1; j < _nDualPointsJ-1; j++) {
            if (i==1 && j==1){
                std::cout << "this should be the problematic one\n";
            }
            nodes(i,j) = (_vertices(i-1,j-1,0) + _vertices(i,j-1,0) + _vertices(i,j,0) + _vertices(i-1,j,0)) / 4.0;
        }
    }
    
    // find the corner dual nodes
    nodes(0,0) = _vertices(0,0,0);
    nodes(0, _nDualPointsJ-1) = _vertices(0, _nPointsJ-1, 0);
    nodes(_nDualPointsI-1, 0) = _vertices(_nPointsI-1, 0, 0);
    nodes(_nDualPointsI-1, _nDualPointsJ-1) = _vertices(_nPointsI-1, _nPointsJ-1, 0);

    // find dual nodes on edges
    for (int i = 1; i < _nDualPointsI-1; i++) {
        nodes(i, 0) = (_vertices(i-1,0,0) + _vertices(i,0,0)) / 2.0;
        nodes(i, _nDualPointsJ-1) = (_vertices(i-1, _nPointsJ-1, 0) + _vertices(i, _nPointsJ-1, 0)) / 2.0;
    }
    for (int j = 1; j < _nDualPointsJ-1; j++) {
        nodes(0, j) = (_vertices(0, j-1, 0) + _vertices(0, j, 0)) / 2.0;    
        nodes(_nDualPointsI-1, j) = (_vertices(_nPointsI-1, j-1, 0) + _vertices(_nPointsI-1, j, 0)) / 2.0;    
    }

    if (_config.getTopology() == Topology::AXISYMMETRIC) {
        _wedgeAngle = 1.0 * M_PI / 180.0;
        for (int i = 0; i < _nDualPointsI; i++) {
            for (int j = 0; j < _nDualPointsJ; j++) {
                _dualNodes(i,j,0).x() = nodes(i,j).x();
                _dualNodes(i,j,1).x() = nodes(i,j).x();

                _dualNodes(i,j,0).y() = nodes(i,j).y() * cos(-_wedgeAngle/2.0);
                _dualNodes(i,j,1).y() = nodes(i,j).y() * cos(+_wedgeAngle/2.0);

                _dualNodes(i,j,0).z() = nodes(i,j).y() * sin(-_wedgeAngle/2.0);
                _dualNodes(i,j,1).z() = nodes(i,j).y() * sin(+_wedgeAngle/2.0);
            }
        }
    }
    else {
        _cellThickness = 1.0;
        for (int i = 0; i < _nDualPointsI; i++) {
            for (int j = 0; j < _nDualPointsJ; j++) {
                _dualNodes(i,j,0).x() = nodes(i,j).x();
                _dualNodes(i,j,1).x() = nodes(i,j).x();

                _dualNodes(i,j,0).y() = nodes(i,j).y();
                _dualNodes(i,j,1).y() = nodes(i,j).y();

                _dualNodes(i,j,0).z() = -_cellThickness / 2.0;
                _dualNodes(i,j,1).z() = +_cellThickness / 2.0;
            }
        }
            
    }
}

void CMesh::getElementEdges(int i, int j, int k, Vector3D &iEdge, Vector3D &jEdge, Vector3D &kEdge) const {
    Vector3D pt0, ptI, ptJ, ptK;
    pt0 = _dualNodes(i, j, k);
    ptI = _dualNodes(i+1, j, k);
    ptJ = _dualNodes(i, j+1, k);
    ptK = _dualNodes(i, j, k+1);
    iEdge = ptI - pt0;
    jEdge = ptJ - pt0;
    kEdge = ptK - pt0;
}
