#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>
#include "CMesh.hpp"
#include "Config.hpp"
#include "commonFunctions.hpp"

CMesh::CMesh(Config& config) : _config(config) {
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
    CInput input(_config.gridFilePath());

    _nPointsI = input.getNumberPointsI();
    _nPointsJ = input.getNumberPointsJ();
    _nPointsK = input.getNumberPointsK();
    _nPointsTotal = _nPointsI * _nPointsJ * _nPointsK;
    _nDualPointsI = _nPointsI + 1;
    _nDualPointsJ = _nPointsJ + 1;
    _nDualPointsK = _nPointsK + 1;
    _nDualPointsTotal = _nDualPointsI * _nDualPointsJ * _nDualPointsK;

    _vertices.resize(_nPointsI, _nPointsJ, _nPointsK);

    for (size_t i = 0; i < _nPointsI; i++) {
        for (size_t j = 0; j < _nPointsJ; j++) {
            for (size_t k = 0; k < _nPointsK; k++) {
                _vertices(i,j,k) = input.getCoordinates(i,j,k);
            }
        }
    }
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
    std::cout << std::endl;
    std::cout << "=========================================\n";
    std::cout << "        INFORMATION OF GRID FILE        \n";
    std::cout << "=========================================\n";
    std::cout << std::endl;

    std::cout << "Number of Dimensions    : " << _nDimensions << "\n";
    std::cout << "Points along Axis I     : " << _nPointsI << "\n";
    std::cout << "Points along Axis J     : " << _nPointsJ << "\n";
    std::cout << "Points along Axis K     : " << _nPointsK << "\n";
    std::cout << std::endl;

    std::cout << "The grid file '" << _config.gridFilePath() << "' contains a total of " << _nPointsTotal << " points.\n";
    std::cout << "Aspect ratio ";
    _aspectRatioStats.printInfo();
    std::cout << "Skewness ";
    _skewnessStats.printInfo();
    std::cout << "Orthogonality ";
    _orthogonalityStats.printInfo();
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
    auto aspectRatios = computeAspectRatio();
    _aspectRatioStats = Statistics(aspectRatios);

    std::vector<FloatType> skewness, orthogonality;
    computeSkewnessAndOrthogonality(skewness, orthogonality);
    _skewnessStats = Statistics(skewness);
    _orthogonalityStats = Statistics(orthogonality);
    
}

void CMesh::computeSkewnessAndOrthogonality(std::vector<FloatType> &skewness, std::vector<FloatType> &orthogonality){
    // internal faces along i --> do it in a single function???
    for (size_t i=1; i<_nPointsI-1; i++){
        for (size_t j=0; j<_nPointsJ; j++){
            for (size_t k=0; k<_nPointsK; k++){
                Vector3D point0 = _vertices(i-1,j,k);
                Vector3D point1 = _vertices(i,j,k);
                Vector3D midPoint = (point0 + point1) / 2.0;
                Vector3D surfaceCenter = _centersI(i,j,k);
                FloatType l1 = (midPoint - surfaceCenter).magnitude();
                FloatType l2 = (point1 - point0).magnitude();
                skewness.push_back(l1 / l2);
                Vector3D surface = _surfacesI(i,j,k);
                Vector3D pointToPoint = point1 - point0;
                FloatType angle = computeAngleBetweenVectors(surface, pointToPoint);
                orthogonality.push_back(angle);
            }
        }
    }

    // internal faces along j
    for (size_t i=0; i<_nPointsI; i++){
        for (size_t j=1; j<_nPointsJ-1; j++){
            for (size_t k=0; k<_nPointsK; k++){
                Vector3D point0 = _vertices(i,j-1,k);
                Vector3D point1 = _vertices(i,j,k);
                Vector3D midPoint = (point0 + point1) / 2.0;
                Vector3D surfaceCenter = _centersJ(i,j,k);
                FloatType l1 = (midPoint - surfaceCenter).magnitude();
                FloatType l2 = (point1 - point0).magnitude();
                skewness.push_back(l1 / l2);
                Vector3D surface = _surfacesJ(i,j,k);
                Vector3D pointToPoint = point1 - point0;
                FloatType angle = computeAngleBetweenVectors(surface, pointToPoint);
                orthogonality.push_back(angle);
            }
        }
    }

    // internal faces along k
    for (size_t i=0; i<_nPointsI; i++){
        for (size_t j=0; j<_nPointsJ; j++){
            for (size_t k=1; k<_nPointsK-1; k++){
                Vector3D point0 = _vertices(i,j,k-1);
                Vector3D point1 = _vertices(i,j,k);
                Vector3D midPoint = (point0 + point1) / 2.0;
                Vector3D surfaceCenter = _centersK(i,j,k);
                FloatType l1 = (midPoint - surfaceCenter).magnitude();
                FloatType l2 = (point1 - point0).magnitude();
                skewness.push_back(l1 / l2);
                Vector3D surface = _surfacesK(i,j,k);
                Vector3D pointToPoint = point1 - point0;
                FloatType angle = computeAngleBetweenVectors(surface, pointToPoint);
                orthogonality.push_back(angle);
            }
        }
    }
}


Matrix3D<FloatType> CMesh::computeAspectRatio() {
    Matrix3D<FloatType> aspectRatio(_nPointsI, _nPointsJ, _nPointsK);
    for (size_t i=0; i<_nPointsI; i++){
        for (size_t j=0; j<_nPointsJ; j++){
            for (size_t k=0; k<_nPointsK; k++){
                Vector3D point0 = _dualNodes(i,j,k);
                Vector3D pointI = _dualNodes(i+1,j,k);
                Vector3D pointJ = _dualNodes(i,j+1,k);
                Vector3D pointK = _dualNodes(i,j,k+1);
                FloatType iEdge = (pointI - point0).magnitude();
                FloatType jEdge = (pointJ - point0).magnitude();
                FloatType kEdge = (pointK - point0).magnitude();
                if (_nDimensions == 2) {
                    FloatType maxEdge = std::max(iEdge, jEdge);
                    FloatType minEdge = std::min(iEdge, jEdge);
                    aspectRatio(i, j, k) = maxEdge / minEdge;
                } else {
                    FloatType maxEdge = std::max({iEdge, jEdge, kEdge});
                    FloatType minEdge = std::min({iEdge, jEdge, kEdge});
                    aspectRatio(i, j, k) = maxEdge / minEdge;
                }
            }
        }
    }
    return aspectRatio;
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

    _boundarySurfaces[BoundaryIndices::I_START] = _surfacesI.getBoundarySlice(BoundaryIndices::I_START);
    _boundarySurfaces[BoundaryIndices::I_END] = _surfacesI.getBoundarySlice(BoundaryIndices::I_END);
    _boundarySurfaces[BoundaryIndices::J_START] = _surfacesJ.getBoundarySlice(BoundaryIndices::J_START);
    _boundarySurfaces[BoundaryIndices::J_END] = _surfacesJ.getBoundarySlice(BoundaryIndices::J_END);
    _boundarySurfaces[BoundaryIndices::K_START] = _surfacesK.getBoundarySlice(BoundaryIndices::K_START);
    _boundarySurfaces[BoundaryIndices::K_END] = _surfacesK.getBoundarySlice(BoundaryIndices::K_END);

    // compute the areas for each boundary
    for (const auto& index : BoundaryIndicesArray) {
        _boundaryAreas[index] = computeSurfaceIntegral(_boundarySurfaces[index]);
    }

}

const Matrix3D<Vector3D>& CMesh::getSurfaces(FluxDirection direction) const {
    switch (direction) {
        case FluxDirection::I: return _surfacesI;
        case FluxDirection::J: return _surfacesJ;
        case FluxDirection::K: return _surfacesK;
        default:
            throw std::runtime_error("Invalid flux direction.");
    }
}

const Matrix3D<Vector3D>& CMesh::getMidPoints(FluxDirection direction) const {
    switch (direction) {
        case FluxDirection::I: return _centersI;
        case FluxDirection::J: return _centersJ;
        case FluxDirection::K: return _centersK;
        default:
            throw std::runtime_error("Invalid flux direction.");
    }
}


void CMesh::computeDualGrid2D() {
    assert (_nDimensions == 2 && "Can only compute a 2D dual grid for a 2D mesh.");

    std::cout << "Computing dual grid coordinates for 2D mesh\n";

    Matrix2D<Vector3D> nodes(_nDualPointsI, _nDualPointsJ);

    // find internal dual nodes
    for (int i = 1; i < _nDualPointsI-1; i++) {
        for (int j = 1; j < _nDualPointsJ-1; j++) {
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


void CMesh::computeDualGrid3D() {
    assert (_nDimensions == 3 && "Can only compute a 3D dual grid for a 3D mesh.");
    std::cout << "Computing dual grid coordinates for 3D mesh\n";

    // internal dual nodes
    for (size_t i = 1; i < _nDualPointsI-1; i++) {
        for (size_t j = 1; j < _nDualPointsJ-1; j++) {
            for (size_t k = 1; k < _nDualPointsK-1; k++){
                _dualNodes(i,j,k) = (_vertices(i-1,j-1,k-1) + _vertices(i,j-1,k-1) + _vertices(i-1,j,k-1) + _vertices(i-1,j-1,k) + \
                                     _vertices(i  ,j,  k-1) + _vertices(i,j-1,k)   + _vertices(i-1,j,k)   + _vertices(i  ,j  ,k) ) / 8.0;
            }
        }
    }
    
    // corner dual nodes
    _dualNodes(0,0,0) = _vertices(0,0,0);
    _dualNodes(_nDualPointsI-1, 0, 0) = _vertices(_nPointsI-1, 0, 0);
    _dualNodes(0, _nDualPointsJ-1, 0) = _vertices(0, _nPointsJ-1, 0);
    _dualNodes(0, 0, _nDualPointsK-1) = _vertices(0, 0, _nPointsK-1);
    _dualNodes(_nDualPointsI-1, _nDualPointsJ-1, 0) = _vertices(_nPointsI-1, _nPointsJ-1, 0);
    _dualNodes(_nDualPointsI-1, 0, _nDualPointsK-1) = _vertices(_nPointsI-1, 0, _nPointsK-1);
    _dualNodes(0, _nDualPointsJ-1, _nDualPointsK-1) = _vertices(0, _nPointsJ-1, _nPointsK-1);
    _dualNodes(_nDualPointsI-1, _nDualPointsJ-1, _nDualPointsK-1) = _vertices(_nPointsI-1, _nPointsJ-1, _nPointsK-1);


    // edges
    // i-oriented
    for (size_t i=1; i<_nDualPointsI-1; i++) {
        _dualNodes(i,0,0) = (_vertices(i-1,0,0) + _vertices(i,0,0)) / 2.0;
        _dualNodes(i,_nDualPointsJ-1,0) = (_vertices(i-1,_nPointsJ-1,0) + _vertices(i,_nPointsJ-1,0)) / 2.0;
        _dualNodes(i,0,_nDualPointsK-1) = (_vertices(i-1,0,_nPointsK-1) + _vertices(i,0,_nPointsK-1)) / 2.0;
        _dualNodes(i,_nDualPointsJ-1,_nDualPointsK-1) = (_vertices(i-1,_nPointsJ-1,_nPointsK-1) + _vertices(i,_nPointsJ-1,_nPointsK-1)) / 2.0;
    }
    // j-oriented
    for (size_t j=1; j<_nDualPointsJ-1; j++) {
        _dualNodes(0,j,0) = (_vertices(0,j-1,0) + _vertices(0,j,0)) / 2.0;
        _dualNodes(_nDualPointsI-1,j,0) = (_vertices(_nPointsI-1,j-1,0) + _vertices(_nPointsI-1,j,0)) / 2.0;
        _dualNodes(0,j,_nDualPointsK-1) = (_vertices(0,j-1,_nPointsK-1) + _vertices(0,j,_nPointsK-1)) / 2.0;
        _dualNodes(_nDualPointsI-1,j,_nDualPointsK-1) = (_vertices(_nPointsI-1,j-1,_nPointsK-1) + _vertices(_nPointsI-1,j,_nPointsK-1)) / 2.0;
    }
    // k-oriented
    for (size_t k=1; k<_nDualPointsK-1; k++) {
        _dualNodes(0,0,k) = (_vertices(0,0,k-1) + _vertices(0,0,k)) / 2.0;
        _dualNodes(_nDualPointsI-1,0,k) = (_vertices(_nPointsI-1,0,k-1) + _vertices(_nPointsI-1,0,k)) / 2.0;
        _dualNodes(0,_nDualPointsJ-1,k) = (_vertices(0,_nPointsJ-1,k-1) + _vertices(0,_nPointsJ-1,k)) / 2.0;
        _dualNodes(_nDualPointsI-1,_nDualPointsJ-1,k) = (_vertices(_nPointsI-1,_nPointsJ-1,k-1) + _vertices(_nPointsI-1,_nPointsJ-1,k)) / 2.0;
    }
    

    //boundaries
    // i-faces
    for (size_t j=1; j<_nDualPointsJ-1; j++) {
        for (size_t k=1; k<_nDualPointsK-1; k++) {
            _dualNodes(0,j,k) = (_vertices(0,j-1,k-1) + _vertices(0,j,k-1) + _vertices(0,j-1,k) + _vertices(0,j,k)) / 4.0;
            _dualNodes(_nDualPointsI-1,j,k) = (_vertices(_nPointsI-1,j-1,k-1) + _vertices(_nPointsI-1,j,k-1) + _vertices(_nPointsI-1,j-1,k) + _vertices(_nPointsI-1,j,k)) / 4.0;
        }
    }
    // j-faces
    for (size_t i=1; i<_nDualPointsI-1; i++) {
        for (size_t k=1; k<_nDualPointsK-1; k++) {
            _dualNodes(i,0,k) = (_vertices(i-1,0,k-1) + _vertices(i,0,k-1) + _vertices(i-1,0,k) + _vertices(i,0,k)) / 4.0;
            _dualNodes(i,_nDualPointsJ-1,k) = (_vertices(i-1,_nPointsJ-1,k-1) + _vertices(i,_nPointsJ-1,k-1) + _vertices(i-1,_nPointsJ-1,k) + _vertices(i,_nPointsJ-1,k)) / 4.0;
        }
    }
    // k-faces
    for (size_t i=1; i<_nDualPointsI-1; i++) {
        for (size_t j=1; j<_nDualPointsJ-1; j++) {
            _dualNodes(i,j,0) = (_vertices(i-1,j-1,0) + _vertices(i,j-1,0) + _vertices(i-1,j,0) + _vertices(i,j,0)) / 4.0;
            _dualNodes(i,j,_nDualPointsK-1) = (_vertices(i-1,j-1,_nPointsK-1) + _vertices(i,j-1,_nPointsK-1) + _vertices(i-1,j,_nPointsK-1) + _vertices(i,j,_nPointsK-1)) / 4.0;
        }
    }

    
}

void CMesh::getElementEdges(size_t i, size_t j, size_t k, Vector3D &iEdge, Vector3D &jEdge, Vector3D &kEdge) const {
    Vector3D pt0, ptI, ptJ, ptK;
    pt0 = _dualNodes(i, j, k);
    ptI = _dualNodes(i+1, j, k);
    ptJ = _dualNodes(i, j+1, k);
    ptK = _dualNodes(i, j, k+1);
    iEdge = ptI - pt0;
    jEdge = ptJ - pt0;
    kEdge = ptK - pt0;
}
