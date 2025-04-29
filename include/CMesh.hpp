#pragma once

#include "types.hpp"
#include "Config.hpp"

/** 
  *  \brief     Class handling mesh capabilities.
  *  \details   It computes the finite volume mesh.
  *  \author    Francesco Neri
  */
class CMesh {
    
    public:
        CMesh(Config &config);
        ~CMesh() {};

        FloatType computeElementVolume(const std::vector<Vector3D> &boundSurfaces, const std::vector<Vector3D> &boundCenters);

        const Matrix3D<Vector3D>& getSurfaces(FluxDirection direction) const;

        const Matrix3D<Vector3D>& getMidPoints(FluxDirection direction) const;

        const Matrix3D<Vector3D>& getSurfacesI() const {return _surfacesI;}
        
        const Matrix3D<Vector3D>& getMidPointsI() const {return _centersI;}

        const Matrix3D<Vector3D>& getSurfacesJ() const {return _surfacesJ;}
        
        const Matrix3D<Vector3D>& getMidPointsJ() const {return _centersJ;}

        const Matrix3D<Vector3D>& getSurfacesK() const {return _surfacesK;}

        const Matrix3D<Vector3D>& getMidPointsK() const {return _centersK;}
        
        const Matrix3D<Vector3D> getVertices() const {return _vertices;}
        
        const Matrix3D<Vector3D> getDualNodes() const {return _dualNodes;}
        
        const Matrix3D<FloatType> getVolumes() const {return _volumes;}

        const FloatType getVolume(size_t i, size_t j, size_t k) const {return _volumes(i,j,k);}

        const size_t getNumberDimensions() const {return _nDimensions;}

        const size_t getNumberPointsI() const {return _nPointsI;}
        
        const size_t getNumberPointsJ() const {return _nPointsJ;}
        
        const size_t getNumberPointsK() const {return _nPointsK;}

        void getElementEdges(size_t i, size_t j, size_t k, Vector3D &iEdge, Vector3D &jEdge, Vector3D &kEdge) const;

        FloatType getBoundaryTotalArea(BoundaryIndices boundIndex) const {return _boundaryAreas.at(boundIndex);}

        const Matrix2D<Vector3D> getMeshBoundary(BoundaryIndices boundIndex) const {return _boundarySurfaces.at(boundIndex); }

        // compute the surface vector and center of surface given 4 points
        void computeSurfaceVectorAndCG(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3, const Vector3D &p4, Vector3D &normal, Vector3D &center);

    private:
        const Config& _config;

        unsigned long int _nPointsI, _nPointsJ, _nPointsK, _nPointsTotal;

        unsigned long int _nDualPointsI, _nDualPointsJ, _nDualPointsK, _nDualPointsTotal;
        
        unsigned short int _nDimensions;

        Matrix3D<Vector3D> _surfacesI, _surfacesJ, _surfacesK;  // array of interfaces

        Matrix3D<FloatType> _volumes; // array of cell volumes

        Matrix3D<Vector3D> _vertices; // array of vertices coordinates

        Matrix3D<Vector3D> _dualNodes; // array of dual nodes coordinates

        Matrix3D<Vector3D> _centersI, _centersJ, _centersK; // array of cell center coordinates

        FloatType _wedgeAngle {0.0}, _cellThickness {0.0}; // angle of the wedge used to compute the cell volume, or the thickness for 2d

        std::map<BoundaryIndices, FloatType> _boundaryAreas;

        std::map<BoundaryIndices, Matrix2D<Vector3D>> _boundarySurfaces;

        Statistics _aspectRatioStats, _skewnessStats, _orthogonalityStats;

        // read the coordinates from the mesh CSV file
        void readPoints();

        // allocate memory for nodes, surfaces, volumes and surface centers
        void allocateMemory();

        // compute the dual grid coordinates for a 2d mesh
        void computeDualGrid2D();
        
        // compute the dual grid coordinates for a 3d mesh
        void computeDualGrid3D();

        // compute the surface normals and centers 
        void computeMeshInterfaces();

        // compute the cell volumes
        void computeMeshVolumes();
        
        // compute the mesh quality
        void computeMeshQuality();

        // compute the areas of the boundaries
        void computeBoundaryAreas();
        
        // print some info on screen
        void printMeshInfo();

        // aspect ratio
        Matrix3D<FloatType> computeAspectRatio();

        void computeSkewnessAndOrthogonality(std::vector<FloatType> &skewness, std::vector<FloatType> &orthogonality);

};