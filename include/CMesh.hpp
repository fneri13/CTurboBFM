#pragma once
#include <string>
#include <vector>
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


    private:
        Config _config;

        unsigned long int _nPointsI, _nPointsJ, _nPointsK, _nPointsTotal;

        unsigned long int _nDualPointsI, _nDualPointsJ, _nDualPointsK, _nDualPointsTotal;
        
        unsigned short int _nDimensions;

        Matrix3D<Vector3D> _surfacesI, _surfacesJ, _surfacesK;  // array of interfaces

        Matrix3D<FloatType> _volumes; // array of cell volumes

        std::string _filename; // name of the CSV file containing the mesg

        Matrix3D<Vector3D> _vertices; // array of vertices coordinates

        Matrix3D<Vector3D> _dualNodes; // array of dual nodes coordinates

        Matrix3D<Vector3D> _centersI, _centersJ, _centersK; // array of cell center coordinates

        void readPoints();

        void allocateMemory();

        void computeDualGrid2D();
        
        void computeDualGrid3D();

        void computeMeshInterfaces();

        void computeMeshVolumes();
        
        void computeMeshQuality();

        void computeBoundaryAreas();
        
        void printMeshInfo();

        void computeSurfaceVectorAndCG(const Vector3D &p1, const Vector3D &p2, const Vector3D &p3, const Vector3D &p4, Vector3D &normal, Vector3D &center);
};