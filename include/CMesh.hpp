#pragma once
#include <string>
#include <vector>
#include "types.hpp"

/** 
  *  \brief     Class handling mesh capabilities.
  *  \details   It computes the finite volume mesh.
  *  \author    Francesco Neri
  */
class CMesh {
    using ArrayPoints = std::vector<std::vector<std::vector<Point3D>>>;
    
    public:
        CMesh(std::string filename);
        ~CMesh() {};


    private:
        unsigned long int _nPointsI, _nPointsJ, _nPointsK, _nPointsTotal;
        
        unsigned short int _nDimensions;

        Matrix4D<FloatType> _surfacesI, _surfacesJ, _surfacesK;  // array of interfaces

        Matrix3D<FloatType> _volumes; // array of cell volumes

        std::string _filename; // name of the CSV file containing the mesg

        ArrayPoints _vertices; // array of vertices coordinates

        Matrix4D<FloatType> _centers; // array of cell center coordinates

        void readPoints();

        void allocateMemory();

        void computeMeshInterfaces();

        void computeMeshVolumes();
        
        void computeMeshQuality();
        void computeBoundaryAreas();
        void printMeshInfo();
};