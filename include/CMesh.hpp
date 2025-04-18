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
        std::string _filename;
        ArrayPoints _points;

        /// Read the points from the CSV file
        void readPoints();

        /// Compute the mesh interfaces based on vertex centered approach
        void computeMeshInterfaces();

        /// Compute the volumes of the mesh cells
        void computeMeshVolumes();
        
        void computeMeshQuality();
        void computeBoundaryAreas();
        void printMeshInfo();
};