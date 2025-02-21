#pragma once
#include <string>
#include <vector>

struct Point3D {
    double x, y, z;
};

class CMesh {
    using arrayPoints = std::vector<std::vector<std::vector<Point3D>>>;
    public:
        CMesh(std::string filename);
        ~CMesh() {};
    private:
        unsigned long int _nPointsI, _nPointsJ, _nPointsK, _nPointsTotal;
        unsigned short int _nDimensions;
        std::string _filename;
        arrayPoints _points;
        void readPoints();
};