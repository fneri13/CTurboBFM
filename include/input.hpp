#pragma once

#include <map>
#include <string>
#include "types.hpp"


/** 
*  \brief Class handling input csv files
*/

class Input {
    
public:
    Input() = default;

    Input(std::string filename);

    Matrix3D<FloatType> getField(FieldNames fieldName) const {
        return _fieldsMap.at(fieldName);
    }

    size_t getNumberPointsI() const {return _ni;};
    
    size_t getNumberPointsJ() const {return _nj;};
    
    size_t getNumberPointsK() const {return _nk;};

    Vector3D getCoordinates(size_t i, size_t j, size_t k) const;

    FloatType getField(FieldNames fieldName, size_t i, size_t j, size_t k) const ;

protected:
    void readCsvFile();

private:
    std::string _filename;
    std::map<FieldNames, Matrix3D<FloatType>> _fieldsMap;
    size_t _ni, _nj, _nk;

};