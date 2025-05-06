#pragma once

#include <map>
#include <string>
#include "types.hpp"


/** 
*  \brief     Class handling base input capabilities.
*  \details   interface for input reading of csv files.
*  \author    Francesco Neri
*/

class CInput {
    
    public:
        CInput() {};

        CInput(std::string filename);

        Matrix3D<FloatType> getField(FieldNames fieldName) {
            return _fieldsMap[fieldName];
        }

        size_t getNumberPointsI() const {return _ni;};
        
        size_t getNumberPointsJ() const {return _nj;};
        
        size_t getNumberPointsK() const {return _nk;};

        /** \brief get the coordinates of a point (i,j,k)*/
        Vector3D getCoordinates(size_t i, size_t j, size_t k) const;

        /** \brief get the value of a field at a point (i,j,k)*/
        FloatType getField(FieldNames fieldName, size_t i, size_t j, size_t k);

    private:
        std::string _filename;
        
        std::map<FieldNames, Matrix3D<FloatType>> _fieldsMap;
        
        size_t _ni, _nj, _nk;
        

        /** \brief read a csv file
         *  \details read a csv file and store the data in a map of Matrix3D<FloatType>
        */
        void readCSVFile();

};