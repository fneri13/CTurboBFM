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

    private:
        std::string _filename;
        std::map<FieldNames, Matrix3D<FloatType>> _fieldsMap;

        void readCSVFile();

};