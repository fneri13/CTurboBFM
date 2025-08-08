#pragma once

#include <map>
#include <string>
#include "types.hpp"


/** 
*  \brief     Class handling base input table capabilities.
*  \details   interface for input reading of csv table input files.
*  \author    Francesco Neri
*/

class CInputTable {
    
    public:
        CInputTable() {};

        CInputTable(std::string filename);

        std::vector<FloatType> getField(FieldNames fieldName) const {
            return _fieldsMap.at(fieldName);
        }

    private:
        std::string _filename;
        
        std::map<FieldNames, std::vector<FloatType>> _fieldsMap;        

        void readCSVFile();

};