#pragma once

#include <map>
#include <string>
#include "types.hpp"


/** 
*  \brief     Class handling input information from a csv table of quantities
*/
class InputTable {
    
public:
    InputTable() = default;

    InputTable(std::string filename);

    std::vector<FloatType> getField(InputField fieldName) const {
        return _fieldsMap.at(fieldName);
    }

protected:    
    void readCsvFile();

private:
    std::string _filename;
    std::map<InputField, std::vector<FloatType>> _fieldsMap;        

};