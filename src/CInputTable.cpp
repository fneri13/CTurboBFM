#include <fstream>
#include <sstream>
#include "CInputTable.hpp"
#include <cassert>



CInputTable::CInputTable(std::string filename){
    _filename = filename;
    readCSVFile();
}


void CInputTable::readCSVFile() {
    std::ifstream file(_filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << _filename << std::endl;
        return;
    }

    std::string line;
    std::vector<FieldNames> fieldOrder;

    // Read header
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, ',')) {
            try {
                fieldOrder.push_back(FieldNameMapper::fromString(field));
            } catch (const std::invalid_argument& e) {
                std::cerr << e.what() << std::endl;
                return;
            }
        }
    }

    // Read data lines
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        size_t index = 0;

        while (std::getline(ss, item, ',')) {
            if (index >= fieldOrder.size()) break;
            try {
                FloatType value = std::stod(item);
                _fieldsMap[fieldOrder[index]].push_back(value);
            } catch (...) {
                std::cerr << "Invalid numeric value: " << item << std::endl;
            }
            ++index;
        }
    }

    file.close();
}

