#include <fstream>
#include <sstream>
#include "CInput.hpp"
#include <cassert>



CInput::CInput(std::string filename){
    _filename = filename;
    readCSVFile();
}


void CInput::readCSVFile() {
    std::ifstream file(_filename);
    if (!file) {
        std::cerr << "Error opening coordinates CSV file '" << _filename << "'!\n";
        std::exit(EXIT_FAILURE);  // EXIT_FAILURE = standard failure code
    }

    // Read header lines
    std::string line;
    std::vector<std::string> headers;
    while (std::getline(file, line)) {
        if (line.rfind("NI=", 0) == 0)
            _ni = std::stoi(line.substr(3));
        else if (line.rfind("NJ=", 0) == 0)
            _nj = std::stoi(line.substr(3));
        else if (line.rfind("NK=", 0) == 0)
            _nk = std::stoi(line.substr(3));
        else if (line.rfind("x,y,z", 0) == 0) {  
            std::istringstream headerStream(line);
            std::string column;
            while (std::getline(headerStream, column, ',')) {
                headers.push_back(column);
            }            
            break;  // Exit the loop after capturing the header
        }
    }

    // allocate the memory for the arrays
    for (auto &head: headers){
        FieldNames field = FieldNameMapper::fromString(head);
        _fieldsMap[field] = Matrix3D<FloatType>(_ni, _nj, _nk);
    }

    // Read data
    size_t nPoint = 0;
    size_t i, j, k;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<FloatType> values;

        // Split line into tokens
        while (std::getline(ss, token, ',')) {
            values.push_back(static_cast<FloatType>(std::stod(token)));
        }

        // Skip incomplete lines
        if (values.size() != headers.size()) {
            continue;
        }

        // Compute i, j, k indices
        i = nPoint / (_nj * _nk);
        j = (nPoint % (_nj * _nk)) / _nk;
        k = nPoint % _nk;
        nPoint++;

        // Assign each value to the correct field matrix
        for (size_t idx = 0; idx < headers.size(); ++idx) {
            FieldNames field = FieldNameMapper::fromString(headers[idx]);
            _fieldsMap[field](i, j, k) = values[idx];
        }
    }

    // check if the number of points read is equal to the total
    assert(nPoint == _ni * _nj * _nk && "Number of points in the input file is not equal to ni*nj*nk");

    file.close();
}


Vector3D CInput::getCoordinates(size_t i, size_t j, size_t k) const {
    return Vector3D(_fieldsMap.at(FieldNames::X_COORDS)(i, j, k),
                    _fieldsMap.at(FieldNames::Y_COORDS)(i, j, k),
                    _fieldsMap.at(FieldNames::Z_COORDS)(i, j, k));
};


FloatType CInput::getField(FieldNames fieldName, size_t i, size_t j, size_t k){
    auto it = _fieldsMap.find(fieldName);
    if (it != _fieldsMap.end()) {
        return it->second(i, j, k);
    } else {
        // Handle missing key case appropriately
        throw std::invalid_argument("FieldName not found in _fieldsMap.");
    }
};
