#include <iostream>
#include "../../include/types.hpp"

int main() {
    Matrix2D<FloatType> mat(4, 3);

    std::cout << "mat(2,1) = " << mat(2, 1) << std::endl;
    mat(2, 1) = 42.5;
    std::cout << "mat(2,1) = " << mat(2, 1) << std::endl;

    try {
        mat(4, 0) = 1.0;  // will throw
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    return 0;
}
