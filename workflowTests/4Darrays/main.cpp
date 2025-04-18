#include <iostream>
#include "../../include/types.hpp"

int main() {
    Matrix4D<FloatType> mat(4, 3, 5, 1);

    std::cout << "mat(2,1,0,0) = " << mat(2, 1, 0, 0) << std::endl;
    mat(2, 1, 0, 0) = 42.5;
    std::cout << "mat(2,1,0,0) = " << mat(2, 1, 0, 0) << std::endl;

    try {
        mat(4, 0, 5, 1) = 1.0;  // will throw
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    return 0;
}
