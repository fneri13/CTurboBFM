#include <iostream>
#include "../../include/types.hpp"

int main() {
    Matrix3D<FloatType> mat(4, 3, 5);

    std::cout << "mat(2,1,0) = " << mat(2, 1, 0) << std::endl;
    mat(2, 1, 0) = 42.5;
    std::cout << "mat(2,1,0) = " << mat(2, 1, 0) << std::endl;

    size_t ni = mat.sizeI();
    size_t nj = mat.sizeJ();
    size_t nk = mat.sizeK();
    std::cout << "Matrix dimensions: " << ni << " x " << nj << " x " << nk << std::endl;
    std::cout << "Total size: " << ni * nj * nk << std::endl;

    for (size_t i = 0; i < ni; i++) {
        for (size_t j = 0; j < nj; j++) {
            for (size_t k = 0; k < nk; k++) {
                std::cout << mat(i, j, k) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    try {
        mat(4, 0, 5) = 1.0;  // will throw
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    return 0;
}
