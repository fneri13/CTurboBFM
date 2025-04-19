#include <iostream>
#include "../../include/types.hpp"

int main() {
    Vector3D point1(1.0, 2.0, 3.0);
    Vector3D point2(4.0, 5.0, 6.0);

    Vector3D result = point1.cross(point2);

    std::cout << "Cross product: " << result << std::endl;
    std::cout << "Magnitude of point1: " << point1.magnitude() << std::endl;
    std::cout << "Normalized point1: " << point1.normalized() << std::endl;
    std::cout << "point1 + point2: " << (point1 + point2) << std::endl;
    std::cout << "point1 - point2: " << (point1 - point2) << std::endl;
    std::cout << "point1 * 2: " << (point1 * 2) << std::endl;
    std::cout << "point1 / 2: " << (point1 / 2) << std::endl;
    return 0;
}
