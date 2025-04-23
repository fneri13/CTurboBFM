#include "../include/commonFunctions.hpp"
#include <iostream>

FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces, const Matrix2D<Vector3D> &vectors){
    auto ni = surfaces.sizeI();
    auto nj = surfaces.sizeJ();
    FloatType sum = 0.0;

    for (int i=0; i<ni; i++){
        for (int j=0; j<nj; j++){
            sum += surfaces(i,j).dot(vectors(i,j));
        }
    }
    return sum;
}

FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces){
    auto ni = surfaces.sizeI();
    auto nj = surfaces.sizeJ();
    FloatType sum = 0.0;

    for (int i=0; i<ni; i++){
        for (int j=0; j<nj; j++){
            sum += surfaces(i,j).magnitude();
        }
    }
    return sum;
}

std::array<FloatType, 5> getEulerPrimitiveFromConservative(std::array<FloatType, 5> conservative){
    std::array<FloatType, 5> primitive;
    primitive[0] = conservative[0];
    primitive[1] = conservative[1] / conservative[0];
    primitive[2] = conservative[2] / conservative[0];
    primitive[3] = conservative[3] / conservative[0];
    primitive[4] = conservative[4] / conservative[0];
    return primitive;
}

std::array<FloatType, 5> getEulerConservativeFromPrimitive(std::array<FloatType, 5> primitive){
    std::array<FloatType, 5> conservative;
    conservative[0] = primitive[0];
    conservative[1] = primitive[1] * conservative[0];
    conservative[2] = primitive[2] * conservative[0];
    conservative[3] = primitive[3] * conservative[0];
    conservative[4] = primitive[4] * conservative[0];
    return conservative;
}
