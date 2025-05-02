#include "commonFunctions.hpp"
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

FloatType computeSurfaceIntegral(const Matrix2D<Vector3D> &surfaces, const Matrix2D<FloatType> &vecX, const Matrix2D<FloatType> &vecY, const Matrix2D<FloatType> &vecZ){
    auto ni = surfaces.sizeI();
    auto nj = surfaces.sizeJ();
    FloatType sum = 0.0;

    for (int i=0; i<ni; i++){
        for (int j=0; j<nj; j++){
            sum += surfaces(i,j).x() * vecX(i,j) + surfaces(i,j).y() * vecY(i,j) + surfaces(i,j).z() * vecZ(i,j);
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


StateVector getEulerPrimitiveFromConservative(StateVector conservative){
    StateVector primitive;
    primitive[0] = conservative[0];
    primitive[1] = conservative[1] / conservative[0];
    primitive[2] = conservative[2] / conservative[0];
    primitive[3] = conservative[3] / conservative[0];
    primitive[4] = conservative[4] / conservative[0];
    return primitive;
}


StateVector getEulerConservativeFromPrimitive(StateVector primitive){
    StateVector conservative;
    conservative[0] = primitive[0];
    conservative[1] = primitive[1] * primitive[0];
    conservative[2] = primitive[2] * primitive[0];
    conservative[3] = primitive[3] * primitive[0];
    conservative[4] = primitive[4] * primitive[0];
    return conservative;
}


StateVector computeEulerFluxFromPrimitive(StateVector primitive, Vector3D surface, CFluid fluid){
    Vector3D normal = surface / surface.magnitude();
    Vector3D velocity = {primitive[1], primitive[2], primitive[3]};
    FloatType normalVelocity = velocity.dot(normal);
    FloatType pressure = fluid.computePressure_rho_u_et(primitive[0], velocity, primitive[4]);
    FloatType totEnthalpy = fluid.computeTotalEnthalpy_rho_u_et(primitive[0], velocity, primitive[4]);

    StateVector flux {};
    flux[0] = primitive[0] * normalVelocity;
    flux[1] = primitive[0] * normalVelocity * velocity.x() + pressure * normal.x();
    flux[2] = primitive[0] * normalVelocity * velocity.y() + pressure * normal.y();
    flux[3] = primitive[0] * normalVelocity * velocity.z() + pressure * normal.z();
    flux[4] = primitive[0] * normalVelocity * totEnthalpy;
    return flux;

}


StateVector computeEulerFluxFromConservative(StateVector conservative, Vector3D surface, CFluid fluid){
    StateVector primitive = getEulerPrimitiveFromConservative(conservative);
    StateVector flux = computeEulerFluxFromPrimitive(primitive, surface, fluid);
    return flux;
}


FloatType computeAngleBetweenVectors(const Vector3D& v1, const Vector3D& v2){
    return acos(v1.dot(v2) / (v1.magnitude() * v2.magnitude()));
}


Vector3D rotateVectorAlongXAxis(const Vector3D& v, FloatType theta){
    Vector3D newVector;
    newVector.x() = v.x();
    newVector.y() = v.y() * cos(theta) - v.z() * sin(theta);
    newVector.z() = v.y() * sin(theta) + v.z() * cos(theta);
    return newVector;
}

Vector3D computeCylindricalVectorFromCartesian(Vector3D vec, FloatType theta){
    Vector3D newVector;
    newVector.x() = vec.x();
    newVector.y() = vec.y() * cos(theta) + vec.z() * sin(theta);
    newVector.z() = - vec.y() * sin(theta) + vec.z() * cos(theta);
    return newVector;
}

Vector3D computeCartesianVectorFromCylindrical(Vector3D vec, FloatType theta){
    Vector3D newVector;
    newVector.x() = vec.x();
    newVector.y() = vec.y() * cos(theta) - vec.z() * sin(theta);
    newVector.z() = vec.y() * sin(theta) + vec.z() * cos(theta);
    return newVector;
}


void integrateRadialEquilibrium(const std::vector<FloatType>& density,
    const std::vector<FloatType>& velocityTang,
    const std::vector<FloatType>& radius,
    const FloatType& hubPressure,
    std::vector<FloatType>& pressure)
{

    const size_t N = radius.size();
    pressure[0] = hubPressure;

    for (size_t i = 1; i < N; ++i) {
        FloatType r1 = radius[i - 1];
        FloatType r2 = radius[i];

        FloatType dr = r2 - r1;

        FloatType rho1 = density[i - 1];
        FloatType rho2 = density[i];

        FloatType ut1 = velocityTang[i - 1];
        FloatType ut2 = velocityTang[i];

        // Evaluate dp/dr = rho * u_theta^2 / r at both points
        FloatType dpdr1 = rho1 * ut1 * ut1 / r1;
        FloatType dpdr2 = rho2 * ut2 * ut2 / r2;

        // Trapezoidal integration step
        pressure[i] = pressure[i - 1] + 0.5 * (dpdr1 + dpdr2) * dr;
    }

}
