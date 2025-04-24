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
    conservative[1] = primitive[1] * conservative[0];
    conservative[2] = primitive[2] * conservative[0];
    conservative[3] = primitive[3] * conservative[0];
    conservative[4] = primitive[4] * conservative[0];
    return conservative;
}

StateVector computeEulerFluxFromPrimitive(StateVector state, Vector3D surface, CFluid fluid){
    Vector3D normal = surface / surface.magnitude();
    Vector3D velocity = {state[1], state[2], state[3]};
    FloatType normalVelocity = velocity.dot(normal);
    FloatType pressure = fluid.computePressure_rho_u_et(state[0], velocity, state[4]);
    FloatType totEnthalpy = fluid.computeTotalEnthalpy_rho_u_et(state[0], velocity, state[4]);

    StateVector flux {};
    flux[0] = state[0] * normalVelocity;
    flux[1] = state[0] * normalVelocity * velocity.x() + pressure * normal.x();
    flux[2] = state[0] * normalVelocity * velocity.y() + pressure * normal.y();
    flux[3] = state[0] * normalVelocity * velocity.z() + pressure * normal.z();
    flux[4] = state[0] * normalVelocity * totEnthalpy;
    return flux;

}