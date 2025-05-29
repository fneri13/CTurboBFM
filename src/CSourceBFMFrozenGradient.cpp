#include "CSourceBFMFrozenGradient.hpp"

StateVector CSourceBFMFrozenGradient::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    computeFlowState(i, j, k, primitive);
    updateState(i, j, k, primitive);

    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);
    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);

    return inviscidComponent + viscousComponent;
}


StateVector CSourceBFMFrozenGradient::computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce) {

    inviscidForce(i, j, k) = _inviscidForceCartesian*_bladeIsPresent;
    
    StateVector source({0,0,0,0,0});
    source[1] = _inviscidForceCartesian.x();
    source[2] = _inviscidForceCartesian.y();
    source[3] = _inviscidForceCartesian.z();
    source[4] = _inviscidForceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0]*_bladeIsPresent;
}

StateVector CSourceBFMFrozenGradient::computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce) {

    viscousForce(i, j, k) = _viscousForceCartesian*_bladeIsPresent;

    StateVector source({0,0,0,0,0});
    source[1] = _viscousForceCartesian.x();
    source[2] = _viscousForceCartesian.y();
    source[3] = _viscousForceCartesian.z();
    source[4] = _viscousForceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0]*_bladeIsPresent;
}

void CSourceBFMFrozenGradient::updateState(size_t i, size_t j, size_t k, const StateVector& primitive) {
    FloatType angularMomentumDerivative = _mesh.getInputFields(FieldNames::ANGULAR_MOMENTUM_DERIVATIVE, i, j, k);
    FloatType entropyDerivative = _mesh.getInputFields(FieldNames::ENTROPY_DERIVATIVE, i, j, k);

    // compute the theta global component
    FloatType meridionalVelMag = std::sqrt(_relativeVelocityCylindric.y() * _relativeVelocityCylindric.y() + _relativeVelocityCylindric.x() * _relativeVelocityCylindric.x());
    FloatType ftheta = angularMomentumDerivative * meridionalVelMag / _radius;

    // compute the loss component
    FloatType relVelMag = _relativeVelocityCylindric.magnitude();
    FloatType temperature = _fluid.computeTemperature_rho_u_et(primitive[0], {primitive[1], primitive[2], primitive[3]}, primitive[4]);
    FloatType floss = temperature * meridionalVelMag / relVelMag * entropyDerivative;

    // reconstruct the viscous force
    _viscousForceCartesian = _viscousForceDirectionCartesian * floss;
    _viscousForceCylindrical = computeCylindricalVectorFromCartesian(_viscousForceCartesian, _theta);
    
    // reconstruct the turning force
    FloatType fnTheta = ftheta - _viscousForceCylindrical.z();
    FloatType fnMag = std::abs(fnTheta) / std::abs(_inviscidForceDirectionCylindrical.z());
    _inviscidForceCylindrical = _inviscidForceDirectionCylindrical * fnMag;
    _inviscidForceCartesian = computeCartesianVectorFromCylindrical(_inviscidForceCylindrical, _theta);
    
}