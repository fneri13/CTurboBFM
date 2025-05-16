#include "CSourceBFMFrozenForce.hpp"

StateVector CSourceBFMFrozenForce::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    computeFlowState(i, j, k, primitive);
    updateState(i, j, k, primitive);

    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);
    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);

    return inviscidComponent + viscousComponent;
}


StateVector CSourceBFMFrozenForce::computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce) {
    Vector3D forceCartesian = _forceCartesian - _viscousForceCartesian;
    Vector3D forceCyl = computeCylindricalVectorFromCartesian(forceCartesian, _theta);

    inviscidForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCyl.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}

StateVector CSourceBFMFrozenForce::computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce) {

    viscousForce(i, j, k) = _viscousForceCartesian;

    StateVector source({0,0,0,0,0});
    source[1] = _viscousForceCartesian.x();
    source[2] = _viscousForceCartesian.y();
    source[3] = _viscousForceCartesian.z();
    source[4] = _viscousForceCyl.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}

void CSourceBFMFrozenForce::updateState(size_t i, size_t j, size_t k, const StateVector& primitive) {
    _forceCylindrical = {_mesh.getInputFields(FieldNames::AXIAL_FORCE, i, j, k), _mesh.getInputFields(FieldNames::RADIAL_FORCE, i, j, k), _mesh.getInputFields(FieldNames::TANGENTIAL_FORCE, i, j, k)};
    _forceCartesian = computeCartesianVectorFromCylindrical(_forceCylindrical, _theta);
    _viscousForceDirectionCartesian = - _relativeVelocityCartesian.normalized();
    _viscousForceCartesian = _viscousForceDirectionCartesian *_forceCartesian.dot(_viscousForceDirectionCartesian);
    _viscousForceCyl = computeCylindricalVectorFromCartesian(_viscousForceCartesian, _theta);
}