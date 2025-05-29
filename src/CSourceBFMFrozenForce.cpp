#include "CSourceBFMFrozenForce.hpp"

StateVector CSourceBFMFrozenForce::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    computeFlowState(i, j, k, primitive);
    updateState(i, j, k, primitive);

    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);
    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);

    return inviscidComponent + viscousComponent;
}


StateVector CSourceBFMFrozenForce::computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce) {

    inviscidForce(i, j, k) = _inviscidForceCartesian*_bladeIsPresent;
    
    StateVector source({0,0,0,0,0});
    source[1] = _inviscidForceCartesian.x();
    source[2] = _inviscidForceCartesian.y();
    source[3] = _inviscidForceCartesian.z();
    source[4] = _inviscidForceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0]*_bladeIsPresent;
}

StateVector CSourceBFMFrozenForce::computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce) {

    viscousForce(i, j, k) = _viscousForceCartesian*_bladeIsPresent;

    StateVector source({0,0,0,0,0});
    source[1] = _viscousForceCartesian.x();
    source[2] = _viscousForceCartesian.y();
    source[3] = _viscousForceCartesian.z();
    source[4] = _viscousForceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0]*_bladeIsPresent;
}

void CSourceBFMFrozenForce::updateState(size_t i, size_t j, size_t k, const StateVector& primitive) {
    Vector3D forceCylindrical = {_mesh.getInputFields(FieldNames::AXIAL_FORCE, i, j, k), _mesh.getInputFields(FieldNames::RADIAL_FORCE, i, j, k), _mesh.getInputFields(FieldNames::TANGENTIAL_FORCE, i, j, k)};
    Vector3D forceCartesian = computeCartesianVectorFromCylindrical(forceCylindrical, _theta);
    _viscousForceCartesian = _viscousForceDirectionCartesian * forceCartesian.dot(_viscousForceDirectionCartesian);
    _viscousForceCylindrical = computeCylindricalVectorFromCartesian(_viscousForceCartesian, _theta);
    _inviscidForceCartesian = forceCartesian - _viscousForceCartesian;
    _inviscidForceCylindrical = computeCylindricalVectorFromCartesian(_inviscidForceCartesian, _theta);
}