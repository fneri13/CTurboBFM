#include "CSourceBFMLiftDrag.hpp"

StateVector CSourceBFMLiftDrag::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    computeFlowState(i, j, k, primitive);
    computeRelativeFlowAngle(i, j, k, primitive);

    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);
    return inviscidComponent + viscousComponent;
}

void CSourceBFMLiftDrag::computeRelativeFlowAngle(size_t i, size_t j, size_t k, const StateVector& primitive) {
    Vector3D relVelMeridional = {_relativeVelocityCylindric.x(), _relativeVelocityCylindric.y(), 0};
    _relativeFlowAngle = std::acos(_relativeVelocityCylindric.dot(relVelMeridional) / _relativeVelocityCylindric.magnitude() / relVelMeridional.magnitude());

    // keep the convenetion used when computing beta_0
    if (_relativeVelocityCylindric.z() < 0) {
        _relativeFlowAngle *= -1;
    }
}


StateVector CSourceBFMLiftDrag::computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce) {
    FloatType solidity = _mesh.getInputFields(FieldNames::SOLIDITY, i, j, k);
    FloatType hParameter = _mesh.getInputFields(FieldNames::LIFT_DRAG_H_PARAMETER, i, j, k);
    FloatType beta0 = _mesh.getInputFields(FieldNames::LIFT_DRAG_BETA_0, i, j, k);


    FloatType forceMag = 2.0 * M_PI * solidity * _relativeVelocityCylindric.dot(_relativeVelocityCylindric) / hParameter;
    
    FloatType omega = _mesh.getInputFields(FieldNames::RPM, i, j, k) * 2 * M_PI / 60;
    if (omega < 0) {
        forceMag *= (_relativeFlowAngle - beta0);
    }
    else {
        forceMag *= (beta0 - _relativeFlowAngle);
    }
    
    Vector3D forceCylindrical = _inviscidForceDirectionCylindrical * forceMag;
    Vector3D forceCartesian = computeCartesianVectorFromCylindrical(forceCylindrical, _theta);
    inviscidForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}

StateVector CSourceBFMLiftDrag::computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce) {
    FloatType kp_etaMax = _mesh.getInputFields(FieldNames::LIFT_DRAG_KP_ETA_MAX, i, j, k);
    FloatType beta_etaMax = _mesh.getInputFields(FieldNames::LIFT_DRAG_BETA_ETA_MAX, i, j, k);
    FloatType solidity = _mesh.getInputFields(FieldNames::SOLIDITY, i, j, k);
    FloatType hParameter = _mesh.getInputFields(FieldNames::LIFT_DRAG_H_PARAMETER, i, j, k);
    
    FloatType kp = kp_etaMax + 2.0 * M_PI * solidity * (_relativeFlowAngle - beta_etaMax) * (_relativeFlowAngle - beta_etaMax);
        
    FloatType forceMag = _relativeVelocityCylindric.dot(_relativeVelocityCylindric) * kp / hParameter;

    Vector3D forceCylindrical = _viscousForceDirectionCylindrical * forceMag;
    Vector3D forceCartesian = computeCartesianVectorFromCylindrical(forceCylindrical, _theta);
    viscousForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}