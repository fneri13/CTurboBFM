#include "CSourceBFMGong.hpp"

StateVector CSourceBFMGong::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    computeFlowState(i, j, k, primitive);
    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);
    return inviscidComponent + viscousComponent;
}


StateVector CSourceBFMGong::computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce) {
    FloatType Kn = computeKnCoefficient(_normalCamberAxial, _normalCamberRadial, _normalCamberTangential);
    
    FloatType forceMag = Kn/_pitch * _relativeVelocityCylindric.dot(_relativeVelocityCylindric) * 0.5 * std::sin(2 * _deviationAngle);
    
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

FloatType CSourceBFMGong::computeKnCoefficient(FloatType normalCamberAxial, FloatType normalCamberRadial, FloatType normalCamberTangential) const {
   FloatType Kn = 0;

   // No idea here how the metal angle should be with sign or no, positive or no, depending on rotation direction or no. Convention not clear
   FloatType metalAngle = std::atan2(std::sqrt(normalCamberAxial * normalCamberAxial + normalCamberRadial * normalCamberRadial), normalCamberTangential);
   Kn = _Kn1 + _Kn2 * metalAngle;
   return Kn;
}

StateVector CSourceBFMGong::computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce) {
    FloatType forceMag = _Kp * _relativeVelocityCylindric.dot(_relativeVelocityCylindric) / _pitch;
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