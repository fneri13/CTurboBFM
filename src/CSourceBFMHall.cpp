#include "CSourceBFMHall.hpp"

StateVector CSourceBFMHall::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive) {
    computeFlowState(i, j, k, primitive);
    
    FloatType volume = _mesh.getVolume(i, j, k);
    FloatType forceMagnitude = _relativeVelocityCylindric.dot(_relativeVelocityCylindric) * M_PI * _deviationAngle / _pitch / std::abs(_normalCamberTangential);
    Vector3D forceCylindrical = _inviscidForceDirection * forceMagnitude;
    Vector3D forceCartesian = computeCartesianVectorFromCylindrical(forceCylindrical, _theta);
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;

    return source*volume;
}