#include "CSourceBFMHall.hpp"

StateVector CSourceBFMHall::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    computeFlowState(i, j, k, primitive);
    
    FloatType volume = _mesh.getVolume(i, j, k);

    // when the deviation becomes negative, the force mag becomes negative the pull the flow back into place aligned with the camber
    FloatType forceMagnitude = _relativeVelocityCylindric.dot(_relativeVelocityCylindric) * M_PI * _deviationAngle / _pitch / std::abs(_normalCamberTangential);
    
    Vector3D forceCylindrical = _inviscidForceDirectionCylindrical * forceMagnitude;
    Vector3D forceCartesian = computeCartesianVectorFromCylindrical(forceCylindrical, _theta);
    
    inviscidForce(i, j, k) = forceCartesian;
    viscousForce(i, j, k) = Vector3D(0,0,0);
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;

    return source*volume*primitive[0];
}