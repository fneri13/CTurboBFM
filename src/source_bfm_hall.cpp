#include "source_bfm_hall.hpp"

StateVector SourceBFMHall::computeBodyForceSource(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &inviscidForce, 
    Matrix3D<Vector3D> &viscousForce, 
    FlowSolution &conservativeVars,
    FloatType &dt) {

    computeFlowState(i, j, k, primitive, conservativeVars);
    
    FloatType volume = _mesh.getVolume(i, j, k);

    FloatType forceMagnitude = _relVelCylindric.dot(_relVelCylindric) * M_PI * _deviationAngle / 
        _tangentialPitch / std::abs(_normalCamberTangential);
    
    Vector3D forceCylindrical = _inviscidForceDirCylindrical * forceMagnitude;
    Vector3D forceCartesian = computeCartesianComponentsFromCylindrical(forceCylindrical, _theta);
    
    inviscidForce(i, j, k) = forceCartesian;
    viscousForce(i, j, k) = Vector3D(0,0,0);
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;

    return source*volume*primitive[0];
}