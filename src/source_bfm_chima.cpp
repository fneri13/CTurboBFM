#include "source_bfm_chima.hpp"

SourceBFMChima::SourceBFMChima(
    const Config &config, 
    const FluidBase &fluid, 
    const Mesh &mesh,
    std::map<TurboPerformance, 
    std::vector<FloatType>> &turboPerformance) 
    : SourceBFMBase(config, fluid, mesh),
    _inputTable(config.getChimaScalingFunctionsFile()), 
    _turboPerformance(turboPerformance), 
    _leadingEdgeIndex(_config.getLeadingEdgeIndex()),
    _trailingEdgeIndex(config.getTrailingEdgeIndex()) {}


StateVector SourceBFMChima::computeBodyForceSource(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &inviscidForce, 
    Matrix3D<Vector3D> &viscousForce, 
    FlowSolution &conservativeVars,
    FloatType &dt,
    FloatType &timePhysical) {
    
    computeFlowState(i, j, k, primitive, conservativeVars, timePhysical);

    FloatType currentMassFlow = _turboPerformance[TurboPerformance::MASS_FLOW].back();

    _leadingEdgeIndex = _config.getLeadingEdgeIndex();
    _trailingEdgeIndex = _config.getTrailingEdgeIndex();

    _scalingTurning = linearInterpolation(
        _inputTable.getField(InputField::CHIMA_MASS_FLOW), 
        _inputTable.getField(InputField::CHIMA_SCALING_TURNING), 
        currentMassFlow);

    _scalingLoss = linearInterpolation(
        _inputTable.getField(InputField::CHIMA_MASS_FLOW), 
        _inputTable.getField(InputField::CHIMA_SCALING_LOSS), 
        currentMassFlow);
    
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);
    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce, dt, conservativeVars);
    
    return inviscidComponent + viscousComponent;
}


StateVector SourceBFMChima::computeInviscidComponent(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &inviscidForce,
    FloatType &dt,
    FlowSolution &conservativeVars) {
    
    FloatType deltaL = _mesh.getInputFields(InputField::DELTA_ANGULAR_MOMENTUM, i, j, k);
    FloatType streamCoordLocal = _mesh.getInputFields(InputField::STREAMWISE_LENGTH, i, j, k);
    FloatType streamCoordTrailing = _mesh.getInputFields(InputField::STREAMWISE_LENGTH, _trailingEdgeIndex, j, k);
    FloatType dLdM = computeDerivativeFunction(deltaL, streamCoordLocal, streamCoordTrailing);
    dLdM *= _scalingTurning;

    _tangentialForce = dLdM * _velMeridional / _radius;
    FloatType forceInviscidTangential = _tangentialForce - _viscousForceCylindrical.z();
    // FloatType forceInviscidMagnitude = std::abs(forceInviscidTangential / _inviscidForceDirCylindrical.z());
    FloatType forceInviscidMagnitude = forceInviscidTangential / _inviscidForceDirCylindrical.z();

    Vector3D forceCylindrical = _inviscidForceDirCylindrical * forceInviscidMagnitude;
    Vector3D forceCartesian = computeCartesianComponentsFromCylindrical(forceCylindrical, _theta);
    inviscidForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}


// StateVector SourceBFMChima::computeInviscidComponent(
//     size_t i, 
//     size_t j, 
//     size_t k, 
//     const StateVector& primitive, 
//     Matrix3D<Vector3D> &inviscidForce,
//     FloatType &dt,
//     FlowSolution &conservativeVars) {
    
//     FloatType deltaAngularMomentum = _mesh.getInputFields(InputField::DELTA_ANGULAR_MOMENTUM, i, j, k);
//     FloatType streamCoordLocal = _mesh.getInputFields(InputField::STREAMWISE_LENGTH, i, j, k);
//     FloatType streamCoordTrailing = _mesh.getInputFields(InputField::STREAMWISE_LENGTH, _trailingEdgeIndex, j, k);
//     FloatType dLdM = computeDerivativeFunction(deltaAngularMomentum, streamCoordLocal, streamCoordTrailing);
//     FloatType deltaLdeltaM = deltaAngularMomentum / (streamCoordTrailing);
//     dLdM *= _scalingTurning;
    
//     if (inviscidForce(i, j, k).magnitude() <= 1) {
//         // first iteration, compute global tangential force as turning model
//         _tangentialForce = dLdM * _velMeridional / _radius;
//     }
//     else {
//         // deviation model turning force
//         FloatType deltaL = deltaLdeltaM * streamCoordLocal/streamCoordTrailing;
//         deltaL *= _scalingTurning;
        
//         // get info at leading edge
//         // StateVector conservativeLeadEdge = conservativeVars.at(_leadingEdgeIndex, j, k);
//         // StateVector primitiveLeadEdge = getPrimitiveVariablesFromConservative(conservativeLeadEdge);
//         // Vector3D velLeadEdge = {primitiveLeadEdge[1], primitiveLeadEdge[2], primitiveLeadEdge[3]};
//         // Vector3D velCylLeadEdge = computeCylindricalComponentsFromCartesian(velLeadEdge, _theta);
//         // FloatType velTanLeadEdge = velCylLeadEdge.z();

//         FloatType velThetaRef = deltaL / _radius;
//         FloatType velThetaCurrent = _velCylindrical.z();
//         FloatType varphiTurn = std::sqrt(streamCoordLocal/streamCoordTrailing);
//         FloatType fthetaOld = inviscidForce(i, j, k).z() + _viscousForceCylindrical.z();
//         _tangentialForce = fthetaOld + (velThetaRef - velThetaCurrent) / * varphiTurn;
//     }

//     FloatType forceInviscidTangential = _tangentialForce - _viscousForceCylindrical.z();
//     FloatType forceInviscidMagnitude = forceInviscidTangential / _inviscidForceDirCylindrical.z();
//     Vector3D forceCylindrical = _inviscidForceDirCylindrical * forceInviscidMagnitude;
//     Vector3D forceCartesian = computeCartesianComponentsFromCylindrical(forceCylindrical, _theta);
//     inviscidForce(i, j, k) = forceCartesian;
    
//     StateVector source({0,0,0,0,0});
//     source[1] = forceCartesian.x();
//     source[2] = forceCartesian.y();
//     source[3] = forceCartesian.z();
//     source[4] = forceCylindrical.z() * _omega * _radius;
//     FloatType volume = _mesh.getVolume(i, j, k);

//     return source*volume*primitive[0];
// }


StateVector SourceBFMChima::computeViscousComponent(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &viscousForce) {
    
    FloatType deltaS = _mesh.getInputFields(InputField::DELTA_ENTROPY, i, j, k);
    FloatType streamCoordLocal = _mesh.getInputFields(InputField::STREAMWISE_LENGTH, i, j, k);
    FloatType streamCoordTrailing = _mesh.getInputFields(InputField::STREAMWISE_LENGTH, _trailingEdgeIndex, j, k);
    FloatType dSdM = computeDerivativeFunction(deltaS, streamCoordLocal, streamCoordTrailing);
    dSdM *= _scalingLoss;
    
    FloatType relVelMag = _relVelCylindric.magnitude();
    FloatType temperature = _fluid.computeTemperature_rho_u_et(
        primitive[0], 
        {primitive[1], primitive[2], primitive[3]}, 
        primitive[4]);

    FloatType forceMag = temperature * _velMeridional / relVelMag * dSdM;
    Vector3D forceCylindrical = _viscousForceDirCylindrical * forceMag;
    _viscousForceCylindrical = forceCylindrical;
    Vector3D forceCartesian = computeCartesianComponentsFromCylindrical(forceCylindrical, _theta);
    viscousForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}


FloatType SourceBFMChima::computeDerivativeFunction(
    FloatType yFinal, 
    FloatType xLocal, 
    FloatType xFinal) {
    
    // // parabolic
    FloatType derivative = -2.0*yFinal*xLocal / (xFinal*xFinal) + 2.0*yFinal/xFinal;

    // linear
    // FloatType derivative = yFinal/xFinal;
    
    return derivative;
}
