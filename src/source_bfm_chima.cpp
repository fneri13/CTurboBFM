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
    _trailingEdgeIndex(config.getTrailingEdgeIndex()) {}


StateVector SourceBFMChima::computeBodyForceSource(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &inviscidForce, 
    Matrix3D<Vector3D> &viscousForce, 
    FlowSolution &conservativeVars) {
    
    computeFlowState(i, j, k, primitive, conservativeVars);

    FloatType currentMassFlow = _turboPerformance[TurboPerformance::MASS_FLOW].back();
    _scalingTurning = linearInterpolation(
        _inputTable.getField(InputField::CHIMA_MASS_FLOW), 
        _inputTable.getField(InputField::CHIMA_SCALING_TURNING), 
        currentMassFlow);
    _scalingLoss = linearInterpolation(
        _inputTable.getField(InputField::CHIMA_MASS_FLOW), 
        _inputTable.getField(InputField::CHIMA_SCALING_LOSS), 
        currentMassFlow);
    
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);
    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    
    return inviscidComponent + viscousComponent;
}


StateVector SourceBFMChima::computeInviscidComponent(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &inviscidForce) {
    
    FloatType deltaTotEnthalpy_dm_ref = _mesh.getInputFields(InputField::DELTA_TOT_ENTHALPY_DM, i, j, k);
    FloatType deltaTotEnthalpy_dm = deltaTotEnthalpy_dm_ref * _scalingTurning;
    _tangentialForce = deltaTotEnthalpy_dm * _velMeridional /( _radius * _omega);
    FloatType forceInviscidTangential = _tangentialForce - _viscousForceCylindrical.z();
    FloatType forceInviscidMagnitude = std::abs(forceInviscidTangential / _inviscidForceDirCylindrical.z());

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




StateVector SourceBFMChima::computeViscousComponent(
    size_t i, 
    size_t j, 
    size_t k, 
    const StateVector& primitive, 
    Matrix3D<Vector3D> &viscousForce) {
    
    FloatType deltaS_deltaM_ref = _mesh.getInputFields(InputField::DELTA_ENTROPY_DM, i, j, k);
    FloatType deltaS_deltaM = deltaS_deltaM_ref * _scalingLoss;
    FloatType relVelMag = _relVelCylindric.magnitude();
    FloatType temperature = _fluid.computeTemperature_rho_u_et(
        primitive[0], 
        {primitive[1], primitive[2], primitive[3]}, 
        primitive[4]);
    FloatType forceMag = temperature * _velMeridional / relVelMag * deltaS_deltaM;

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


