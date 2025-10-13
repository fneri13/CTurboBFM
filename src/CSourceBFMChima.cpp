#include "CSourceBFMChima.hpp"

CSourceBFMChima::CSourceBFMChima(const Config &config, const CFluidBase &fluid, const CMesh &mesh,
                                  std::map<TurboPerformance, std::vector<FloatType>> &turboPerformance) 
    : CSourceBFMBase(config, fluid, mesh),
      _inputTable(config.getChimaScalingFunctionsFile()), _turboPerformance(turboPerformance), _trailingEdgeIndex(config.getTrailingEdgeIndex())
{}


StateVector CSourceBFMChima::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    computeFlowState(i, j, k, primitive);

    FloatType currentMassFlow = _turboPerformance[TurboPerformance::MASS_FLOW].back();
    
    _scalingTurning = interpolateLinear(_inputTable.getField(FieldNames::CHIMA_MASS_FLOW), _inputTable.getField(FieldNames::CHIMA_SCALING_TURNING), currentMassFlow);
    
    _scalingLoss = interpolateLinear(_inputTable.getField(FieldNames::CHIMA_MASS_FLOW), _inputTable.getField(FieldNames::CHIMA_SCALING_LOSS), currentMassFlow);
    
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);

    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    
    return inviscidComponent + viscousComponent;
}


StateVector CSourceBFMChima::computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce) {
    
    FloatType deltaTotEnthalpy_dm = _mesh.getInputFields(FieldNames::DELTA_TOT_ENTHALPY_DM, i, j, k);
    
    deltaTotEnthalpy_dm *= _scalingTurning;
    
    _tangentialForce = deltaTotEnthalpy_dm * _velMeridional /( _radius * _omega);

    FloatType forceInviscidTangential = _tangentialForce - _viscousForceCylindrical.z();

    FloatType forceInviscidMagnitude = forceInviscidTangential / _inviscidForceDirectionCylindrical.z();

    Vector3D forceCylindrical = _inviscidForceDirectionCylindrical * forceInviscidMagnitude;

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




StateVector CSourceBFMChima::computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce) {
    FloatType deltaS_deltaM_ref = _mesh.getInputFields(FieldNames::DELTA_ENTROPY_DM, i, j, k);
    
    deltaS_deltaM_ref *= _scalingLoss;
    
    FloatType relVelMag = _relativeVelocityCylindric.magnitude();

    FloatType temperature = _fluid.computeTemperature_rho_u_et(primitive[0], {primitive[1], primitive[2], primitive[3]}, primitive[4]);
    
    FloatType forceMag = temperature * _velMeridional / relVelMag * deltaS_deltaM_ref;

    Vector3D forceCylindrical = _viscousForceDirectionCylindrical * forceMag;
    
    _viscousForceCylindrical = forceCylindrical;

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


