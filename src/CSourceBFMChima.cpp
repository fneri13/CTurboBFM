#include "CSourceBFMChima.hpp"

CSourceBFMChima::CSourceBFMChima(const Config &config, const CFluid &fluid, const CMesh &mesh,
                                  std::map<TurboPerformance, std::vector<FloatType>> &turboPerformance) 
    : CSourceBFMBase(config, fluid, mesh),
      _inputTable(config.getChimaScalingFunctionsFile()), _turboPerformance(turboPerformance)
{}


StateVector CSourceBFMChima::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    computeFlowState(i, j, k, primitive);

    FloatType currentMassFlow = _turboPerformance[TurboPerformance::MASS_FLOW].back();
    
    _scalingTurning = interpolateLinear(_inputTable.getField(FieldNames::CHIMA_MASS_FLOW), _inputTable.getField(FieldNames::CHIMA_SCALING_TURNING), currentMassFlow);
    _scalingLoss = interpolateLinear(_inputTable.getField(FieldNames::CHIMA_MASS_FLOW), _inputTable.getField(FieldNames::CHIMA_SCALING_LOSS), currentMassFlow);
    
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);

    computeGlobalTangentialForce(i, j, k, primitive);

    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    
    return inviscidComponent + viscousComponent;
}


StateVector CSourceBFMChima::computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce) {
    FloatType fnTan = _tangentialForce - _viscousForceCylindrical.z();
    fnTan *= _scalingTurning;
    FloatType fnMag = fnTan / std::abs(_inviscidForceDirectionCylindrical.z());
    // fnMag *= _bladeIsPresent;

    Vector3D forceCylindrical = _inviscidForceDirectionCylindrical * fnMag;
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


void CSourceBFMChima::computeGlobalTangentialForce(size_t i, size_t j, size_t k, const StateVector& primitive) {
    FloatType deltaAngularMomentum_dm = _mesh.getInputFields(FieldNames::DELTA_ANGULAR_MOMENTUM_DM, i, j, k);
    _tangentialForce = deltaAngularMomentum_dm * primitive[0] * _velMeridional / _radius;
}


StateVector CSourceBFMChima::computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce) {

    FloatType deltaS_deltaM = _mesh.getInputFields(FieldNames::DELTA_ENTROPY_DM, i, j, k);
    
    FloatType relVelMag = _relativeVelocityCylindric.magnitude();
    
    _velMeridional = std::sqrt(_relativeVelocityCylindric.x() * _relativeVelocityCylindric.x() + 
                               _relativeVelocityCylindric.y() * _relativeVelocityCylindric.y());

    FloatType temperature = _fluid.computeTemperature_rho_u_et(primitive[0], {primitive[1], primitive[2], primitive[3]}, primitive[4]);
    
    FloatType forceMag = temperature * _velMeridional / relVelMag * deltaS_deltaM;

    // forceMag *= _bladeIsPresent;
    forceMag *= _scalingLoss;
    
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


