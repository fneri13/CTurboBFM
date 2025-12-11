#include "CSourceBFMNeri.hpp"

CSourceBFMNeri::CSourceBFMNeri( const Config &config, const CFluidBase &fluid, const CMesh &mesh, FlowSolution &conservativeVars, 
                                std::map<TurboPerformance, std::vector<FloatType>> &turboPerformance) 
                                : CSourceBFMBase(config, fluid, mesh), _conservativeSolution(conservativeVars), _inputTable(config.getChimaScalingFunctionsFile()), 
                                _turboPerformance(turboPerformance), _trailingEdgeIndex(config.getTrailingEdgeIndex()), _leadingEdgeIndex(config.getLeadingEdgeIndex()) {}



StateVector CSourceBFMNeri::computeBodyForceSource( size_t i, size_t j, size_t k, const StateVector& primitive,
                                                    Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    // update flow state
    computeFlowState(i, j, k, primitive);
    computeStreamwiseCoefficient(i, j, k);
    
    // compute scaling coefficients for turning and loss profiles
    FloatType currentMassFlow = _turboPerformance[TurboPerformance::MASS_FLOW].back();
    _scalingTurning = interpolateLinear(_inputTable.getField(FieldNames::CHIMA_MASS_FLOW), _inputTable.getField(FieldNames::CHIMA_SCALING_TURNING), currentMassFlow);
    _scalingLoss = interpolateLinear(_inputTable.getField(FieldNames::CHIMA_MASS_FLOW), _inputTable.getField(FieldNames::CHIMA_SCALING_LOSS), currentMassFlow);
    
    // compute body force contributions
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);
    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    
    return inviscidComponent + viscousComponent;
}



void CSourceBFMNeri::computeStreamwiseCoefficient(size_t i, size_t j, size_t k) {
    // streamwise distribution function taken from Chima article 2006
    FloatType mLocal{0.0};
    mLocal = _mesh.getInputFields(FieldNames::STREAMWISE_LENGTH, i, j, k) + 1e-9;
    _mTrailingEdge = _mesh.getInputFields(FieldNames::STREAMWISE_LENGTH, _trailingEdgeIndex, j, k);
    _mLeadingEdge = _mesh.getInputFields(FieldNames::STREAMWISE_LENGTH, _leadingEdgeIndex, j, k);
    _streamwiseCoeff = std::sqrt((mLocal-_mLeadingEdge)/(_mTrailingEdge-_mLeadingEdge)) * 0.001; // 0.001 for r4 works
}




StateVector CSourceBFMNeri::computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce) {
    
    // compute the magnitude
    FloatType deltaS_deltaM_ref = _mesh.getInputFields(FieldNames::DELTA_ENTROPY_DM, i, j, k);
    FloatType deltaS_deltaM = deltaS_deltaM_ref * _scalingLoss;
    // FloatType relVelMag = _relativeVelocityCylindric.magnitude();
    // FloatType temperature = _fluid.computeTemperature_rho_u_et(primitive[0], {primitive[1], primitive[2], primitive[3]}, primitive[4]);
    FloatType forceMagOld = viscousForce(i, j, k).magnitude();
    StateVector conservativeAtTrailEdge = _conservativeSolution.at(_trailingEdgeIndex, j, k);
    StateVector primitiveAtTrailEdge = getEulerPrimitiveFromConservative(conservativeAtTrailEdge);
    FloatType entropyAtTrailEdge = _fluid.computeEntropy_rho_u_et(primitiveAtTrailEdge[0], {primitiveAtTrailEdge[1], primitiveAtTrailEdge[2], primitiveAtTrailEdge[3]}, primitiveAtTrailEdge[4]);
    FloatType deltaS_deltaM_actual = entropyAtTrailEdge / (_mTrailingEdge - _mLeadingEdge);    

    // compute the new magnitude, assuming proportional control over deltaS/deltaM
    FloatType forceMag = forceMagOld - (deltaS_deltaM_actual - deltaS_deltaM) / deltaS_deltaM * _streamwiseCoeff;
    
    // contribution to bfm gov. equations
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




StateVector CSourceBFMNeri::computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce) {
    
    // update the state
    FloatType deltaTotEnthalpy_dm_ref = _mesh.getInputFields(FieldNames::DELTA_TOT_ENTHALPY_DM, i, j, k);
    FloatType deltaTotEnthalpy_dm = deltaTotEnthalpy_dm_ref * _scalingTurning;
    StateVector conservativeTrailEdge = _conservativeSolution.at(_trailingEdgeIndex, j, k);
    StateVector primitiveTrailEdge = getEulerPrimitiveFromConservative(conservativeTrailEdge);
    FloatType totEnthalpyTrailEdge = _fluid.computeTotalEnthalpy_rho_u_et(primitiveTrailEdge[0], {primitiveTrailEdge[1], primitiveTrailEdge[2], primitiveTrailEdge[3]}, primitiveTrailEdge[4]); 
    FloatType deltaTotEnthalpy_dm_actual = (totEnthalpyTrailEdge - 1005.0*288.15) / (_mTrailingEdge - _mLeadingEdge);

    // compute the magnitude
    FloatType tangForceOld = inviscidForce(i, j, k).z() + _viscousForceCylindrical.z();
    _tangentialForce = tangForceOld - (deltaTotEnthalpy_dm_actual - deltaTotEnthalpy_dm) / deltaTotEnthalpy_dm * _streamwiseCoeff;
    FloatType forceInviscidTangential = _tangentialForce - _viscousForceCylindrical.z();
    FloatType forceInviscidMagnitude = std::abs(forceInviscidTangential / _inviscidForceDirectionCylindrical.z());

    // contribution to flow equations
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
