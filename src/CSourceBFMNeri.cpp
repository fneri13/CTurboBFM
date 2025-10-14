#include "CSourceBFMNeri.hpp"

CSourceBFMNeri::CSourceBFMNeri( const Config &config, const CFluidBase &fluid, const CMesh &mesh, FlowSolution &conservativeVars, 
                                std::map<TurboPerformance, std::vector<FloatType>> &turboPerformance) 
                                : CSourceBFMBase(config, fluid, mesh), _conservativeSolution(conservativeVars), _inputTable(config.getChimaScalingFunctionsFile()), 
                                _turboPerformance(turboPerformance), _trailingEdgeIndex(config.getTrailingEdgeIndex()), _leadingEdgeIndex(config.getLeadingEdgeIndex()) {}



StateVector CSourceBFMNeri::computeBodyForceSource( size_t i, size_t j, size_t k, const StateVector& primitive,
                                                    Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce) {
    computeFlowState(i, j, k, primitive);

    computeStreamwiseCoefficient(i, j, k);
                                                    
    FloatType currentMassFlow = _turboPerformance[TurboPerformance::MASS_FLOW].back();
    
    _scalingTurning = interpolateLinear(_inputTable.getField(FieldNames::CHIMA_MASS_FLOW), _inputTable.getField(FieldNames::CHIMA_SCALING_TURNING), currentMassFlow);
    
    _scalingLoss = interpolateLinear(_inputTable.getField(FieldNames::CHIMA_MASS_FLOW), _inputTable.getField(FieldNames::CHIMA_SCALING_LOSS), currentMassFlow);
    
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);

    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    
    return inviscidComponent + viscousComponent;
}



void CSourceBFMNeri::computeStreamwiseCoefficient(size_t i, size_t j, size_t k) {
    
    FloatType mLocal{0.0};

    mLocal = _mesh.getInputFields(FieldNames::STREAMWISE_LENGTH, i, j, k) + 1e-8;

    _mTrailingEdge = _mesh.getInputFields(FieldNames::STREAMWISE_LENGTH, _trailingEdgeIndex, j, k);

    _mLeadingEdge = _mesh.getInputFields(FieldNames::STREAMWISE_LENGTH, _leadingEdgeIndex, j, k);

    _streamwiseCoeff = std::sqrt((mLocal-_mLeadingEdge)/(_mTrailingEdge-_mLeadingEdge)) * 0.001; // 0.001 for r4 works
}




StateVector CSourceBFMNeri::computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce) {
    
    FloatType deltaS_deltaM_ref = _mesh.getInputFields(FieldNames::DELTA_ENTROPY_DM, i, j, k);
    
    deltaS_deltaM_ref *= _scalingLoss;
    
    FloatType relVelMag = _relativeVelocityCylindric.magnitude();

    FloatType temperature = _fluid.computeTemperature_rho_u_et(primitive[0], {primitive[1], primitive[2], primitive[3]}, primitive[4]);

    FloatType forceMagOld = viscousForce(i, j, k).magnitude();

    StateVector conservativeTrailEdge = _conservativeSolution.at(_trailingEdgeIndex, j, k);

    StateVector primitiveTrailEdge = getEulerPrimitiveFromConservative(conservativeTrailEdge);

    FloatType entropyTrailEdge = _fluid.computeEntropy_rho_u_et(primitiveTrailEdge[0], {primitiveTrailEdge[1], primitiveTrailEdge[2], primitiveTrailEdge[3]}, primitiveTrailEdge[4]);
    
    FloatType deltaS_deltaM_actual = entropyTrailEdge / (_mTrailingEdge - _mLeadingEdge);    
    
    FloatType forceMag = forceMagOld - (deltaS_deltaM_actual - deltaS_deltaM_ref) * temperature * _velMeridional / relVelMag * _streamwiseCoeff;

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
    
    FloatType deltaTotEnthalpy_dm_ref = _mesh.getInputFields(FieldNames::DELTA_TOT_ENTHALPY_DM, i, j, k);
    
    deltaTotEnthalpy_dm_ref *= _scalingTurning;

    StateVector conservativeTrailEdge = _conservativeSolution.at(_trailingEdgeIndex, j, k);

    StateVector primitiveTrailEdge = getEulerPrimitiveFromConservative(conservativeTrailEdge);

    FloatType totEnthalpyTrailEdge = _fluid.computeTotalEnthalpy_rho_u_et(primitiveTrailEdge[0], {primitiveTrailEdge[1], primitiveTrailEdge[2], primitiveTrailEdge[3]}, primitiveTrailEdge[4]); 

    FloatType deltaTotEnthalpy_dm_actual = (totEnthalpyTrailEdge - 1005.0*288.15) / (_mTrailingEdge - _mLeadingEdge);

    FloatType tangForceOld = inviscidForce(i, j, k).z() + _viscousForceCylindrical.z();
    
    _tangentialForce = tangForceOld - (deltaTotEnthalpy_dm_actual - deltaTotEnthalpy_dm_ref) * _velMeridional / _omega / _radius * _streamwiseCoeff;

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
