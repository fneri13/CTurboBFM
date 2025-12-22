#include "CSourceBFMCorrelations.hpp"

StateVector CSourceBFMCorrelations::computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive, 
            Matrix3D<Vector3D> &inviscidForce, Matrix3D<Vector3D> &viscousForce, FlowSolution &conservativeVars) {

    computeFlowState(i, j, k, primitive, conservativeVars);


    StateVector inviscidComponent = computeInviscidComponent(i, j, k, primitive, inviscidForce);
    StateVector viscousComponent = computeViscousComponent(i, j, k, primitive, viscousForce);
    return inviscidComponent + viscousComponent;
}


StateVector CSourceBFMCorrelations::computeInviscidComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &inviscidForce) {
    
    // formulation found in Aungier book, corrected with Cetin compressibility effects
    FloatType kappa1Deg = _mesh.getInputFields(FieldNames::BLADE_METAL_ANGLE, _leadingEdgeIdx, j, k) * 180 / M_PI; // with sign
    FloatType kappa2Deg = _mesh.getInputFields(FieldNames::BLADE_METAL_ANGLE, _trailingEdgeIdx, j, k) * 180 / M_PI; // with sign
    
    FloatType incidenceInletDeg = _flowAngleInletDeg - kappa1Deg; // with sign

    // Keep in mind that correlations formula and values are all supposed to have positive values -> change them later according to rotation direction
    FloatType _flowAngleInletDegAbs = std::abs(_flowAngleInletDeg);
    FloatType n = 0.025*_solidity - 0.06 - (std::pow(_flowAngleInletDegAbs/90.0, 1.0+1.2*_solidity)) / (1.5 + 0.43*_solidity);
    FloatType p = 0.914 + _solidity*_solidity*_solidity / 160.0;
    _incidenceZeroTenThkStar = std::pow(_flowAngleInletDegAbs,p) / (5.0 + 46.0 * std::exp(-2.3*_solidity)) - 0.1*std::pow(_solidity, 3.0) * std::exp((_flowAngleInletDegAbs-70.0) / 4.0);

    // find thickness max along the stream
    FloatType thkMax = 0.0;
    FloatType thkTmp = 0.0;
    for (size_t i = _leadingEdgeIdx; i < _trailingEdgeIdx; i++) {
        thkTmp = _pitch * (1.0 - _blockage);
        if (thkTmp > thkMax) {
            thkMax = thkTmp;
        }
    }
    
    FloatType tcRatio = thkMax / _streamBladeLength;
    FloatType Kti = std::pow(10.0*tcRatio, 0.28/(0.1 + std::pow(tcRatio,0.3)));
    FloatType Ksh = 0.8;// shape factor, for now 1, i don't know which one is it

    FloatType istar = Ksh*Kti*_incidenceZeroTenThkStar + n*_camberAngleDegAbs;
    
    // correction by Cetin
    FloatType inletMach = _relativeVelocityCartesian.magnitude() / _fluid.computeSoundSpeed_rho_u_et(primitive[0], _relativeVelocityCartesian, primitive[4]);
    // FloatType istarCorrected = istar; // no correction
    FloatType istarCorrected = istar + 1.3026*inletMach + 5.7380;

    // compute deviation now
    FloatType deviation10OutStar = 0.01*_solidity*_flowAngleInletDegAbs + (0.74*std::pow(_solidity,1.9) + 3.0*_solidity) * std::pow(_flowAngleInletDegAbs/90.0, 1.67+1.09*_solidity);
    FloatType Ktd = 6.25*tcRatio + 37.5*tcRatio*tcRatio;
    FloatType acRatio = 0.5; // max thickness location / chord length
    FloatType deviationStar = (0.92*acRatio*acRatio + 0.002*std::abs(kappa2Deg))/(1-0.002*_camberAngleDegAbs/std::sqrt(_solidity)) * _camberAngleDegAbs/std::sqrt(_solidity) + (Ksh*Ktd-1.0)*deviation10OutStar;
    
    // correction by Cetin
    // FloatType deviationStarCorrected = deviationStar + 0.0*inletMach + 0.0; // no correction
    FloatType deviationStarCorrected = -1.099379 + 3.0186*deviationStar - 0.1988*deviationStar*deviationStar;

    // off-design by aungier
    FloatType dd_distar = (1.0 + (_solidity+0.25*std::pow(_solidity, 4.0)) * std::pow((_flowAngleInletDegAbs/53.0), 2.5)) / std::exp(3.1*_solidity);
    FloatType umInlet = std::sqrt(_relativeVelocityInlet.x()*_relativeVelocityInlet.x() + _relativeVelocityInlet.y()*_relativeVelocityInlet.y());
    FloatType umOutlet = std::sqrt(_relativeVelocityOutlet.x()*_relativeVelocityOutlet.x() + _relativeVelocityOutlet.y()*_relativeVelocityOutlet.y());

    FloatType incidenceInletDegMag = 0.0;
    if (_omega>0.0) {
        incidenceInletDegMag = incidenceInletDeg*(-1.0); // incidence positive when positive flow turning needs to be applied
    }
    else {
        incidenceInletDegMag = incidenceInletDeg;
    }
    FloatType deviationOutletDeg = deviationStarCorrected + dd_distar*(incidenceInletDegMag - istarCorrected) + 10.0*(1.0 - umOutlet/umInlet);

    // curvature radius
    FloatType curvatureRadius = std::abs(1.0/_mesh.getInputFields(FieldNames::D_BLADE_METAL_ANGLE_DM, i, j, k)); // this always positive i would say
    FloatType meridionalNormalized = (_mesh.getInputFields(FieldNames::STREAMWISE_LENGTH, i,j,k) - _mesh.getInputFields(FieldNames::STREAMWISE_LENGTH, _leadingEdgeIdx,j,k)) /(
                                      _mesh.getInputFields(FieldNames::STREAMWISE_LENGTH, _trailingEdgeIdx,j,k) - _mesh.getInputFields(FieldNames::STREAMWISE_LENGTH, _leadingEdgeIdx,j,k));
    
    
    FloatType deviationExpectedDeg = incidenceInletDegMag + meridionalNormalized * (deviationOutletDeg - incidenceInletDegMag);
    FloatType deviationExpectedRad = deviationExpectedDeg * M_PI / 180.0;

    FloatType uThetaCurrent = _relativeVelocityCylindric.z();
    FloatType kappa = _mesh.getInputFields(FieldNames::BLADE_METAL_ANGLE, i, j, k); // this is with sign
    FloatType Wn, uThetaExpected;
    
    // think that Wn decides the push-pull effect of the blade
    if (_omega>0.0) {
        uThetaExpected = _velMeridional * std::tan(-deviationExpectedRad + kappa);
        Wn = uThetaExpected - uThetaCurrent; // error term in the force model
    }
    else {
        uThetaExpected = _velMeridional * std::tan(deviationExpectedRad + kappa);
        Wn = -uThetaExpected + uThetaCurrent;
    }
    
    FloatType KN = _config.getKnCorrelationBfmCoefficient(); // coeff to change for feeeback control
    FloatType fnMag = KN * primitive[0] * _velMeridional * Wn / _pitch  + primitive[0] * _relativeVelocityCylindric.magnitude() * _relativeVelocityCylindric.magnitude() / curvatureRadius;
    

    if (std::isnan(fnMag)) {
        FloatType a = 0.0;   
    }
    else if (fnMag < 0.0) {
        FloatType a = 0.0;   
    }
    else {
        FloatType a = 0.0;   
    }

    // compute forces now
    Vector3D forceCartesian = _inviscidForceDirectionCartesian * fnMag;
    Vector3D forceCylindrical = _inviscidForceDirectionCylindrical * fnMag;

    // update the force in memory
    inviscidForce(i, j, k) = forceCartesian;
    
    // compute the gov equations source
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume;
}


StateVector CSourceBFMCorrelations::computeViscousComponent(size_t i, size_t j, size_t k, const StateVector& primitive, Matrix3D<Vector3D> &viscousForce) {
    Vector3D forceCartesian = {0,0,0};
    Vector3D forceCylindrical = {0,0,0};
    
    // update the force in memory
    viscousForce(i, j, k) = forceCartesian;
    
    StateVector source({0,0,0,0,0});
    source[1] = forceCartesian.x();
    source[2] = forceCartesian.y();
    source[3] = forceCartesian.z();
    source[4] = forceCylindrical.z() * _omega * _radius;
    
    FloatType volume = _mesh.getVolume(i, j, k);

    return source*volume*primitive[0];
}
