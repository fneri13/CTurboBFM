#include "CAdvectionSchemeROE.hpp"
#include "commonFunctions.hpp"

StateVector CAdvectionSchemeROE::computeFlux(
    const StateVector& Ull,
    const StateVector& Ul,
    const StateVector& Ur,
    const StateVector& Urr,
    const Vector3D& S)
{

    StateVector Wl = getEulerPrimitiveFromConservative(Ul);
    StateVector Wr = getEulerPrimitiveFromConservative(Ur);

    computeNormalTriad(S);

    StateVector WnormL = computeRotatedPrimitive(Wl);
    StateVector WnormR = computeRotatedPrimitive(Wr);

    computeRoeAVG(WnormL, WnormR);
    computeEigenvalues();
    computeEigenvectors();
    computeWaveStrengths(WnormL, WnormR);

    StateVector flux({0.0, 0.0, 0.0, 0.0, 0.0});

    computeRoeFlux(S, WnormL, WnormR, flux);

    return flux;
}


void CAdvectionSchemeROE::computeNormalTriad(const Vector3D& S){
    
    _x1 = {1.0, 0.0, 0.0}; // default x1 vector
    _x2 = {0.0, 1.0, 0.0}; // default x2 vector
    _x3 = {0.0, 0.0, 1.0}; // default x3 vector

    // the first is the surface normal direction
    _n1 = S.normalized();
    
    // the second is simply taken orthogonal to the first

    if (_n1.x() == 0 && _n1.y() == 0){
        _n2 = Vector3D(1, 0, 0);
    }
    else{
        _n2 = Vector3D(-_n1.y(), _n1.x(), 0);
    }
    _n2 = _n2.normalized();
    
    // the third is the cross product
    _n3 = _n1.cross(_n2);
    
}


StateVector CAdvectionSchemeROE::computeRotatedPrimitive(const StateVector& W) const {
    Vector3D velocity ({W[1], W[2], W[3]});
    StateVector rotatedState ({W[0], velocity.dot(_n1), velocity.dot(_n2), velocity.dot(_n3), W[4]});
    return rotatedState;
}


void CAdvectionSchemeROE::computeRoeAVG(const StateVector& Wl, const StateVector& Wr) {

    FloatType rhoL = Wl[0];
    FloatType rhoR = Wr[0];
    FloatType u1L = Wl[1];
    FloatType u1R = Wr[1];
    FloatType u2L = Wl[2];
    FloatType u2R = Wr[2];
    FloatType u3L = Wl[3];
    FloatType u3R = Wr[3];
    FloatType htL = _fluid.computeTotalEnthalpy_rho_u_et(Wl[0], {Wl[1], Wl[2], Wl[3]}, Wl[4]);
    FloatType htR = _fluid.computeTotalEnthalpy_rho_u_et(Wr[0], {Wr[1], Wr[2], Wr[3]}, Wr[4]);

    _rhoAVG = std::sqrt(rhoL * rhoR);
    _u1AVG = roeAVG(rhoL, rhoR, u1L, u1R);
    _u2AVG = roeAVG(rhoL, rhoR, u2L, u2R);
    _u3AVG = roeAVG(rhoL, rhoR, u3L, u3R);
    _htAVG = roeAVG(rhoL, rhoR, htL, htR);
    _aAVG = std::sqrt((_fluid.getGamma() -1.0) * (_htAVG - 0.5 * (_u1AVG*_u1AVG + _u2AVG*_u2AVG + _u3AVG*_u3AVG)));
}


FloatType CAdvectionSchemeROE::roeAVG(FloatType& rhoL, FloatType& rhoR, FloatType& phiL, FloatType& phiR) const {
    FloatType avg = (std::sqrt(rhoL) * phiL + std::sqrt(rhoR) * phiR) / (std::sqrt(rhoL) + std::sqrt(rhoR));
    return avg;
}


void CAdvectionSchemeROE::computeEigenvalues() {
    _eigenvalues[0] = _u1AVG - _aAVG;
    _eigenvalues[1] = _u1AVG;
    _eigenvalues[2] = _u1AVG;
    _eigenvalues[3] = _u1AVG;
    _eigenvalues[4] = _u1AVG + _aAVG;
}

void CAdvectionSchemeROE::computeEigenvectors() {
    _eigenvectors[0] = StateVector({1.0, _u1AVG - _aAVG, _u2AVG, _u3AVG, _htAVG - _aAVG * _u1AVG});
    _eigenvectors[1] = StateVector({1.0, _u1AVG, _u2AVG, _u3AVG, 0.5 * (_u1AVG*_u1AVG + _u2AVG*_u2AVG + _u3AVG*_u3AVG)});
    _eigenvectors[2] = StateVector({0.0, 0.0, 1.0, 0.0, _u2AVG});
    _eigenvectors[3] = StateVector({0.0, 0.0, 0.0, 1.0, _u3AVG});
    _eigenvectors[4] = StateVector({1.0, _u1AVG + _aAVG, _u2AVG, _u3AVG, _htAVG + _aAVG * _u1AVG});
}

void CAdvectionSchemeROE::computeWaveStrengths(const StateVector& WnormL, const StateVector& WnormR) {
    FloatType deltaRho = WnormR[0] - WnormL[0];
    FloatType deltaU1 = WnormR[1] - WnormL[1];
    FloatType deltaU2 = WnormR[2] - WnormL[2];
    FloatType deltaU3 = WnormR[3] - WnormL[3];
    FloatType pL = _fluid.computePressure_rho_u_et(WnormL[0], {WnormL[1], WnormL[2], WnormL[3]}, WnormL[4]);
    FloatType pR = _fluid.computePressure_rho_u_et(WnormR[0], {WnormR[1], WnormR[2], WnormR[3]}, WnormR[4]);
    FloatType deltaP = pR - pL;
    

    _waveStrengths[0] = 1.0 / (2.0 * _aAVG * _aAVG) * (deltaP - _rhoAVG * _aAVG * deltaU1);
    _waveStrengths[1] = deltaRho - deltaP / (_aAVG * _aAVG);
    _waveStrengths[2] = _rhoAVG * deltaU2;
    _waveStrengths[3] = _rhoAVG * deltaU3;
    _waveStrengths[4] = 1.0 / (2.0 * _aAVG * _aAVG) * (deltaP + _rhoAVG * _aAVG * deltaU1);
}


void CAdvectionSchemeROE::computeRoeFlux(const Vector3D& S, const StateVector& WnormL, const StateVector& WnormR, StateVector& flux) const {

    // fluxes in normal reference frame
    StateVector fluxL = computeEulerFluxFromPrimitive(WnormL, {1.0, 0.0, 0.0}, _fluid);
    StateVector fluxR = computeEulerFluxFromPrimitive(WnormR, {1.0, 0.0, 0.0}, _fluid);

    StateVector fluxRoe = (fluxL + fluxR) * 0.5;
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            fluxRoe[i] -= 0.5 * _waveStrengths[j] * std::abs(_eigenvalues[j]) *_eigenvectors[j][i];
        }
    }

    // compute the flux in the original cartesian reference frame
    flux[0] = fluxRoe[0];
    flux[1] = fluxRoe[1] * (_n1.dot(_x1)) + fluxRoe[2] * (_n2.dot(_x1)) + fluxRoe[3] * (_n3.dot(_x1));
    flux[2] = fluxRoe[1] * (_n1.dot(_x2)) + fluxRoe[2] * (_n2.dot(_x2)) + fluxRoe[3] * (_n3.dot(_x2));
    flux[3] = fluxRoe[1] * (_n1.dot(_x3)) + fluxRoe[2] * (_n2.dot(_x3)) + fluxRoe[3] * (_n3.dot(_x3));
    flux[4] = fluxRoe[4];
}