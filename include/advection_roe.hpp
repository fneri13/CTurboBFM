#pragma once

#include "advection_base.hpp"


class AdvectionRoe : public AdvectionBase {
public:

    AdvectionRoe(const Config &config, FluidBase& fluid)  : AdvectionBase(config, fluid) {}

    StateVector computeFlux(
        const StateVector& Ull,
        const StateVector& Ul,
        const StateVector& Ur,
        const StateVector& Urr,
        const Vector3D& S) override;
    
    /** 
     * @brief Computes the orthonormal triad (normal, tangential, binormal) associated with the surface vector S.
     * The normal is aligned with S (left-to-right orientation).
     * */ 
    void computeNormalTriad(const Vector3D& S);

    /**
     * @brief Rotate the primitive variables into the local coordinate system defined by the normal triad.
     */
    StateVector computeRotatedPrimitive(const StateVector& W) const;

    void computeRoeAvgVariables(const StateVector& WnormL, const StateVector& WnormR);

    /**
     * @brief Compute Roe avg values of variable phi.
     */
    FloatType roeAverage(FloatType& rhoL, FloatType& rhoR, FloatType& phiL, FloatType& phiR) const ;

    void computeEigenvalues();

    void computeEigenvectors();

    void computeWaveStrengths(const StateVector& WnormL, const StateVector& WnormR);

    void assembleTotalFlux(
        const Vector3D& S, 
        const StateVector& WnormL, 
        const StateVector& WnormR, 
        StateVector& flux) const;

private:
        Vector3D _xVersor {1.0, 0.0, 0.0}; 
        Vector3D _yVersor {0.0, 1.0, 0.0}; 
        Vector3D _zVersor {0.0, 0.0, 1.0};
        Vector3D _n1, _n2, _n3; // normal, tangential, binormal versors
        FloatType _rhoAVG, _u1AVG, _u2AVG, _u3AVG, _htAVG, _aAVG;
        StateVector _eigenvalues;
        std::array<StateVector, 5> _eigenvectors;
        StateVector _waveStrengths;
        FloatType _entropyFixCoefficient {0.01}; // default value
        
};
