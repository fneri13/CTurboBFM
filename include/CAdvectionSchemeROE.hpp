#pragma once

#include "CAdvectionSchemeBase.hpp"

/**
 * @brief Class for the ROE advection scheme.
 */
class CAdvectionSchemeROE : public CAdvectionSchemeBase {
    public:

        /**
         * @brief Constructs the Roe scheme with a given fluid reference. Original Roe formulation, taken from Toro.
         * @param fluid The fluid object.
         */
        CAdvectionSchemeROE(const CFluid& fluid) : CAdvectionSchemeBase(fluid) {}


        /**
         * @brief Compute the advective flux across a cell face.
         * This function uses the provided state vectors and surface normal to 
         * compute the flux contribution at a cell face using a specific numerical scheme.
         * @param Ull Second cell to the left of the interface.
         * @param Ul  First cell to the left of the interface.
         * @param Ur  First cell to the right of the interface.
         * @param Urr Second cell to the right of the interface.
         * @param surface The surface normal vector at the cell interface (oriented from left to right).
         * @return The computed flux as a StateVector.
         */
        StateVector computeFlux(
            const StateVector& Ull,
            const StateVector& Ul,
            const StateVector& Ur,
            const StateVector& Urr,
            const Vector3D& S) override;
        
        /**
         * @brief compute the three versors in the normal direction, tangential, and binormal
         * @param S The surface normal vector at the cell interface (oriented from left to right).
         */
        void computeNormalTriad(const Vector3D& S) ;

        /**
         * @brief compute the rotated primitve vector in the normal triad frame
         * @param S The surface normal vector at the cell interface (oriented from left to right).
         * @return The rotated state vector
         */
        StateVector computeRotatedPrimitive(const StateVector& W) const;

        void computeRoeAVG(const StateVector& WnormL, const StateVector& WnormR);


        FloatType roeAVG(FloatType& rhoL, FloatType& rhoR, FloatType& phiL, FloatType& phiR) const ;

        void computeEigenvalues();

        void computeEigenvectors();

        void computeWaveStrengths(const StateVector& WnormL, const StateVector& WnormR);

        void computeRoeFlux(const Vector3D& S, const StateVector& WnormL, const StateVector& WnormR, StateVector& flux) const;

    private:
            
            Vector3D _x1, _x2, _x3; 
            Vector3D _n1, _n2, _n3; // normal, tangential, binormal versors
            FloatType _rhoAVG, _u1AVG, _u2AVG, _u3AVG, _htAVG, _aAVG;
            StateVector _eigenvalues;
            std::array<StateVector, 5> _eigenvectors;
            StateVector _waveStrengths;
        
};
