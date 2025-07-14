#pragma once

#include "CAdvectionSchemeBase.hpp"

/**
 * @brief Class for the JST advection scheme.
 */
class CAdvectionSchemeJST : public CAdvectionSchemeBase {
    public:

        /**
         * @brief Constructs the JST scheme with a given fluid reference. Original Jameson formulation.
         * @param fluid The fluid object.
         */
        CAdvectionSchemeJST(const CFluid& fluid) : CAdvectionSchemeBase(fluid) {}


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
         * @brief Compute the r-factor.
         * @param primitive vector of primitive variables.
         * @return The r-factor needed by the JST scheme.
         */    
        FloatType computeRFactor(const StateVector primitive) const ;
        

        /**
         * @brief Compute the s-factor.
         * @param prim1 first vector of primitive variables.
         * @param prim2 second vector of primitive variables.
         * @param prim3 third vector of primitive variables.
         * @return The s-factor needed by the JST scheme.
         */   
        FloatType computeSFactor(const StateVector prim1, const StateVector prim2, const StateVector prim3) const;

    private:
        FloatType _kappa2 {0.5}, _kappa4 {0.02}, _c4 {2.0}; // coefficients of the JST scheme (same values of SU2)
        
            
};
