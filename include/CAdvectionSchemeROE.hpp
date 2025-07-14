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
            const Vector3D& S) const override;

    private:
        
        
            
};
