#pragma once

#include "types.hpp"
#include "CFluidBase.hpp"
#include "CFluidIdeal.hpp"
#include "Config.hpp"

/**
 * @brief Abstract base class for advection flux calculation schemes.
 * Provides a common interface for numerical schemes used to compute
 * intercell fluxes in fluid dynamics simulations.
 */
class CAdvectionSchemeBase {
    public:

        /**
         * @brief Constructs the advection scheme with a given fluid reference.
         * @param fluid The fluid object.
         */
        CAdvectionSchemeBase(const Config &config, const CFluidBase &fluid);

        virtual ~CAdvectionSchemeBase() = default;  

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
        virtual StateVector computeFlux(const StateVector &Ull, 
            const StateVector &Ul, 
            const StateVector &Ur, 
            const StateVector &Urr,
            const Vector3D &surface) = 0;
        
        void reconstructMUSCL(StateVector& Wll, StateVector& Wl, StateVector& Wr, StateVector& Wrr, FluxLimiter fluxLimiter) const;

        StateVector computeSmoothness(const StateVector& Wl, const StateVector& Wc, const StateVector& Wr) const;

        StateVector computeLimiter(const StateVector& smoothness, FluxLimiter fluxLimiter) const;

    protected:
        const CFluidBase& _fluid;  // fluid object used to compute the thermodynamic properties
        const Config& _config; // configuration object for accessing simulation parameters
        bool _MUSCL = false; // flag for MUSCL scheme, default is false
        FluxLimiter _fluxLimiter = FluxLimiter::NONE; // flag for flux limiter

    
};
