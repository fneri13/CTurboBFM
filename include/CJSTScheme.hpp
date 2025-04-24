#pragma once

#include "CAdvectionScheme.hpp"

class CJSTScheme : public CAdvectionScheme {
    public:
        using CAdvectionScheme::StateVector;
        using CAdvectionScheme::SurfaceVector;

        CJSTScheme(const CFluid& fluid) : CAdvectionScheme(fluid) {}

        StateVector computeFlux(
            const StateVector& Ull,
            const StateVector& Ul,
            const StateVector& Ur,
            const StateVector& Urr,
            const SurfaceVector& S) const override;

    private:
        FloatType _kappa2 {0.2}, _kappa4 {0.5}, _c4 {2.0};
            
};
