#pragma once

#include "CAdvectionScheme.hpp"

// JST scheme advection flux formulation by Jameson
class CJSTScheme : public CAdvectionScheme {
    public:

        CJSTScheme(const CFluid& fluid) : CAdvectionScheme(fluid) {}

        StateVector computeFlux(
            const StateVector& Ull,
            const StateVector& Ul,
            const StateVector& Ur,
            const StateVector& Urr,
            const Vector3D& S) const override;
        
        FloatType computeRFactor(const StateVector primitive) const ;
        
        FloatType computeSFactor(const StateVector prim1, const StateVector prim2, const StateVector prim3) const;

    private:
        FloatType _kappa2 {0.2}, _kappa4 {0.5}, _c4 {2.0};
        
            
};
