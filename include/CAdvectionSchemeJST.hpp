#pragma once

#include "CAdvectionSchemeBase.hpp"

// JST scheme advection flux formulation by Jameson
class CAdvectionSchemeJST : public CAdvectionSchemeBase {
    public:

        CAdvectionSchemeJST(const CFluid& fluid) : CAdvectionSchemeBase(fluid) {}

        StateVector computeFlux(
            const StateVector& Ull,
            const StateVector& Ul,
            const StateVector& Ur,
            const StateVector& Urr,
            const Vector3D& S) const override;
        
        FloatType computeRFactor(const StateVector primitive) const ;
        
        FloatType computeSFactor(const StateVector prim1, const StateVector prim2, const StateVector prim3) const;

    private:
        FloatType _kappa2 {0.5}, _kappa4 {0.02}, _c4 {2.0};
        
            
};
