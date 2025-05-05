#pragma once

#include "types.hpp"
#include "CFluid.hpp"

class CAdvectionSchemeBase {
public:

    CAdvectionSchemeBase(const CFluid &fluid) : _fluid(fluid) {};

    virtual ~CAdvectionSchemeBase() = default;  

    virtual StateVector computeFlux(const StateVector &Ull, 
        const StateVector &Ul, 
        const StateVector &Ur, 
        const StateVector &Urr,
        const Vector3D &surface) const = 0;

protected:
    const CFluid& _fluid;

    
};
