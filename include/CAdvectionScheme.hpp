#pragma once
#include <vector>
#include "types.hpp"
#include "CFluid.hpp"

class CAdvectionScheme {
public:

    CAdvectionScheme(const CFluid &fluid) : _fluid(fluid) {};

    virtual ~CAdvectionScheme() = default;  

    virtual StateVector computeFlux(const StateVector &Ull, 
        const StateVector &Ul, 
        const StateVector &Ur, 
        const StateVector &Urr,
        const Vector3D &surface) const = 0;

protected:
    const CFluid& _fluid;

    
};
