#pragma once
#include <vector>
#include "types.hpp"
#include "CFluid.hpp"

class CAdvectionScheme {
public:
    using StateVector = std::array<double, 5>;
    using SurfaceVector = Vector3D;

    CAdvectionScheme(const CFluid &fluid) : _fluid(fluid) {};

    virtual StateVector computeFlux(const StateVector &Ull, 
        const StateVector &Ul, 
        const StateVector &Ur, 
        const StateVector &Urr,
        const SurfaceVector &surface) const = 0;

protected:
    const CFluid& _fluid;
    
};
