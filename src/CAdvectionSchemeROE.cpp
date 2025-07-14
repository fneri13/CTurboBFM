#include "CAdvectionSchemeROE.hpp"
#include "commonFunctions.hpp"

StateVector CAdvectionSchemeROE::computeFlux(
    const StateVector& Ull,
    const StateVector& Ul,
    const StateVector& Ur,
    const StateVector& Urr,
    const Vector3D& S) const
{
    auto Wll = getEulerPrimitiveFromConservative(Ull);
    auto Wl = getEulerPrimitiveFromConservative(Ul);
    auto Wr = getEulerPrimitiveFromConservative(Ur);
    auto Wrr = getEulerPrimitiveFromConservative(Urr);
    auto Wavg = (Wl + Wr) * 0.5; // element-wise addition

    
    StateVector flux({0.0, 0.0, 0.0, 0.0, 0.0});
    return flux;
}
