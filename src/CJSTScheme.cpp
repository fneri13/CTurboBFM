#include "../include/CJSTScheme.hpp"
#include "../include/commonFunctions.hpp"

CJSTScheme::StateVector CJSTScheme::computeFlux(
    const StateVector& Ull,
    const StateVector& Ul,
    const StateVector& Ur,
    const StateVector& Urr,
    const SurfaceVector& S) const
{

    auto area = S.magnitude();
    auto S_dir = S / area;
    auto U_avg = 0.5 * (Ul + Ur); // element-wise addition

    auto Wll = getEulerPrimitiveFromConservative(Ull);
    auto Wl = getEulerPrimitiveFromConservative(Ul);
    auto Wr = getEulerPrimitiveFromConservative(Ur);
    auto Wrr = getEulerPrimitiveFromConservative(Urr);

    // float r_factors[2] = {
    //     _fluid.computeSoundSpeed_rho_u_et(Wl[0], {Wl[1], Wl[2], Wl[3]}, Wl[4]),
    //     _fluid.computeSoundSpeed_rho_u_et(Wr[0], {Wr[1], Wr[2], Wr[3]}, Wr[4])
    // };

    // float s_factors[2] = {
    //     std::abs((Wll[4] - 2 * Wl[4] + Wr[4]) / (Wll[4] + 2 * Wl[4] + Wr[4])),
    //     std::abs((Wl[4] - 2 * Wr[4] + Wrr[4]) / (Wl[4] + 2 * Wr[4] + Wrr[4]))
    // };

    // float r = std::max(r_factors[0], r_factors[1]);
    // float s = std::max(s_factors[0], s_factors[1]);

    // float psi2 = _kappa2 * s * r;
    // float psi4 = std::max(0.0f, _kappa4 * r - _c4 * psi2);

    // auto flux = EulerFluxFromConservatives(U_avg, S, _fluid);
    // auto dissipation = psi2 * (Ur - Ul) - psi4 * ((Urr - Ur) - 2 * (Ur - Ul) + (Ul - Ull));

    // for (size_t i = 0; i < 5; ++i) {
    //     flux[i] -= dissipation[i];
    // }

    return flux;
}
