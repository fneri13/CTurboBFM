#include "CAdvectionSchemeBase.hpp"

CAdvectionSchemeBase::CAdvectionSchemeBase(const Config &config, const CFluid &fluid)
    : _fluid(fluid), _config(config) 
{
    _MUSCL = _config.getMUSCLreconstruction();
    _fluxLimiter = _config.getFluxLimiter();
}

void CAdvectionSchemeBase::reconstructMUSCL(StateVector& Wll, StateVector& Wl, StateVector& Wr, StateVector& Wrr, FluxLimiter fluxLimiter) const {
    StateVector smoothnessL = computeSmoothness(Wll, Wl, Wr);
    StateVector smoothnessR = computeSmoothness(Wr, Wl, Wrr);

    StateVector limiterL = computeLimiter(smoothnessL, fluxLimiter);
    StateVector limiterR = computeLimiter(smoothnessR, fluxLimiter);

    // if (std::abs(limiterL[0]) > 1e-2 || std::abs(limiterR[0]) > 1e-2) {
    //     std::cout << "MUSCL reconstruction is applied with flux limiter: " << static_cast<int>(fluxLimiter) << std::endl;
    // }

    Wl += limiterL * (Wr - Wl) * 0.5;
    Wr -= limiterR * (Wrr - Wr) * 0.5;
}


StateVector CAdvectionSchemeBase::computeSmoothness(const StateVector& Wl, const StateVector& Wc, const StateVector& Wr) const {
    StateVector smoothness;
    
    for (std::size_t i = 0; i < 5; ++i) {
        smoothness[i] = (Wc[i] - Wl[i]) * (Wr[i] - Wc[i] + 1E-8);
    }

    return smoothness;
}

StateVector CAdvectionSchemeBase::computeLimiter(const StateVector& smoothness, FluxLimiter fluxLimiter) const{
    StateVector limiter({0.0, 0.0, 0.0, 0.0, 0.0});

    if (fluxLimiter == FluxLimiter::NONE) {
        return limiter; // No limiter applied
    }
    else if (fluxLimiter == FluxLimiter::VAN_ALBADA) {
        for (std::size_t i = 0; i < 5; ++i) {
            limiter[i] = (smoothness[i]*smoothness[i] + smoothness[i]) / (1.0 + smoothness[i]*smoothness[i]);
        }
    }
    else if (fluxLimiter == FluxLimiter::VAN_LEER) {
        for (std::size_t i = 0; i < 5; ++i) {
            limiter[i] = (smoothness[i] + std::abs(smoothness[i])) / (1.0 + std::abs(smoothness[i]));
        }
    }
    else if (fluxLimiter == FluxLimiter::MIN_MOD) {
        for (std::size_t i = 0; i < 5; ++i) {
            limiter[i] = std::max(0.0, std::min(1.0, smoothness[i]));
        }
    }

    return limiter;
}