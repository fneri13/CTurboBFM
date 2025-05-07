#pragma once

#include "CSourceBFMBase.hpp"

// Base class handling BFM Hall model source terms computation
class CSourceBFMHall : public CSourceBFMBase {
public:

    CSourceBFMHall(const Config &config, const CFluid &fluid, const CMesh &mesh) 
        : CSourceBFMBase(config, fluid, mesh) {}

    virtual ~CSourceBFMHall() = default;  

protected:

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive) override;
    
};
