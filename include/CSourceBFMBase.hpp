#pragma once

#include "types.hpp"
#include "CMesh.hpp"
#include "Config.hpp"
#include "CFluid.hpp"
#include "commonFunctions.hpp"

// Base class handling BFM source terms computation
class CSourceBFMBase {
public:

    CSourceBFMBase(const Config &config, const CFluid &fluid, const CMesh &mesh) 
        : _config(config), _fluid(fluid), _mesh(mesh) {};

    virtual ~CSourceBFMBase() = default;  

    virtual StateVector computeSource(size_t i, size_t j, size_t k, const StateVector& primitive) const;

    virtual StateVector computeBlockageSource(size_t i, size_t j, size_t k, const StateVector& primitive) const;

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k, const StateVector& primitive) const;

protected:
    const Config& _config;
    const CFluid& _fluid;
    const CMesh& _mesh;
    
    
};
