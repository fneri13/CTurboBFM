#pragma once

#include "types.hpp"
#include "CMesh.hpp"
#include "Config.hpp"
#include "CFluid.hpp"
#include "commonFunctions.hpp"

// Base class handling BFM source terms computation
class CSourceBFMBase {
public:

    CSourceBFMBase(const Config &config, const CMesh &mesh, const CFluid &fluid, const FlowSolution &solution) 
        : _config(config), _mesh(mesh), _fluid(fluid), _solution(solution) {};

    virtual ~CSourceBFMBase() = default;  

    virtual StateVector computeSource(size_t i, size_t j, size_t k);

    virtual StateVector computeBlockageSource(size_t i, size_t j, size_t k);

    virtual StateVector computeBodyForceSource(size_t i, size_t j, size_t k);

protected:
    const Config& _config;
    const CMesh& _mesh;
    const CFluid& _fluid;
    const FlowSolution& _solution;
    
};
