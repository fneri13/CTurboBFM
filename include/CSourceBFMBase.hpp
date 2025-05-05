#pragma once

#include "types.hpp"
#include "CMesh.hpp"
#include "Config.hpp"
#include "CFluid.hpp"

// Base class handling BFM source terms computation
class CSourceBFMBase {
public:

    CSourceBFMBase(const Config &config, const CMesh &mesh, const CFluid &fluid, const FlowSolution &solution) 
        : _config(config), _mesh(mesh), _solution(solution), _fluid(fluid) {};

    virtual ~CSourceBFMBase() = default;  

    virtual StateVector computeSource(size_t i, size_t j, size_t k) const = 0;

protected:
    const CFluid& _fluid;
    const FlowSolution& _solution;
    const Config& _config;
    const CMesh& _mesh;
    
};
