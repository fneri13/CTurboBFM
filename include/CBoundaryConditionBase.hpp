#pragma once

#include "types.hpp"
#include "Config.hpp"
#include "CMesh.hpp"
#include "CFluid.hpp"

/** 
  *  \brief     Class handling base boundary condition capabilities.
  *  \details   It computes the boundary fluxes.
  *  \author    Francesco Neri
  */
class CBoundaryConditionBase {
    
    public:
        CBoundaryConditionBase(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex) 
            : _config(config), _mesh(mesh), _fluid(fluid), _boundaryIndex(boundIndex){};
            
        virtual ~CBoundaryConditionBase() {}

        virtual StateVector computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices) = 0;
        
    protected:
        const Config& _config;
        const CMesh& _mesh;
        const CFluid& _fluid;
        const BoundaryIndices _boundaryIndex;
        
};