#pragma once

#include "CBoundaryConditionBase.hpp"

/** 
  *  \brief     Class handling Euler wall boundary condition capabilities.
  *  \details   It computes the boundary fluxes.
  *  \author    Francesco Neri
  */
class CBoundaryConditionEulerWall : public CBoundaryConditionBase {
    
    public:
        CBoundaryConditionEulerWall(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex);
            
        virtual ~CBoundaryConditionEulerWall() {}

        virtual StateVector computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices) override;
        
    protected:
        
};