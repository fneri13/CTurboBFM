#pragma once

#include "CBoundaryConditionBase.hpp"

/** 
  *  \brief     Class handling wedge boundary condition capabilities for axisymmetric simulations.
  *  \details   It computes the boundary fluxes.
  *  \author    Francesco Neri
  */
class CBoundaryConditionWedge : public CBoundaryConditionBase {
    
    public:
        CBoundaryConditionWedge(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex);
            
        virtual ~CBoundaryConditionWedge() {}

        virtual StateVector computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint) override;
        
    private:
        
};