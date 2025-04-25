#pragma once

#include "CBoundaryConditionBase.hpp"

/** 
  *  \brief     Class handling inlet boundary condition capabilities.
  *  \details   It computes the boundary fluxes.
  *  \author    Francesco Neri
  */
class CBoundaryConditionInlet : public CBoundaryConditionBase {
    
    public:
        CBoundaryConditionInlet(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex, std::vector<FloatType> inletValues);
            
        virtual ~CBoundaryConditionInlet() {}

        virtual StateVector computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint) override;
        
    protected:
        std::vector<FloatType> _boundaryValues;
};