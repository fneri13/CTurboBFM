#pragma once

#include "CBoundaryConditionBase.hpp"

/** 
  *  \brief     Class handling outlet boundary condition capabilities.
  *  \details   It computes the boundary fluxes.
  *  \author    Francesco Neri
  */
class CBoundaryConditionOutlet : public CBoundaryConditionBase {
    
    public:
        CBoundaryConditionOutlet(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex, std::vector<FloatType> bcValues);
            
        virtual ~CBoundaryConditionOutlet() {}

        virtual StateVector computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices) override;
        
    protected:
        std::vector<FloatType> _boundaryValues;
};