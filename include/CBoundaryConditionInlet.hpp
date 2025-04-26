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

        // Formulation taken from 'Formulation and Implementation of Inflow/Outflow Boundary Conditions to Simulate Propulsive Effects', Rodriguez et al.
        virtual StateVector computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint) override;
        
    protected:
        std::vector<FloatType> _boundaryValues;
};