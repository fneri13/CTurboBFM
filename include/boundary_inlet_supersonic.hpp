#pragma once

#include "boundary_base.hpp"

/** 
  *  \brief     Class handling supersonic inlet boundary condition capabilities.
  *  \details   It computes the boundary fluxes.
  *  \author    Francesco Neri
  */
class BoundaryInletSupersonic : public BoundaryBase {
    
    public:
        BoundaryInletSupersonic(const Config &config, const Mesh &mesh, FluidBase &fluid, BoundaryIndices boundIndex, std::vector<FloatType> inletValues);
            
        virtual ~BoundaryInletSupersonic() {}

        // Formulation taken from 'Formulation and Implementation of Inflow/Outflow Boundary Conditions to Simulate Propulsive Effects', Rodriguez et al.
        virtual StateVector computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &flowSolution, const size_t iterCounter) override;
        
    protected:
        std::vector<FloatType> _boundaryValues;
};