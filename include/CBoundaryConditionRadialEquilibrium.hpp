#pragma once

#include "CBoundaryConditionBase.hpp"

/** 
  *  \brief     Class handling outlet boundary with radial equilibrium.
  *  \details   It computes the boundary fluxes.
  *  \author    Francesco Neri
  */
class CBoundaryConditionRadialEquilibrium : public CBoundaryConditionBase {
    
    public:
        
        /** 
         * @brief Constructs the object with given references.
         * @param config The configuration object.
         * @param mesh The mesh object.
         * @param fluid The fluid object.
         * @param boundIndex The boundary index.
         * @param radialProfile reference to the radial profile of pressure used in the boundary condition flux evaluation.
         */
        CBoundaryConditionRadialEquilibrium(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex, std::vector<FloatType>& radialProfile);
            
        virtual ~CBoundaryConditionRadialEquilibrium() {}
        
        /**
         * @brief Compute the boundary flux.
         * Formulation taken from 'Formulation and Implementation of Inflow/Outflow Boundary Conditions to Simulate Propulsive Effects', Rodriguez et al.
         * @param internalConservative The internal point conservative variables.
         * @param surface The surface normal vector (also not normalized).
         * @param midPoint The midpoint of the boundary face.
         * @param indices The indices (i,j,k) of the boundary face.
         * @return The boundary flux.
         */
        virtual StateVector computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &flowSolution, const size_t iterCounter) override;
        
    protected:
        std::vector<FloatType>& _radialPressureProfile; // Reference to the radial profile of the static pressure (which is updated at every iteration)
};