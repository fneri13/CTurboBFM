#pragma once

#include "CBoundaryConditionBase.hpp"

/** 
  *  \brief     Class handling radial equilibrium boundary condition capabilities.
  *  \details   It computes the boundary fluxes.
  *  \author    Francesco Neri
  */
class CBoundaryConditionRadialEquilibrium : public CBoundaryConditionBase {
    
    public:
        CBoundaryConditionRadialEquilibrium(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex, std::vector<FloatType>& radialProfile);
            
        virtual ~CBoundaryConditionRadialEquilibrium() {}

        virtual StateVector computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices) override;
        
    protected:
        std::vector<FloatType>& _radialPressureProfile;
};