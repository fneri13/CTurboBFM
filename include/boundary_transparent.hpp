#pragma once

#include "boundary_base.hpp"
#include "advection_base.hpp"
#include "advection_roe.hpp"
#include "advection_jst.hpp"

/** 
  *  \brief     Class handling wedge (periodic boundaries) boundary condition capabilities for axisymmetric simulations.
  *  \details   It computes the boundary fluxes.
  *  \author    Francesco Neri
  */
class BoundaryTransparent : public BoundaryBase {
    
    public:

        /**
         * @brief Constructs the object with given references.
         * @param config The configuration object.
         * @param mesh The mesh object.
         * @param fluid The fluid object.
         * @param boundIndex The boundary index.
         */
        BoundaryTransparent(const Config &config, const Mesh &mesh, FluidBase &fluid, BoundaryIndices boundIndex, AdvectionBase &advScheme);
            
        virtual ~BoundaryTransparent() {}

        /**
         * @brief Compute the boundary flux.
         * Approach using ghost points, where the ghost cell values are the same as the internal cell values
         * @param internalConservative The internal point conservative variables.
         * @param surface The surface normal vector (also not normalized).
         * @param midPoint The midpoint of the boundary face.
         * @param indices The indices (i,j,k) of the boundary face.
         * @return The boundary flux.
         */
        virtual StateVector computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &flowSolution, const size_t iterCounter) override;

    private:
        AdvectionBase &_advScheme;
        
};