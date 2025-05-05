#pragma once

#include "CBoundaryConditionBase.hpp"

/**
 * @brief Class for Euler walls (inviscid) boundary condition flux calculation.
 */
class CBoundaryConditionEulerWall : public CBoundaryConditionBase {
    
    public:
        
        /**
         * @brief Constructs the object with given references.
         * @param config The configuration object.
         * @param mesh The mesh object.
         * @param fluid The fluid object.
         * @param boundIndex The boundary index.
         */
        CBoundaryConditionEulerWall(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex);
            
        virtual ~CBoundaryConditionEulerWall() {}

        /**
         * @brief Compute the boundary flux.
         * @param internalConservative The internal point conservative variables.
         * @param surface The surface normal vector (also not normalized).
         * @param midPoint The midpoint of the boundary face.
         * @param indices The indices (i,j,k) of the boundary face.
         * @return The boundary flux.
         */
        virtual StateVector computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices) override;
        
        
};