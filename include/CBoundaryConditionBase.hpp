#pragma once

#include "types.hpp"
#include "Config.hpp"
#include "CMesh.hpp"
#include "CFluid.hpp"

/**
 * @brief Abstract base class for boundary condition flux calculation.
 * Provides a common interface for it.
 */
class CBoundaryConditionBase {
    
    public:
        
        /**
         * @brief Constructs the object with given references.
         * @param config The configuration object.
         * @param mesh The mesh object.
         * @param fluid The fluid object.
         * @param boundIndex The boundary index.
         */
        CBoundaryConditionBase(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex) 
            : _config(config), _mesh(mesh), _fluid(fluid), _boundaryIndex(boundIndex){};
            
        virtual ~CBoundaryConditionBase() {}


        /**
         * @brief Compute the boundary flux.
         * @param internalConservative The internal point conservative variables.
         * @param surface The surface normal vector (also not normalized).
         * @param midPoint The midpoint of the boundary face.
         * @param indices The indices (i,j,k) of the boundary face.
         */
        virtual StateVector computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &flowSolution) = 0;
        
    protected:
        const Config& _config;
        const CMesh& _mesh;
        const CFluid& _fluid;
        const BoundaryIndices _boundaryIndex;
        
};