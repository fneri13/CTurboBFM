#pragma once

#include "types.hpp"
#include "config.hpp"
#include "mesh.hpp"
#include "fluid_base.hpp"
#include "fluid_ideal.hpp"

/**
 * @brief Abstract base class for boundary condition flux calculation.
 * Provides a common interface for it.
 */
class BoundaryBase {
    
    public:
        
        /**
         * @brief Constructs the object with given references.
         * @param config The configuration object.
         * @param mesh The mesh object.
         * @param fluid The fluid object.
         * @param boundIndex The boundary index.
         */
        BoundaryBase(const Config &config, const Mesh &mesh, FluidBase &fluid, BoundaryIndices boundIndex) 
            : _config(config), _mesh(mesh), _fluid(fluid), _boundaryIndex(boundIndex){};
            
        virtual ~BoundaryBase() {}


        /**
         * @brief Compute the boundary flux.
         * @param internalConservative The internal point conservative variables.
         * @param surface The surface normal vector (also not normalized).
         * @param midPoint The midpoint of the boundary face.
         * @param indices The indices (i,j,k) of the boundary face.
         * @param flowSolution The entire flow solution.
         */
        virtual StateVector computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &flowSolution, const size_t iterCounter) = 0;


        StateVector computeSubsonicsInletFlux(StateVector internalConservative, Vector3D surface, 
                                        FloatType totPressure, FloatType totTemperature, Vector3D flowDirection);
        
    protected:
        const Config& _config;
        const Mesh& _mesh;
        const FluidBase& _fluid;
        const BoundaryIndices _boundaryIndex;
        
};