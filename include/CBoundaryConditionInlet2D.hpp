#pragma once

#include "CBoundaryConditionInlet.hpp"
#include "CBoundaryConditionBase.hpp"
#include "CInput.hpp"

/**
 * @brief Class for inlet defined by 2D grids for boundary condition flux calculation.
 * Provides a common interface for it.
 */
class CBoundaryConditionInlet2D : public CBoundaryConditionBase {
    
    public:

        /**
         * @brief Constructs the object with given references.
         * @param config The configuration object.
         * @param mesh The mesh object.
         * @param fluid The fluid object.
         * @param boundIndex The boundary index.
         * @param inletValues The inlet values (total pressure, total temperature, flow direction).
         */
        CBoundaryConditionInlet2D(const Config &config, const CMesh &mesh, CFluidBase &fluid, BoundaryIndices boundIndex, std::string inletFilePath);
            
        virtual ~CBoundaryConditionInlet2D() {}

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
        std::string _inletFilePath;
        
        ReferenceFrame _referenceFrame;

        CInput _inputGrid;

        size_t _inputNi, _inputNj, _inputNk;
};