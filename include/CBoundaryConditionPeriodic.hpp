#pragma once

#include "CBoundaryConditionBase.hpp"

/**
 * @brief Class for periodic boundary condition flux calculation.
 */
class CBoundaryConditionPeriodic : public CBoundaryConditionBase {
    
    public:

        /**
         * @brief Constructs the object with given references.
         * @param config The configuration object.
         * @param mesh The mesh object.
         * @param fluid The fluid object.
         * @param boundIndex The boundary index.
         * @param inletValues The inlet values (total pressure, total temperature, flow direction).
         */
        CBoundaryConditionPeriodic(const Config &config, const CMesh &mesh, CFluid &fluid, BoundaryIndices boundIndex, std::vector<FloatType> inletValues);
            
        virtual ~CBoundaryConditionPeriodic() {}

        /**
         * @brief Compute the boundary flux.
         * @param internalConservative The internal point conservative variables.
         * @param surface The surface normal vector (also not normalized).
         * @param midPoint The midpoint of the boundary face.
         * @param indices The indices (i,j,k) of the boundary face.
         * @return The boundary flux.
         */
        virtual StateVector computeBoundaryFlux(StateVector internalConservative, Vector3D surface, Vector3D midPoint, std::array<size_t, 3> indices, const FlowSolution &flowSolution, const size_t iterCounter) override;
        
    protected:
        std::vector<FloatType> _boundaryValues;
        FloatType _periodicityAngle = 0.0;
};