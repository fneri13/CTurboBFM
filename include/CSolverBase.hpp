#pragma once

#include "types.hpp"
#include "Config.hpp"
#include "CMesh.hpp"
#include "CFluid.hpp"
#include "CAdvectionScheme.hpp"
#include "CJSTScheme.hpp"
#include "CBoundaryConditionBase.hpp"
#include "CBoundaryConditionEulerWall.hpp"
#include "CBoundaryConditionInlet.hpp"
#include "CBoundaryConditionOutlet.hpp"
#include "CBoundaryConditionInletSupersonic.hpp"
#include "CBoundaryConditionOutletSupersonic.hpp"

/** 
  *  \brief     Class handling base solver capabilities.
  *  \details   It computes the finite volume mesh.
  *  \author    Francesco Neri
  */
class CSolverBase {
    
    public:
        CSolverBase(Config &config, CMesh &mesh);
        
        virtual ~CSolverBase() {}

        void readBoundaryConditions();

        virtual void initializeSolutionArrays() = 0;

        virtual void solve() = 0;

        const std::array<int, 3> getStepMask(FluxDirection direction) const;

    protected:
        const Config& _config;
        
        const CMesh& _mesh;
        
        int _nDimensions {0}, _nPointsI, _nPointsJ, _nPointsK;
        
        Matrix3D<FloatType> _timeStep;
        
        std::unique_ptr<CFluid> _fluid;

        std::unique_ptr<CAdvectionScheme> _advectionScheme;
        
        std::map<BoundaryIndices, FloatType> _massFlows;
        
        std::map<BoundaryIndices, BoundaryType> _boundaryTypes;

        std::map<BoundaryIndices, std::vector<FloatType>> _boundaryValues;

        std::map<BoundaryIndices, std::unique_ptr<CBoundaryConditionBase>> _boundaryConditions;

        Topology _topology;

        

};