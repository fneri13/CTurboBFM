#pragma once
#include <string>
#include <vector>
#include "types.hpp"
#include "Config.hpp"
#include "CMesh.hpp"

/** 
  *  \brief     Class handling base solver capabilities.
  *  \details   It computes the finite volume mesh.
  *  \author    Francesco Neri
  */
class CSolverBase {
    
    public:
        CSolverBase(Config &config, CMesh &mesh);
        
        virtual ~CSolverBase() {}

        virtual void instantiateSolutionArrays() = 0;

        // virtual void initializeSolutionArrays() = 0;

        // virtual void solve() = 0;

        // virtual void computeFluxes() = 0;

        // virtual void spatialIntegration() = 0;

        // virtual void computeTimestepArray() = 0;

        // virtual void checkConvergence() = 0;

        // virtual void printInfoResiduals() = 0;

        // virtual void saveSolution() = 0;

        // virtual void computeGradient() = 0;

        // virtual void interpolateScalar() = 0;

    protected:
        const Config& _config;
        const CMesh& _mesh;
        int _nDimensions {0}, _nPointsI, _nPointsJ, _nPointsK;

};