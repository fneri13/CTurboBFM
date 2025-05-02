#pragma once

#include "CSolverBase.hpp"
#include "fstream"

/** 
  *  \brief     Class handling base output capabilities.
  *  \details   interface for output of solutions.
  *  \author    Francesco Neri
  */

class COutputBase {
    
    public:
        COutputBase(const Config &config, const CMesh &mesh, const FlowSolution &solution, const CFluid &fluid) : _config(config), _mesh(mesh), _solution(solution),  _fluid(fluid){};
        
        virtual ~COutputBase() = default;

        virtual void writeSolution() = 0;

        void getScalarFieldsMap(std::map<std::string, Matrix3D<FloatType>>& scalarsMap);

    protected:
        const Config& _config;
        const CMesh& _mesh;
        const FlowSolution& _solution;
        const CFluid& _fluid;
};