#pragma once

#include "COutputBase.hpp"

/** 
  *  \brief     Class handling base output capabilities.
  *  \details   interface for output of solutions.
  *  \author    Francesco Neri
  */

class COutputCSV : public COutputBase {
    
    public:
        COutputCSV(const Config &config, const CMesh &mesh, const FlowSolution &solution, const CFluid &fluid, const Matrix3D<Vector3D> &inviscidForce, const Matrix3D<Vector3D> &viscousForce) 
                : COutputBase(config, mesh, solution, fluid, inviscidForce, viscousForce) {};
        
        virtual ~COutputCSV() = default;

        virtual void writeSolution(size_t iterationCounter) override;

};