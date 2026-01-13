#pragma once

#include "output_base.hpp"

/** 
  *  \brief     Class handling base output capabilities.
  *  \details   interface for output of solutions.
  *  \author    Francesco Neri
  */

class OutputCSV : public OutputBase {
    
    public:
        OutputCSV(const Config &config, const Mesh &mesh, const FlowSolution &solution, const FluidBase &fluid, const Matrix3D<Vector3D> &inviscidForce, const Matrix3D<Vector3D> &viscousForce, const Matrix3D<FloatType> &deviationAngle) 
                : OutputBase(config, mesh, solution, fluid, inviscidForce, viscousForce, deviationAngle) {};
        
        virtual ~OutputCSV() = default;

        virtual void writeSolution(size_t iterationCounter, bool alsoGradients=false) override;

};