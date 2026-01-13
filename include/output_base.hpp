#pragma once

#include "solver_base.hpp"
#include "fstream"
#include <filesystem>

/** 
  *  \brief     Class handling base output capabilities.
  *  \details   interface for output of solutions.
  *  \author    Francesco Neri
  */

class OutputBase {
    
    public:
        OutputBase(const Config &config, const Mesh &mesh, const FlowSolution &solution, const FluidBase &fluid, 
                    const Matrix3D<Vector3D> &inviscidForce, const Matrix3D<Vector3D> &viscousForce, 
                    const Matrix3D<FloatType> &deviationAngle);
        
        virtual ~OutputBase() = default;

        virtual void writeSolution(size_t iterationCounter, bool alsoGradients=false) = 0;

        void getScalarFieldsMap(std::map<std::string, Matrix3D<FloatType>>& scalarsMap, bool alsoGradients=false) const;

    protected:
        const Config& _config;
        
        const Mesh& _mesh;
        
        const FlowSolution& _solution;
        
        const FluidBase& _fluid;
        
        const Matrix3D<Vector3D>& _inviscidForce;
        
        const Matrix3D<Vector3D>& _viscousForce;
        
        const Matrix3D<FloatType>& _deviationAngle;
        
        bool _isUnsteadyOutput;

        std::string getOutputFilename(size_t iterationCounter);
        
        std::string _outputDirectory = "Volume_CSV";

        OutputFields _outputFields;
};