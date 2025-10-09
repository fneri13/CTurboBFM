#pragma once

#include "CSolverBase.hpp"
#include "fstream"
#include <filesystem>

/** 
  *  \brief     Class handling base output capabilities.
  *  \details   interface for output of solutions.
  *  \author    Francesco Neri
  */

class COutputBase {
    
    public:
        COutputBase(const Config &config, const CMesh &mesh, const FlowSolution &solution, const CFluidBase &fluid, 
                    const Matrix3D<Vector3D> &inviscidForce, const Matrix3D<Vector3D> &viscousForce, 
                    const Matrix3D<FloatType> &deviationAngle);
        
        virtual ~COutputBase() = default;

        virtual void writeSolution(size_t iterationCounter) = 0;

        void getScalarFieldsMap(std::map<std::string, Matrix3D<FloatType>>& scalarsMap) const;

    protected:
        const Config& _config;
        const CMesh& _mesh;
        const FlowSolution& _solution;
        const CFluidBase& _fluid;
        const Matrix3D<Vector3D>& _inviscidForce;
        const Matrix3D<Vector3D>& _viscousForce;
        const Matrix3D<FloatType>& _deviationAngle;
        bool _isUnsteadyOutput;

        std::string getOutputFilename(size_t iterationCounter);
        std::string _outputDirectory = "volumeData";
};