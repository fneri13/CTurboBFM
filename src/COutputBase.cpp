#include "COutputBase.hpp"
#include "commonFunctions.hpp"


COutputBase::COutputBase(const Config &config, const CMesh &mesh, const FlowSolution &solution, const CFluid &fluid, 
                        const Matrix3D<Vector3D> &inviscidForce, const Matrix3D<Vector3D> &viscousForce,
                        const Matrix3D<FloatType> &deviationAngle)

    : _config(config), _mesh(mesh), _solution(solution), _fluid(fluid), _inviscidForce(inviscidForce), _viscousForce(viscousForce), _deviationAngle(deviationAngle) {
        _isUnsteadyOutput = _config.saveUnsteadySolution();
        std::filesystem::create_directory(_outputDirectory);
    }



void COutputBase::getScalarFieldsMap(std::map<std::string, Matrix3D<FloatType>>& scalarsMap) {
    size_t ni = _mesh.getNumberPointsI();
    size_t nj = _mesh.getNumberPointsJ();
    size_t nk = _mesh.getNumberPointsK();
    
    // primary solution variables
    scalarsMap["Density"] = _solution.getDensity();
    scalarsMap["Velocity X"] = _solution.getVelocityX();
    scalarsMap["Velocity Y"] = _solution.getVelocityY();
    scalarsMap["Velocity Z"] = _solution.getVelocityZ();
    scalarsMap["Total Energy"] = _solution.getTotalEnergy();
    
    // auxillary solution variables
    scalarsMap["Pressure"] = Matrix3D<FloatType>(ni, nj, nk);
    scalarsMap["Temperature"] = Matrix3D<FloatType>(ni, nj, nk);
    scalarsMap["Mach"] = Matrix3D<FloatType>(ni, nj, nk);
    scalarsMap["Total Pressure"] = Matrix3D<FloatType>(ni, nj, nk);
    scalarsMap["Total Temperature"] = Matrix3D<FloatType>(ni, nj, nk);
    scalarsMap["Entropy"] = Matrix3D<FloatType>(ni, nj, nk);

    // BFM simulation additional variables
    bool isBFMActive = _config.isBFMActive();
    if (isBFMActive){
        scalarsMap["Blockage"] = _mesh.getInputFields(FieldNames::BLOCKAGE);
        scalarsMap["Relative Mach"] = Matrix3D<FloatType>(ni, nj, nk);
        scalarsMap["Grid Velocity X"] = Matrix3D<FloatType>(ni, nj, nk);
        scalarsMap["Grid Velocity Y"] = Matrix3D<FloatType>(ni, nj, nk);
        scalarsMap["Grid Velocity Z"] = Matrix3D<FloatType>(ni, nj, nk);
        scalarsMap["Relative Velocity X"] = Matrix3D<FloatType>(ni, nj, nk);
        scalarsMap["Relative Velocity Y"] = Matrix3D<FloatType>(ni, nj, nk);
        scalarsMap["Relative Velocity Z"] = Matrix3D<FloatType>(ni, nj, nk);
        scalarsMap["Viscous Body Force X"] = Matrix3D<FloatType>(ni, nj, nk);
        scalarsMap["Viscous Body Force Y"] = Matrix3D<FloatType>(ni, nj, nk);
        scalarsMap["Viscous Body Force Z"] = Matrix3D<FloatType>(ni, nj, nk);
        scalarsMap["Inviscid Body Force X"] = Matrix3D<FloatType>(ni, nj, nk);
        scalarsMap["Inviscid Body Force Y"] = Matrix3D<FloatType>(ni, nj, nk);
        scalarsMap["Inviscid Body Force Z"] = Matrix3D<FloatType>(ni, nj, nk);
        scalarsMap["Deviation Angle"] = Matrix3D<FloatType>(ni, nj, nk);
    }

    FloatType rho, ux, uy, uz, et;
    for (size_t i = 0; i < ni; ++i) {
        for (size_t j = 0; j < nj; ++j) {
            for (size_t k = 0; k < nk; ++k) {
                rho = scalarsMap["Density"](i, j, k);
                ux = scalarsMap["Velocity X"](i, j, k);
                uy = scalarsMap["Velocity Y"](i, j, k);
                uz = scalarsMap["Velocity Z"](i, j, k);
                et = scalarsMap["Total Energy"](i, j, k);

                scalarsMap["Pressure"](i, j, k) = _fluid.computePressure_rho_u_et(rho, {ux, uy, uz}, et);
                scalarsMap["Temperature"](i, j, k) = _fluid.computeTemperature_rho_u_et(rho, {ux, uy, uz}, et);
                scalarsMap["Mach"](i, j, k) = _fluid.computeMachNumber_rho_u_et(rho, {ux, uy, uz}, et);
                scalarsMap["Total Pressure"](i, j, k) = _fluid.computeTotalPressure_rho_u_et(rho, {ux, uy, uz}, et);
                scalarsMap["Total Temperature"](i, j, k) = _fluid.computeTotalTemperature_rho_u_et(rho, {ux, uy, uz}, et);
                scalarsMap["Entropy"](i, j, k) = _fluid.computeEntropy_rho_u_et(rho, {ux, uy, uz}, et);

                if (isBFMActive){
                    FloatType omega = _mesh.getInputFields(FieldNames::RPM)(i, j, k) * 2.0 * M_PI / 60.0;
                    FloatType radius = std::sqrt(_mesh.getVertex(i, j, k).z() * _mesh.getVertex(i, j, k).z() + _mesh.getVertex(i, j, k).y() * _mesh.getVertex(i, j, k).y());
                    FloatType theta = std::atan2(_mesh.getVertex(i, j, k).z(), _mesh.getVertex(i, j, k).y());
                    Vector3D gridVelocityCylindrical = {0.0, 0.0, omega * radius};
                    Vector3D gridVelocityCartesian = computeCartesianVectorFromCylindrical(gridVelocityCylindrical, theta);

                    scalarsMap["Grid Velocity X"](i, j, k) = gridVelocityCartesian.x();
                    scalarsMap["Grid Velocity Y"](i, j, k) = gridVelocityCartesian.y();
                    scalarsMap["Grid Velocity Z"](i, j, k) = gridVelocityCartesian.z();

                    scalarsMap["Relative Velocity X"](i, j, k) = scalarsMap["Velocity X"](i, j, k) - gridVelocityCartesian.x();
                    scalarsMap["Relative Velocity Y"](i, j, k) = scalarsMap["Velocity Y"](i, j, k) - gridVelocityCartesian.y();
                    scalarsMap["Relative Velocity Z"](i, j, k) = scalarsMap["Velocity Z"](i, j, k) - gridVelocityCartesian.z();

                    scalarsMap["Relative Mach"](i, j, k) = _fluid.computeMachNumber_rho_u_et(rho, {scalarsMap["Relative Velocity X"](i, j, k), scalarsMap["Relative Velocity Y"](i, j, k), scalarsMap["Relative Velocity Z"](i, j, k)}, et);
                    scalarsMap["Viscous Body Force X"](i, j, k) = _viscousForce(i, j, k).x();
                    scalarsMap["Viscous Body Force Y"](i, j, k) = _viscousForce(i, j, k).y();
                    scalarsMap["Viscous Body Force Z"](i, j, k) = _viscousForce(i, j, k).z();
                    scalarsMap["Inviscid Body Force X"](i, j, k) = _inviscidForce(i, j, k).x();
                    scalarsMap["Inviscid Body Force Y"](i, j, k) = _inviscidForce(i, j, k).y();
                    scalarsMap["Inviscid Body Force Z"](i, j, k) = _inviscidForce(i, j, k).z();

                    scalarsMap["Deviation Angle"](i, j, k) = _deviationAngle(i, j, k);
                }

            }
        }
    }
}


std::string COutputBase::getOutputFilename(size_t iterationCounter) {
    std::string filename;
    if (_isUnsteadyOutput){
        std::ostringstream oss;
        oss << _config.getSolutionName() << "_" << std::setw(6) << std::setfill('0') << iterationCounter;
        filename = oss.str();
    }
    else {
        filename = _config.getSolutionName();
    }

    return filename;
}