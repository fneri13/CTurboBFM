#include "COutputBase.hpp"
#include "commonFunctions.hpp"

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
                }

            }
        }
    }
}
