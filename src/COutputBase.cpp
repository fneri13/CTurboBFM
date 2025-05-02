#include "COutputBase.hpp"

void COutputBase::getScalarFieldsMap(std::map<std::string, Matrix3D<FloatType>>& scalarsMap) {

    scalarsMap["Density"] = _solution.getDensity();
    scalarsMap["Velocity X"] = _solution.getVelocityX();
    scalarsMap["Velocity Y"] = _solution.getVelocityY();
    scalarsMap["Velocity Z"] = _solution.getVelocityZ();
    scalarsMap["Total Energy"] = _solution.getTotalEnergy();
    
    size_t ni = _mesh.getNumberPointsI();
    size_t nj = _mesh.getNumberPointsJ();
    size_t nk = _mesh.getNumberPointsK();
    scalarsMap["Pressure"] = Matrix3D<FloatType>(ni, nj, nk);
    scalarsMap["Temperature"] = Matrix3D<FloatType>(ni, nj, nk);

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
            }
        }
    }
}
