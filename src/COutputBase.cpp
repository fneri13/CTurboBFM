#include "COutputBase.hpp"
#include "commonFunctions.hpp"


COutputBase::COutputBase(const Config &config, const CMesh &mesh, const FlowSolution &solution, const CFluid &fluid, 
                        const Matrix3D<Vector3D> &inviscidForce, const Matrix3D<Vector3D> &viscousForce,
                        const Matrix3D<FloatType> &deviationAngle)

    : _config(config), _mesh(mesh), _solution(solution), _fluid(fluid), _inviscidForce(inviscidForce), _viscousForce(viscousForce), _deviationAngle(deviationAngle) {
        _isUnsteadyOutput = _config.saveUnsteadySolution();
        std::filesystem::create_directory(_outputDirectory);
    }



void COutputBase::getScalarFieldsMap(std::map<std::string, Matrix3D<FloatType>>& scalarsMap) const {
    size_t ni = _mesh.getNumberPointsI();
    size_t nj = _mesh.getNumberPointsJ();
    size_t nk = _mesh.getNumberPointsK();
    
    // primary solution variables
    scalarsMap["Density"] = _solution.getDensity();
    scalarsMap["Velocity X"] = _solution.getVelocityX();
    scalarsMap["Velocity Y"] = _solution.getVelocityY();
    scalarsMap["Velocity Z"] = _solution.getVelocityZ();
    scalarsMap["Total Energy"] = _solution.getTotalEnergy();
    
    // Auxiliary solution variables
    scalarsMap.emplace("Pressure",          Matrix3D<FloatType>(ni, nj, nk));
    scalarsMap.emplace("Temperature",       Matrix3D<FloatType>(ni, nj, nk));
    scalarsMap.emplace("Mach",              Matrix3D<FloatType>(ni, nj, nk));
    scalarsMap.emplace("Total Pressure",    Matrix3D<FloatType>(ni, nj, nk));
    scalarsMap.emplace("Total Temperature", Matrix3D<FloatType>(ni, nj, nk));
    scalarsMap.emplace("Entropy",           Matrix3D<FloatType>(ni, nj, nk));

    // Cache references for speed
    auto& density     = scalarsMap["Density"];
    auto& velX        = scalarsMap["Velocity X"];
    auto& velY        = scalarsMap["Velocity Y"];
    auto& velZ        = scalarsMap["Velocity Z"];
    auto& etot        = scalarsMap["Total Energy"];
    auto& pressure    = scalarsMap["Pressure"];
    auto& temperature = scalarsMap["Temperature"];
    auto& mach        = scalarsMap["Mach"];
    auto& pTotal      = scalarsMap["Total Pressure"];
    auto& tTotal      = scalarsMap["Total Temperature"];
    auto& entropy     = scalarsMap["Entropy"];

    // BFM simulation additional variables
    const bool isBFMActive = _config.isBFMActive();
    if (isBFMActive){
        scalarsMap["Blockage"] = _mesh.getInputFields(FieldNames::BLOCKAGE);
        scalarsMap.emplace("Relative Mach",          Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Grid Velocity X",        Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Grid Velocity Y",        Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Grid Velocity Z",        Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Relative Velocity X",    Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Relative Velocity Y",    Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Relative Velocity Z",    Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Viscous Body Force X",   Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Viscous Body Force Y",   Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Viscous Body Force Z",   Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Inviscid Body Force X",  Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Inviscid Body Force Y",  Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Inviscid Body Force Z",  Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Deviation Angle",        Matrix3D<FloatType>(ni, nj, nk));

    }

    

    // Reusable temporaries
    Vector3D vel;
    FloatType rho, et;

    for (size_t i = 0; i < ni; ++i) {
        for (size_t j = 0; j < nj; ++j) {
            for (size_t k = 0; k < nk; ++k) {
                rho = density(i, j, k);
                vel.x() = velX(i, j, k);
                vel.y() = velY(i, j, k);
                vel.z() = velZ(i, j, k);
                et  = etot(i, j, k);

                pressure(i, j, k)    = _fluid.computePressure_rho_u_et(rho, vel, et);
                temperature(i, j, k) = _fluid.computeTemperature_rho_u_et(rho, vel, et);
                mach(i, j, k)        = _fluid.computeMachNumber_rho_u_et(rho, vel, et);
                pTotal(i, j, k)      = _fluid.computeTotalPressure_rho_u_et(rho, vel, et);
                tTotal(i, j, k)      = _fluid.computeTotalTemperature_rho_u_et(rho, vel, et);
                entropy(i, j, k)     = _fluid.computeEntropy_rho_u_et(rho, vel, et);

                if (isBFMActive){
                    auto& gvx = scalarsMap["Grid Velocity X"];
                    auto& gvy = scalarsMap["Grid Velocity Y"];
                    auto& gvz = scalarsMap["Grid Velocity Z"];
                    auto& rvx = scalarsMap["Relative Velocity X"];
                    auto& rvy = scalarsMap["Relative Velocity Y"];
                    auto& rvz = scalarsMap["Relative Velocity Z"];
                    auto& rmach = scalarsMap["Relative Mach"];
                    auto& vbfX = scalarsMap["Viscous Body Force X"];
                    auto& vbfY = scalarsMap["Viscous Body Force Y"];
                    auto& vbfZ = scalarsMap["Viscous Body Force Z"];
                    auto& ibfX = scalarsMap["Inviscid Body Force X"];
                    auto& ibfY = scalarsMap["Inviscid Body Force Y"];
                    auto& ibfZ = scalarsMap["Inviscid Body Force Z"];
                    auto& devAng = scalarsMap["Deviation Angle"];

                    FloatType omega = _mesh.getInputFields(FieldNames::RPM)(i, j, k) * 2.0 * M_PI / 60.0;
                    FloatType radius = std::sqrt(_mesh.getVertex(i, j, k).z() * _mesh.getVertex(i, j, k).z() +
                                                 _mesh.getVertex(i, j, k).y() * _mesh.getVertex(i, j, k).y());
                    FloatType theta = std::atan2(_mesh.getVertex(i, j, k).z(), _mesh.getVertex(i, j, k).y());
                    Vector3D gridVelCyl = {0.0, 0.0, omega * radius};
                    Vector3D gridVelCart = computeCartesianVectorFromCylindrical(gridVelCyl, theta);

                    gvx(i, j, k) = gridVelCart.x();
                    gvy(i, j, k) = gridVelCart.y();
                    gvz(i, j, k) = gridVelCart.z();

                    rvx(i, j, k) = velX(i, j, k) - gridVelCart.x();
                    rvy(i, j, k) = velY(i, j, k) - gridVelCart.y();
                    rvz(i, j, k) = velZ(i, j, k) - gridVelCart.z();

                    vel.x() = rvx(i, j, k);
                    vel.y() = rvy(i, j, k);
                    vel.z() = rvz(i, j, k);
                    rmach(i, j, k) = _fluid.computeMachNumber_rho_u_et(rho, vel, et);

                    vbfX(i, j, k) = _viscousForce(i, j, k).x();
                    vbfY(i, j, k) = _viscousForce(i, j, k).y();
                    vbfZ(i, j, k) = _viscousForce(i, j, k).z();
                    ibfX(i, j, k) = _inviscidForce(i, j, k).x();
                    ibfY(i, j, k) = _inviscidForce(i, j, k).y();
                    ibfZ(i, j, k) = _inviscidForce(i, j, k).z();

                    devAng(i, j, k) = _deviationAngle(i, j, k);
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