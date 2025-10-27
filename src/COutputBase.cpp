#include "COutputBase.hpp"
#include "commonFunctions.hpp"


COutputBase::COutputBase(const Config &config, const CMesh &mesh, const FlowSolution &solution, const CFluidBase &fluid, 
                        const Matrix3D<Vector3D> &inviscidForce, const Matrix3D<Vector3D> &viscousForce,
                        const Matrix3D<FloatType> &deviationAngle)

    : _config(config), _mesh(mesh), _solution(solution), _fluid(fluid), _inviscidForce(inviscidForce), _viscousForce(viscousForce), _deviationAngle(deviationAngle) {
        _isUnsteadyOutput = _config.saveUnsteadySolution();
        std::filesystem::create_directory(_outputDirectory);
    }



void COutputBase::getScalarFieldsMap(std::map<std::string, Matrix3D<FloatType>>& scalarsMap, bool alsoGradients) const {
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


    // allocate space also for gradients if needed
    if (alsoGradients){
        scalarsMap.emplace("Density Gradient X",          Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Density Gradient Y",          Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Density Gradient Z",          Matrix3D<FloatType>(ni, nj, nk));

        scalarsMap.emplace("Velocity X Gradient X",       Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Velocity X Gradient Y",       Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Velocity X Gradient Z",       Matrix3D<FloatType>(ni, nj, nk));

        scalarsMap.emplace("Velocity Y Gradient X",       Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Velocity Y Gradient Y",       Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Velocity Y Gradient Z",       Matrix3D<FloatType>(ni, nj, nk));

        scalarsMap.emplace("Velocity Z Gradient X",       Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Velocity Z Gradient Y",       Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Velocity Z Gradient Z",       Matrix3D<FloatType>(ni, nj, nk));

        scalarsMap.emplace("Pressure Gradient X",     Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Pressure Gradient Y",     Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Pressure Gradient Z",     Matrix3D<FloatType>(ni, nj, nk));
    }


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

    // compute gradients
    Matrix3D<Vector3D> rhoGrad(ni, nj, nk), velXGrad(ni, nj, nk), velYGrad(ni, nj, nk), velZGrad(ni, nj, nk), pressGrad(ni, nj, nk);
    if (alsoGradients){
        computeGradientGreenGauss(  _mesh.getSurfacesI(), _mesh.getSurfacesJ(), _mesh.getSurfacesK(), 
                                    _mesh.getMidPointsI(), _mesh.getMidPointsJ(), _mesh.getMidPointsK(), 
                                    _mesh.getVertices(), _mesh.getVolumes(), scalarsMap["Density"], rhoGrad);
        
        computeGradientGreenGauss(  _mesh.getSurfacesI(), _mesh.getSurfacesJ(), _mesh.getSurfacesK(), 
                                    _mesh.getMidPointsI(), _mesh.getMidPointsJ(), _mesh.getMidPointsK(), 
                                    _mesh.getVertices(), _mesh.getVolumes(), scalarsMap["Velocity X"], velXGrad);
        
        computeGradientGreenGauss(  _mesh.getSurfacesI(), _mesh.getSurfacesJ(), _mesh.getSurfacesK(), 
                                    _mesh.getMidPointsI(), _mesh.getMidPointsJ(), _mesh.getMidPointsK(), 
                                    _mesh.getVertices(), _mesh.getVolumes(), scalarsMap["Velocity Y"], velYGrad);
        
        computeGradientGreenGauss(  _mesh.getSurfacesI(), _mesh.getSurfacesJ(), _mesh.getSurfacesK(), 
                                    _mesh.getMidPointsI(), _mesh.getMidPointsJ(), _mesh.getMidPointsK(), 
                                    _mesh.getVertices(), _mesh.getVolumes(), scalarsMap["Velocity Z"], velZGrad);
        
        computeGradientGreenGauss(  _mesh.getSurfacesI(), _mesh.getSurfacesJ(), _mesh.getSurfacesK(), 
                                    _mesh.getMidPointsI(), _mesh.getMidPointsJ(), _mesh.getMidPointsK(), 
                                    _mesh.getVertices(), _mesh.getVolumes(), scalarsMap["Pressure"], pressGrad);
    }

    // plug them in the map
    for (size_t i = 0; i < ni; ++i) {
        for (size_t j = 0; j < nj; ++j) {
            for (size_t k = 0; k < nk; ++k) {
                scalarsMap["Density Gradient X"](i, j, k) = rhoGrad(i, j, k).x();
                scalarsMap["Density Gradient Y"](i, j, k) = rhoGrad(i, j, k).y();
                scalarsMap["Density Gradient Z"](i, j, k) = rhoGrad(i, j, k).z();

                scalarsMap["Velocity X Gradient X"](i, j, k) = velXGrad(i, j, k).x();
                scalarsMap["Velocity X Gradient Y"](i, j, k) = velXGrad(i, j, k).y();
                scalarsMap["Velocity X Gradient Z"](i, j, k) = velXGrad(i, j, k).z();

                scalarsMap["Velocity Y Gradient X"](i, j, k) = velYGrad(i, j, k).x();
                scalarsMap["Velocity Y Gradient Y"](i, j, k) = velYGrad(i, j, k).y();
                scalarsMap["Velocity Y Gradient Z"](i, j, k) = velYGrad(i, j, k).z();

                scalarsMap["Velocity Z Gradient X"](i, j, k) = velZGrad(i, j, k).x();
                scalarsMap["Velocity Z Gradient Y"](i, j, k) = velZGrad(i, j, k).y();
                scalarsMap["Velocity Z Gradient Z"](i, j, k) = velZGrad(i, j, k).z();

                scalarsMap["Pressure Gradient X"](i, j, k) = pressGrad(i, j, k).x();
                scalarsMap["Pressure Gradient Y"](i, j, k) = pressGrad(i, j, k).y();
                scalarsMap["Pressure Gradient Z"](i, j, k) = pressGrad(i, j, k).z();
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