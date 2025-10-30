#include "COutputBase.hpp"
#include "commonFunctions.hpp"


COutputBase::COutputBase(const Config &config, const CMesh &mesh, const FlowSolution &solution, const CFluidBase &fluid, 
                        const Matrix3D<Vector3D> &inviscidForce, const Matrix3D<Vector3D> &viscousForce,
                        const Matrix3D<FloatType> &deviationAngle)

    : _config(config), _mesh(mesh), _solution(solution), _fluid(fluid), _inviscidForce(inviscidForce), _viscousForce(viscousForce), _deviationAngle(deviationAngle) {
        
        _isUnsteadyOutput = _config.saveUnsteadySolution();
        
        std::filesystem::create_directory(_outputDirectory);

        _outputFields = _config.getOutputFields();
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
    if (_outputFields == OutputFields::SECONDARY || _outputFields == OutputFields::TURBO_BFM){
        scalarsMap.emplace("Pressure",          Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Temperature",       Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Mach",              Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Total Pressure",    Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Total Temperature", Matrix3D<FloatType>(ni, nj, nk));
        scalarsMap.emplace("Entropy",           Matrix3D<FloatType>(ni, nj, nk));
    }
    
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


    
    // BFM simulation additional variables
    const bool isBFMActive = _config.isBFMActive();
    if (isBFMActive && _outputFields == OutputFields::TURBO_BFM){
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

    

    if (_outputFields == OutputFields::SECONDARY || _outputFields == OutputFields::TURBO_BFM){
        Vector3D vel, gridVelCyl, gridVelCart, relVel;
        FloatType rho, et, omega, radius, theta;
        
        for (size_t i = 0; i < ni; ++i) {
            for (size_t j = 0; j < nj; ++j) {
                for (size_t k = 0; k < nk; ++k) {
                
                    rho = scalarsMap["Density"](i, j, k);
                    vel.x() = scalarsMap["Velocity X"](i, j, k);
                    vel.y() = scalarsMap["Velocity Y"](i, j, k);
                    vel.z() = scalarsMap["Velocity Z"](i, j, k);
                    et  = scalarsMap["Total Energy"](i, j, k);

                    scalarsMap["Pressure"](i, j, k) = _fluid.computePressure_rho_u_et(rho, vel, et);
                    scalarsMap["Temperature"](i, j, k) = _fluid.computeTemperature_rho_u_et(rho, vel, et);
                    scalarsMap["Mach"](i, j, k) = _fluid.computeMachNumber_rho_u_et(rho, vel, et);
                    scalarsMap["Total Pressure"](i, j, k) = _fluid.computeTotalPressure_rho_u_et(rho, vel, et);
                    scalarsMap["Total Temperature"](i, j, k) = _fluid.computeTotalTemperature_rho_u_et(rho, vel, et);
                    scalarsMap["Entropy"](i, j, k) = _fluid.computeEntropy_rho_u_et(rho, vel, et);

                    if (isBFMActive && _outputFields == OutputFields::TURBO_BFM){
                        omega = _mesh.getInputFields(FieldNames::RPM)(i, j, k) * 2.0 * M_PI / 60.0;
                        radius = std::sqrt(_mesh.getVertex(i, j, k).z() * _mesh.getVertex(i, j, k).z() +
                                                    _mesh.getVertex(i, j, k).y() * _mesh.getVertex(i, j, k).y());
                        theta = std::atan2(_mesh.getVertex(i, j, k).z(), _mesh.getVertex(i, j, k).y());
                        gridVelCyl = {0.0, 0.0, omega * radius};
                        gridVelCart = computeCartesianVectorFromCylindrical(gridVelCyl, theta);

                        scalarsMap["Grid Velocity X"](i, j, k) = gridVelCart.x();
                        scalarsMap["Grid Velocity Y"](i, j, k) = gridVelCart.y();
                        scalarsMap["Grid Velocity Z"](i, j, k) = gridVelCart.z();
                        
                        relVel = vel - gridVelCart;
                        scalarsMap["Relative Velocity X"](i, j, k) = relVel.x();
                        scalarsMap["Relative Velocity Y"](i, j, k) = relVel.y();
                        scalarsMap["Relative Velocity Z"](i, j, k) = relVel.z();

                        scalarsMap["Relative Mach"](i, j, k) = _fluid.computeMachNumber_rho_u_et(rho, relVel, et);

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


    

    
    if (alsoGradients){
        // compute gradients
        Matrix3D<Vector3D> rhoGrad(ni, nj, nk), velXGrad(ni, nj, nk), velYGrad(ni, nj, nk), velZGrad(ni, nj, nk), pressGrad(ni, nj, nk);
        
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