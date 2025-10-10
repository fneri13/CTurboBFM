#pragma once
#include "types.hpp"
#include "CFluidBase.hpp"
#include "CoolProp.h"
#include "AbstractState.h"
#include <memory>
#include <iostream>


// Class for ideal gases
class CFluidReal : public CFluidBase {
public:

    /** Constructor for real gas
     * @param fluidName Name of the fluid according to CoolProp
    */
    CFluidReal(std::string fluidName);

    /** Compute static energy from pressure and density
     * @param p Pressure [Pa]
     * @param rho Density [kg/m³]
     * @return Static energy [J/kg]
    */
    FloatType computeStaticEnergy_p_rho(FloatType p, FloatType rho) const override;
    
    /** Compute pressure from density and static energy
     * @param rho Density [kg/m³]
     * @param e Static energy [J/kg]
     * @return Pressure [Pa]
    */
    FloatType computePressure_rho_e(FloatType rho, FloatType e) const override;

    /** Compute pressure from density and temperature
     * @param rho Density [kg/m³]
     * @param e Static energy [J/kg]
     * @return Pressure [Pa]
    */
    FloatType computePressure_rho_T(FloatType rho, FloatType Temp) const override;
    
    /** Compute speed of sound from pressure and density
     * @param p Pressure [Pa]
     * @param rho Density [kg/m³]
     * @return Speed of sound [m/s]
    */
    FloatType computeSoundSpeed_p_rho(FloatType p, FloatType rho) const override;

    /** Compute speed of sound from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Speed of sound [m/s]
    */
    FloatType computeSoundSpeed_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;
    
    /** Compute pressure from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Pressure [Pa]
    */
    FloatType computePressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;
    
    /** Compute total pressure from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Total pressure [Pa]
    */
    FloatType computeTotalPressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;
    
    /** Compute total temperature from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Total temperature [K]
    */
    FloatType computeTotalTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;

    /** Compute static temperature from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Static temperature [K]
    */
    FloatType computeTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;
    
    /** Compute total enthalpy from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Total enthalpy [J/kg]
    */
    FloatType computeTotalEnthalpy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;

    /** Compute static pressure from total pressure and Mach number
     * @param pt Total pressure [Pa]
     * @param M Mach number [-]
     * @return Static pressure [Pa]
    */
    FloatType computeStaticPressure_pt_M(FloatType pt, FloatType M) const override;
        
    /** Compute static temperature from total temperature and Mach number
     * @param Tt Total temperature [K]
     * @param M Mach number [-]
     * @return Static temperature [K]
    */
    FloatType computeStaticTemperature_Tt_M(FloatType Tt, FloatType M) const override;
    
    /** Compute entropy from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Entropy [J/kgK]
    */
    FloatType computeEntropy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;

    /** Compute density from pressure and temperature
     * @param p Pressure [Pa]
     * @param T Temperature [K]
     * @return Density [kg/m³]
    */
    FloatType computeDensity_p_T(FloatType p, FloatType T) const override;

    /** Compute Mach number from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Mach number [-]
    */
    FloatType computeMachNumber_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const override;
    
    /** Compute total pressure from static pressure and Mach number
     * @param pressure Static pressure [Pa]
     * @param mach Mach number [-]
     * @return Total pressure [Pa]
    */
    FloatType computeTotalPressure_p_M(FloatType pressure, FloatType mach) const override;
    
    /** Compute total temperature from static temperature and Mach number
     * @param temperature Static temperature [K]
     * @param mach Mach number [-]
     * @return Total temperature [K]
    */
    FloatType computeTotalTemperature_T_M(FloatType temperature, FloatType mach) const override;

    /** Compute entropy from pressure and density
     * @param pressure Pressure [Pa]
     * @param density Density [kg/m³]
     * @return Entropy [J/kgK]
    */
    FloatType computeEntropy_p_rho(FloatType pressure, FloatType density) const override;

    /** Compute entropy from pressure and density
     * @param pressure Pressure [Pa]
     * @param density Temperature [J/kgK]
     * @return Entropy [J/kgK]
    */
    FloatType computeEntropy_p_T(FloatType pressure, FloatType temperature) const override;

    /** Compute initial primitive variables from given freestream conditions
     * @param initMach Freestream Mach number [-]
     * @param initTemperature Freestream static temperature [K]
     * @param initPressure Freestream static pressure [Pa]
     * @param flowDirection Array of unit vector for flow direction [-]
     * @param density Output: initial density [kg/m³]
     * @param velocity Output: initial velocity vector [m/s]
     * @param totEnergy Output: initial total energy [J/kg]
    */
    void computeInitFields(FloatType initMach, FloatType initTemperature, FloatType initPressure, Vector3D flowDirection, FloatType &density, Vector3D &velocity, FloatType &totEnergy) override;

    /** Compute pressure from primitive state vector
     * @param primitive State vector (rho, u, v, w, et) [SI units]
     * @return Pressure [Pa]
    */
    FloatType computePressure_primitive(StateVector primitive) const override;

    /** Compute total efficiency from total pressure and temperature ratios
     * @param pressureRatio Total pressure ratio P0_exit / P0_inlet [-]
     * @param temperatureRatio Total temperature ratio T0_exit / T0_inlet [-]
     * @return Total efficiency [-]
    */
    FloatType computeTotalEfficiency_PRtt_TRt(FloatType pressureRatio, FloatType temperatureRatio) const override;

private:
    std::string _fluidName; ///< Name of the fluid according to CoolProp
};
