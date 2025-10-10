#pragma once
#include "types.hpp"
#include "CoolProp.h"

// Class for ideal gases
class CFluidBase {
public:

    /** Constructor for base class
    */
    CFluidBase(){};
    
    virtual ~CFluidBase() = default;

    /** Compute static energy from pressure and density
     * @param p Pressure [Pa]
     * @param rho Density [kg/m³]
     * @return Static energy [J/kg]
    */
    virtual FloatType computeStaticEnergy_p_rho(FloatType p, FloatType rho) const = 0;
    
    /** Compute pressure from density and static energy
     * @param rho Density [kg/m³]
     * @param e Static energy [J/kg]
     * @return Pressure [Pa]
    */
    virtual FloatType computePressure_rho_e(FloatType rho, FloatType e) const = 0;

    /** Compute pressure from density and temperature
     * @param rho Density [kg/m³]
     * @param e Static energy [J/kg]
     * @return Pressure [Pa]
    */
    virtual FloatType computePressure_rho_T(FloatType rho, FloatType Temp) const = 0;
    
    /** Compute speed of sound from pressure and density
     * @param p Pressure [Pa]
     * @param rho Density [kg/m³]
     * @return Speed of sound [m/s]
    */
    virtual FloatType computeSoundSpeed_p_rho(FloatType p, FloatType rho) const = 0;
    
    /** Compute static energy from velocity magnitude and total energy
     * @param velMag Velocity magnitude [m/s]
     * @param et Total energy [J/kg]
     * @return Static energy [J/kg]
    */
    inline FloatType computeStaticEnergy_u_et(FloatType velMag, FloatType et) const {
        return et - 0.5 * velMag * velMag;
    };

    /** Compute speed of sound from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Speed of sound [m/s]
    */
    virtual FloatType computeSoundSpeed_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;
    
    /** Compute pressure from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Pressure [Pa]
    */
    virtual FloatType computePressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;
    
    /** Compute total pressure from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Total pressure [Pa]
    */
    virtual FloatType computeTotalPressure_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;
    
    /** Compute total temperature from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Total temperature [K]
    */
    virtual FloatType computeTotalTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;

    /** Compute static temperature from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Static temperature [K]
    */
    virtual FloatType computeTemperature_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;
    
    /** Compute total enthalpy from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Total enthalpy [J/kg]
    */
    virtual FloatType computeTotalEnthalpy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;

    /** Compute static pressure from total pressure and Mach number
     * @param pt Total pressure [Pa]
     * @param M Mach number [-]
     * @return Static pressure [Pa]
    */
    virtual FloatType computeStaticPressure_pt_M(FloatType pt, FloatType M) const = 0;
        
    /** Compute static temperature from total temperature and Mach number
     * @param Tt Total temperature [K]
     * @param M Mach number [-]
     * @return Static temperature [K]
    */
    virtual FloatType computeStaticTemperature_Tt_M(FloatType Tt, FloatType M) const = 0;
    
    /** Compute entropy from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Entropy [J/kgK]
    */
    virtual FloatType computeEntropy_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;

    /** Compute density from pressure and temperature
     * @param p Pressure [Pa]
     * @param T Temperature [K]
     * @return Density [kg/m³]
    */
    virtual FloatType computeDensity_p_T(FloatType p, FloatType T) const = 0;

    /** Compute Mach number from density, velocity vector, and total energy
     * @param rho Density [kg/m³]
     * @param u Velocity vector [m/s]
     * @param et Total energy [J/kg]
     * @return Mach number [-]
    */
    virtual FloatType computeMachNumber_rho_u_et(FloatType rho, const Vector3D& u, FloatType et) const = 0;
    
    /** Compute total pressure from static pressure and Mach number
     * @param pressure Static pressure [Pa]
     * @param mach Mach number [-]
     * @return Total pressure [Pa]
    */
    virtual FloatType computeTotalPressure_p_M(FloatType pressure, FloatType mach) const = 0;
    
    /** Compute total temperature from static temperature and Mach number
     * @param temperature Static temperature [K]
     * @param mach Mach number [-]
     * @return Total temperature [K]
    */
    virtual FloatType computeTotalTemperature_T_M(FloatType temperature, FloatType mach) const = 0;

    /** Compute entropy from pressure and density
     * @param pressure Pressure [Pa]
     * @param density Density [kg/m³]
     * @return Entropy [J/kgK]
    */
    virtual FloatType computeEntropy_p_rho(FloatType pressure, FloatType density) const = 0;

    /** Compute entropy from pressure and density
     * @param pressure Pressure [Pa]
     * @param density Temperature [J/kgK]
     * @return Entropy [J/kgK]
    */
    virtual FloatType computeEntropy_p_T(FloatType pressure, FloatType temperature) const = 0;

    /** Compute initial primitive variables from given freestream conditions
     * @param initMach Freestream Mach number [-]
     * @param initTemperature Freestream static temperature [K]
     * @param initPressure Freestream static pressure [Pa]
     * @param flowDirection Array of unit vector for flow direction [-]
     * @param density Output: initial density [kg/m³]
     * @param velocity Output: initial velocity vector [m/s]
     * @param totEnergy Output: initial total energy [J/kg]
    */
    virtual void computeInitFields(FloatType initMach, FloatType initTemperature, FloatType initPressure, Vector3D flowDirection, FloatType &density, Vector3D &velocity, FloatType &totEnergy) = 0;

    /** Compute pressure from primitive state vector
     * @param primitive State vector (rho, u, v, w, et) [SI units]
     * @return Pressure [Pa]
    */
    virtual FloatType computePressure_primitive(StateVector primitive) const = 0;

    /** Compute total efficiency from total pressure and temperature ratios
     * @param pressureRatio Total pressure ratio P0_exit / P0_inlet [-]
     * @param temperatureRatio Total temperature ratio T0_exit / T0_inlet [-]
     * @return Total efficiency [-]
    */
    virtual FloatType computeTotalEfficiency_PRtt_TRt(FloatType pressureRatio, FloatType temperatureRatio) const = 0;

    /** Get specific heat ratio gamma
     * @return Gamma [-]
    */
    virtual FloatType getGamma() const { return 1.4; };

private:
    
};
