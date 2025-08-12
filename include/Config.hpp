#pragma once
#include <string>
#include <map>
#include <sstream>
#include "types.hpp"

class Config {

private:

    size_t massFlowUpdateFrequency = {100};

    std::map<std::string, std::string> _data;
    
    bool load(const std::string& filename);
    
    bool parseBool(const std::string& value, bool defaultValue) const;
    

    /**
     * @brief Return the string file corresponding to a certain entry.
     * @param value option to read from config file.
     * @param defaultNoneValue if true return None string if the option is not present, otherwise throw an error
     */
    std::string parseString(const std::string& value, bool defaultNoneValue=false) const;
    
    FloatType parseFloat(const std::string& value) const;

    FloatType parseFloat(const std::string& value, FloatType defaultValue) const;

    int parseInt(const std::string& value) const;

    int parseInt(const std::string& value, int defaultValue) const;

    template <typename T>
    std::vector<T> parseVector(const std::string& key) const {
        if (!has(key)) {
            throw std::runtime_error("Missing value for key \"" + key + "\" in configuration.");
        }

        std::string value = get(key);
        std::vector<T> result;
        std::stringstream ss(value);
        std::string item;

        while (std::getline(ss, item, ',')) {
            // Trim leading/trailing whitespace
            item.erase(0, item.find_first_not_of(" \t\r\n"));
            item.erase(item.find_last_not_of(" \t\r\n") + 1);

            if (item.empty()) continue;

            std::stringstream itemStream(item);
            T val;
            itemStream >> val;

            if (itemStream.fail() || !itemStream.eof()) {
                throw std::runtime_error("Invalid vector element for key \"" + key + "\": \"" + item + "\"");
            }

            result.push_back(val);
        }

        return result;
    }

public:
    Config () = default;
    Config (const std::string& filename);
    ~Config () = default;
    
    std::string get(const std::string& key) const;
    
    bool has(const std::string& key) const;
    
    std::string gridFilePath() const;

    bool restartSolution() const;
    
    std::string restartSolutionFilepath() const;
    
    bool isBFMActive() const;
    
    bool isBlockageActive() const;
    
    BFM_Model getBFMModel() const;

    KindSolver getKindSolver() const;

    Topology getTopology() const;

    FloatType getFluidGamma() const {return parseFloat("FLUID_GAMMA");} 

    FloatType getFluidGasConstant() const {return parseFloat("FLUID_R_CONSTANT");} 
    
    FloatType getFluidKinematicViscosity() const {return parseFloat("FLUID_KINEMATIC_VISCOSITY");} 

    BoundaryType getBoundaryType(const BoundaryIndices) const;

    std::vector<FloatType> getInletBCValues() const { return parseVector<FloatType>("INLET_VALUE");}

    std::vector<FloatType> getOutletBCValues() const { return parseVector<FloatType>("OUTLET_VALUE");}

    FloatType getInitMachNumber() const { return parseFloat("INIT_MACH_NUMBER");}
    
    FloatType getInitPressure() const { return parseFloat("INIT_PRESSURE");}
    
    FloatType getInitTemperature() const { return parseFloat("INIT_TEMPERATURE");}
    
    Vector3D getInitDirection() const ;
    
    FloatType getCFL() const {return parseFloat("CFL");}

    long unsigned int getMaxIterations() const {return static_cast<long unsigned int>(parseFloat("N_ITERATIONS"));}

    bool saveUnsteadySolution() const {return parseBool("SAVE_UNSTEADY", true);}

    bool isSimulationSteady() const {return parseBool("SIMULATION_IS_STEADY", true);}

    int getSaveIterationsInterval() const {return static_cast<int>(parseFloat("SAVE_ITERATIONS_INTERVAL"));}

    std::string getSolutionName() const {return parseString("SOLUTION_NAME");}

    std::string getChimaScalingFunctionsFile() const {return parseString("CHIMA_SCALING_FUNCTIONS_FILEPATH");}

    FloatType getChimaReferenceMassFlow() const {return parseFloat("CHIMA_REFERENCE_MASS_FLOW");}

    std::string getRestartFilepath() const {return parseString("RESTART_SOLUTION_FILEPATH");}

    std::string getRestartType() const {return parseString("RESTART_TYPE", true);}

    bool saveTurboOutput() const {return parseBool("TURBO_OUTPUT", false);}

    bool isMonitorPointsActive() const {return parseBool("MONITOR_POINTS_ACTIVE", false);}

    size_t getResidualsDropConvergence() const {return static_cast<size_t>(parseInt("RESIDUALS_DROP_CONVERGENCE", 16));}

    TimeIntegration getTimeIntegration() const;

    TimeStepMethod getTimeStepMethod() const;

    ConvectionScheme getConvectionScheme() const;

    // get the coefficients for the Runge-Kutta time integration, taken from Simon thesis (page 64)
    std::vector<FloatType> getTimeIntegrationCoeffs() const ;

    size_t getMassFlowUpdateFrequency() const {return massFlowUpdateFrequency;}

    size_t getSolutionOutputFrequency() const {return getSaveIterationsInterval();}

    FloatType getPeriodicityAngleRad() const {return parseFloat("PERIODICITY_ANGLE")*M_PI/180.0;} 

    FloatType getPeriodicityAngleDeg() const {return parseFloat("PERIODICITY_ANGLE");} 

    void printAllConfigValues() const;

    FloatType computeRampedOutletPressure(const size_t iterCounter, const FloatType outletPressure) const;

    FloatType getOutletPressureRampMaxIterations() const {return parseFloat("OUTLET_PRESSURE_RAMP_ITERATIONS", 0.0);}

    FloatType getHallTholletCoefficient_KN() const {return parseFloat("HALL_THOLLET_COEFFICIENT_KN");} 

    FloatType getHallTholletCoefficient_KF() const {return parseFloat("HALL_THOLLET_COEFFICIENT_KF");} 

    FloatType getHallTholletCoefficient_KD() const {return parseFloat("HALL_THOLLET_COEFFICIENT_KD");} 

    std::string getHallTholletOffDesignActive() const {return parseString("HALL_THOLLET_OFF_DESIGN_ACTIVE");}

    int getMonitorPointsCoordsI() const { return parseInt("MONITOR_POINTS_I_COORDS");}

    int getMonitorPointsCoordsJ() const { return parseInt("MONITOR_POINTS_J_COORDS");}

    int getCircumferentialNumberMonitorPoints() const { return parseInt("MONITOR_POINTS_NUMBER");}

    bool isPerturbationBodyForceActive() const {return parseBool("PERTURBATION_BODY_FORCE", false);}

    bool getMUSCLreconstruction() const {return parseBool("MUSCL_RECONSTRUCTION", false);}

    FluxLimiter getFluxLimiter() const ;

    std::vector<FloatType> getPerturbationCenter() const {return parseVector<FloatType>("PERTURBATION_CENTER");}

    FloatType getPerturbationScalingFactor() const {return parseFloat("PERTURBATION_SCALING_FACTOR");} 

    FloatType getPerturbationTimeDuration() const {return parseFloat("PERTURBATION_TIME_DURATION");} 
    
    FloatType getPerturbationStartTime() const {return parseFloat("PERTURBATION_TIME_START");} 

    FloatType getPerturbationRadialExtension() const {return parseFloat("PERTURBATION_RADIAL_EXTENSION");} 

    FloatType getFixedTimeStep() const {return parseFloat("FIXED_TIME_STEP");} 
};
