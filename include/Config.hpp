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

    int getSaveUnsteadyInterval() const {return static_cast<int>(parseFloat("SAVE_UNSTEADY_INTERVAL"));}

    std::string getSolutionName() const {return parseString("SOLUTION_NAME");}

    std::string getRestartFilepath() const {return parseString("RESTART_SOLUTION_FILEPATH");}

    std::string getRestartType() const {return parseString("RESTART_TYPE", true);}

    bool saveTurboOutput() const {return parseBool("TURBO_OUTPUT", false);}

    TimeIntegration getTimeIntegration() const;

    TimeStepMethod getTimeStepMethod() const;

    ConvectionScheme getConvectionScheme() const;

    // get the coefficients for the Runge-Kutta time integration, taken from Simon thesis (page 64)
    std::vector<FloatType> getTimeIntegrationCoeffs() const ;

    size_t getMassFlowUpdateFrequency() const {return massFlowUpdateFrequency;}

    size_t getSolutionOutputFrequency() const {return getSaveUnsteadyInterval();}

    FloatType getPeriodicityAngleRad() const {return parseFloat("PERIODICITY_ANGLE")*M_PI/180.0;} 

    FloatType getPeriodicityAngleDeg() const {return parseFloat("PERIODICITY_ANGLE");} 

    void printAllConfigValues() const;

};
