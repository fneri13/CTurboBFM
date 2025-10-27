#include "Config.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

Config::Config(const std::string& filename) {
    if (!load(filename)) {
        throw std::runtime_error("Failed to load configuration file: " + filename);
    }
    printAllConfigValues();
}

bool Config::load(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::string line;
    while (std::getline(file, line)) {
        // Trim line
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        // Skip comments and empty lines
        if (line.empty() || line[0] == ';' || line[0] == '[') continue;

        auto delimiter_pos = line.find('=');
        if (delimiter_pos == std::string::npos) continue;

        std::string key = line.substr(0, delimiter_pos);
        std::string value = line.substr(delimiter_pos + 1);

        // Trim key and value
        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t") + 1);

        _data[key] = value;
    }

    return true;
}

std::string Config::get(const std::string& key) const {
    auto it = _data.find(key);
    if (it != _data.end()) {
        return it->second;
    }
    return "";
}

bool Config::has(const std::string& key) const {
    return _data.count(key) > 0;
}

string Config::parseString(const std::string& key, bool defaultNoneValue) const {
    if (!has(key)) {
        if (defaultNoneValue) {
            return "None";
        }
        else {
            throw std::runtime_error("Not found value for key \"" + key + "\" in configuration.");
        }
    }
    return get(key);
}

FloatType Config::parseFloat(const std::string& key) const {
    if (!has(key)) {
        throw std::runtime_error("Not found value for key \"" + key + "\" in configuration.");
    }
    std::string value = get(key);
    try {
        return static_cast<FloatType>(std::stod(value));  // or std::stof if FloatType is float
    } catch (const std::invalid_argument& e) {
        throw std::runtime_error("Invalid float format for key \"" + key + "\": \"" + value + "\"");
    } catch (const std::out_of_range& e) {
        throw std::runtime_error("Float value out of range for key \"" + key + "\": \"" + value + "\"");
    }
}

int Config::parseInt(const std::string& key) const {
    if (!has(key)) {
        throw std::runtime_error("Not found value for key \"" + key + "\" in configuration.");
    }
    std::string value = get(key);
    try {
        return static_cast<int>(std::stod(value));  // or std::stof if FloatType is float
    } catch (const std::invalid_argument& e) {
        throw std::runtime_error("Invalid float format for key \"" + key + "\": \"" + value + "\"");
    } catch (const std::out_of_range& e) {
        throw std::runtime_error("Float value out of range for key \"" + key + "\": \"" + value + "\"");
    }
}

int Config::parseInt(const std::string& key, int defaultValue) const {
    if (!has(key)) {
        return defaultValue;
    }
    std::string value = get(key);
    try {
        return static_cast<int>(std::stod(value));  // or std::stof if FloatType is float
    } catch (const std::invalid_argument& e) {
        throw std::runtime_error("Invalid float format for key \"" + key + "\": \"" + value + "\"");
    } catch (const std::out_of_range& e) {
        throw std::runtime_error("Float value out of range for key \"" + key + "\": \"" + value + "\"");
    }
}

FloatType Config::parseFloat(const std::string& key, FloatType defaultValue) const {
    if (!has(key)) {
        return defaultValue;
    }
    std::string value = get(key);
    try {
        return static_cast<FloatType>(std::stod(value));  // or std::stof if FloatType is float
    } catch (const std::invalid_argument& e) {
        throw std::runtime_error("Invalid float format for key \"" + key + "\": \"" + value + "\"");
    } catch (const std::out_of_range& e) {
        throw std::runtime_error("Float value out of range for key \"" + key + "\": \"" + value + "\"");
    }
}

bool Config::parseBool(const std::string& key, bool defaultValue) const {
    if (!has(key)) {
        return defaultValue;
    } else {
        std::string value = get(key);
        
        // Make it case-insensitive
        std::transform(value.begin(), value.end(), value.begin(),
                       [](unsigned char c){ return std::tolower(c); });

        if (value == "true" || value == "yes" || value == "1") {
            return true;
        } else if (value == "false" || value == "no" || value == "0") {
            return false;
        } else {
            throw std::runtime_error("Invalid value for key \"" + key + "\" in configuration.");
        }
    }
}

string Config::gridFilePath() const {
    string value = parseString("GRID_FILE");
    return value;
}

bool Config::restartSolution() const {
    bool value = parseBool("RESTART_SOLUTION", false);
    return value;
}

string Config::restartSolutionFilepath() const {
    string value = parseString("RESTART_SOLUTION_FILEPATH");
    return value;
}

bool Config::isBFMActive() const {
    return parseBool("BFM_ACTIVE", false);
}

bool Config::isBlockageActive() const {
    return parseBool("BLOCKAGE_ACTIVE", false);
}

BFM_Model Config::getBFMModel() const {
    string value = parseString("BFM_MODEL", true);
    BFM_Model model = BFM_Model::NONE;
    if (value == "Hall-Thollet") {
        model = BFM_Model::HALL_THOLLET;
    } else if (value == "Hall") {
        model = BFM_Model::HALL;
    } else if (value == "Lift-Drag") {
        model = BFM_Model::LIFT_DRAG;
    } else if (value == "Frozen-Force") {
        model = BFM_Model::FROZEN_FORCE;
    } else if (value == "Frozen-Gradient") {
        model = BFM_Model::FROZEN_GRADIENT;
    } else if (value == "Gong") {
        model = BFM_Model::GONG;
    } else if (value == "Neri") {
        model = BFM_Model::NERI;
    } else if (value == "Chima") {
        model = BFM_Model::CHIMA;
    } else if (value == "Lamprakis") {
        model = BFM_Model::LAMPRAKIS;
    } else if (value == "None") {
        model = BFM_Model::NONE;
    } else if (value == "Blockage") {
        model = BFM_Model::ONLY_BLOCKAGE;
    } 
    else {
        throw std::runtime_error("Invalid value for key \"BFM_MODEL\" in configuration.");
    }
    return model;
}

KindSolver Config::getKindSolver() const {
    string value = parseString("KIND_SOLVER");
    KindSolver solver = KindSolver::EULER;
    if (value == "Euler") {
        solver = KindSolver::EULER;
    } else if (value == "RANS") {
        solver = KindSolver::RANS;
    } else {
        throw std::runtime_error("Invalid value for key \"KIND_SOLVER\" in configuration.");
    }
    return solver;
}

Topology Config::getTopology() const {
    string value = parseString("TOPOLOGY");
    Topology topology = Topology::TWO_DIMENSIONAL;
    if (value == "2D") {
        topology = Topology::TWO_DIMENSIONAL;
    } else if (value == "3D") {
        topology = Topology::THREE_DIMENSIONAL;
    } else if (value == "axisymmetric_2d") {
        topology = Topology::AXISYMMETRIC_2D;
    } else if (value == "axisymmetric_3d" || value == "axisymmetric") {
        topology = Topology::AXISYMMETRIC_3D;
    } else {
        throw std::runtime_error("Invalid value for key \"TOPOLOGY\" in configuration.");
    }
    return topology;
}

BoundaryType Config::getBoundaryType(BoundaryIndices bcIndex) const{
    BoundaryType boundaryType;
    std::string boundaryString;

    if (bcIndex == BoundaryIndices::I_START){
        boundaryString = parseVector<std::string>("BOUNDARY_TYPE_I")[0];
    } 
    else if (bcIndex == BoundaryIndices::I_END){
        boundaryString = parseVector<std::string>("BOUNDARY_TYPE_I")[1];
    } 
    else if (bcIndex == BoundaryIndices::J_START){
        boundaryString = parseVector<std::string>("BOUNDARY_TYPE_J")[0];
    } 
    else if (bcIndex == BoundaryIndices::J_END){
        boundaryString = parseVector<std::string>("BOUNDARY_TYPE_J")[1];
    }
    else if (bcIndex == BoundaryIndices::K_START){
        boundaryString = parseVector<std::string>("BOUNDARY_TYPE_K")[0];
    }
    else if (bcIndex == BoundaryIndices::K_END){
        boundaryString = parseVector<std::string>("BOUNDARY_TYPE_K")[1];
    } 
    else {
        throw std::runtime_error("Invalid direction Bounday Index");
    }
    
    if (boundaryString == "inlet") {
        boundaryType = BoundaryType::INLET;
    } 
    else if (boundaryString == "inlet_supersonic") {
        boundaryType = BoundaryType::INLET_SUPERSONIC;
    }
    else if (boundaryString == "outlet") {
        boundaryType = BoundaryType::OUTLET;
    } 
    else if (boundaryString == "outlet_supersonic") {
        boundaryType = BoundaryType::OUTLET_SUPERSONIC;
    } 
    else if (boundaryString == "radial_equilibrium") {
        boundaryType = BoundaryType::RADIAL_EQUILIBRIUM;
    }
    else if (boundaryString == "throttle") {
        boundaryType = BoundaryType::THROTTLE;
    } 
    else if (boundaryString == "wall") {
        boundaryType = BoundaryType::WALL;
    } 
    else if (boundaryString == "wedge") {
        boundaryType = BoundaryType::WEDGE;
    }
    else if (boundaryString == "periodic") {
        boundaryType = BoundaryType::PERIODIC;
    } 
    else {
        throw std::runtime_error("Invalid value for key \"BOUNDARY_TYPES\" in configuration.");
    }

    return boundaryType;
}

TimeIntegration Config::getTimeIntegration() const {
    int value = std::stoi(get("TIME_INTEGRATION_TYPE"));
    TimeIntegration integration = TimeIntegration::RUNGE_KUTTA_4;
    if (value == 0) {
        integration = TimeIntegration::RUNGE_KUTTA_4;
    } 
    else if (value == 1) {
        integration = TimeIntegration::RUNGE_KUTTA_3;
    }
    else {
        throw std::runtime_error("Invalid value for key \"TIME_INTEGRATION\" in configuration.");
    }
    return integration;
}

TimeStepMethod Config::getTimeStepMethod() const {
    std::string value = parseString("TIME_STEP_METHOD");
    TimeStepMethod method = TimeStepMethod::LOCAL;

    if (value == "local") {
        method = TimeStepMethod::LOCAL;
    } 
    else if (value == "global") {
        method = TimeStepMethod::GLOBAL;
    }
    else if (value == "fixed") {
        method = TimeStepMethod::FIXED;
    }
    else {
        throw std::runtime_error("Invalid value for key \"TIME_STEP_METHOD\" in configuration.");
    }
    return method;
}

ConvectionScheme Config::getConvectionScheme() const {
    std::string value = parseString("CONVECTION_SCHEME");
    ConvectionScheme scheme = ConvectionScheme::JST;
    if (value == "JST") {
        scheme = ConvectionScheme::JST;
    } 
    else if (value == "ROE" || value == "roe" || value == "Roe") {
        scheme = ConvectionScheme::ROE;
    }
    else {
        throw std::runtime_error("Invalid value for key \"CONVECTION_SCHEME\" in configuration.");
    }
    return scheme;
}


Vector3D Config::getInitDirection() const { 
    std::string raw = get("INIT_DIRECTION"); // getOption retrieves the string value from config
    
    if (raw == "adaptive") {
        // Compute or return a default adaptive direction
        return Vector3D{0.0, 0.0, 0.0}; 
    }

    std::vector<FloatType> direction = parseVector<FloatType>("INIT_DIRECTION");
    if (direction.size() != 3) {
        throw std::runtime_error("INIT_DIRECTION must have 3 components or be 'adaptive'");
    }

    return Vector3D{direction[0], direction[1], direction[2]};
}

std::vector<FloatType> Config::getTimeIntegrationCoeffs() const {
    TimeIntegration integration = getTimeIntegration();
    if (integration == TimeIntegration::RUNGE_KUTTA_4) {
        return {1.0/4.0, 5.0/14.0, 14.0/25.0, 1.0};
    } 
    else if (integration == TimeIntegration::RUNGE_KUTTA_3) {
        return {8.0/17.0, 17.0/20.0, 1.0};
    }
    else {
        throw std::runtime_error("Invalid value for key \"TIME_INTEGRATION\" in configuration.");
    }
}

void Config::printAllConfigValues() const {
    std::cout << std::endl;
    std::cout << "=========================================\n";
    std::cout << "       INFORMATION OF CONFIG FILE        \n";
    std::cout << "=========================================\n";
    std::cout << std::endl;

    for (const auto& pair : _data) {
        std::cout << pair.first << " = " << pair.second << std::endl;
    }
    
    std::cout << "=========================================\n";

    
}


    FloatType Config::computeRampedOutletPressure(const size_t iterCounter, const FloatType outletPressure) const{
        FloatType maxIter = getOutletPressureRampMaxIterations();
        FloatType iterRatio = static_cast<FloatType>(iterCounter) / maxIter;
        if (iterRatio >= 1.0) {
            return outletPressure;
        }
        
        FloatType initialPressure = getInitPressure();
        return initialPressure + (outletPressure - initialPressure) * iterRatio;
    }


FluxLimiter Config::getFluxLimiter() const {
        std::string value = parseString("FLUX_LIMITER", true);
        if (value == "None") {
            return FluxLimiter::NONE;
        } else if (value == "van albada") {
            return FluxLimiter::VAN_ALBADA;
        } else if (value == "van leer") {
            return FluxLimiter::VAN_LEER;
        } else if (value == "min mod") {
            return FluxLimiter::MIN_MOD;
        } else {
            throw std::runtime_error("Invalid value for key \"FLUX_LIMITER\" in configuration.");
        }
    }

    FluidModel Config::getFluidModel() const {
    std::string value = parseString("FLUID_MODEL", true);
    FluidModel model = FluidModel::IDEAL;
    if (value == "None" || value == "Ideal" || value == "ideal" || value == "IDEAL") {
        model = FluidModel::IDEAL;
    } 
    else if (value == "real" || value == "Real" || value == "REAL") {
        model = FluidModel::REAL;
    }
    else {
        throw std::runtime_error("Invalid value for key \"FLUID_MODEL\" in configuration.");
    }
    return model;
}


FloatType Config::getBfmRelaxationFactor() const {
    bool bfmLagActive = isBfmLagActive();
    if (bfmLagActive) {
        return parseFloat("BFM_RELAXATION_FACTOR", 0.001);
    }
    else {
        return 1.0;
    }
}

ReferenceFrame Config::getInletReferenceFrame() const {
    std::string value = parseString("INLET_REFERENCE_FRAME", true);
    ReferenceFrame frame = ReferenceFrame::CARTESIAN;
    if (value == "Cartesian" || value == "cartesian" || value == "CARTESIAN" || value == "None") {
        frame = ReferenceFrame::CARTESIAN;
    } 
    else if (value == "Cylindrical" || value == "cylindrical" || value == "CYLINDRICAL") {
        frame = ReferenceFrame::CYLINDRICAL;
    }
    else {
        throw std::runtime_error("Invalid value for key \"INLET_REFERENCE_FRAME\" in configuration.");
    }
    return frame;
}