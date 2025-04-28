#pragma once
#include <iostream>
#include <vector>
#include <algorithm> // for minmax_element
#include <numeric>   // for accumulate, inner_product
#include <cmath>

using FloatType = double;


// This class represents a 3D vector with basic operations. FloatType is always the individual data type
class Vector3D {
private:
    std::array<FloatType, 3> _data;

public:
    Vector3D() : _data{0.0, 0.0, 0.0} {} // default constructor
    Vector3D(FloatType x, FloatType y, FloatType z) : _data{x, y, z} {} //overloaded

    // Accessors
    FloatType& operator()(size_t index) {
        if (index >= 3) throw std::out_of_range("Index out of bounds");
        return _data[index];
    }

    const FloatType& operator()(size_t index) const {
        if (index >= 3) throw std::out_of_range("Index out of bounds");
        return _data[index];
    }

    FloatType& x() { return _data[0]; }
    FloatType& y() { return _data[1]; }
    FloatType& z() { return _data[2]; }

    const FloatType& x() const { return _data[0]; }
    const FloatType& y() const { return _data[1]; }
    const FloatType& z() const { return _data[2]; }

    // Operator Overloads
    Vector3D operator+(const Vector3D& other) const {
        return Vector3D(x() + other.x(), y() + other.y(), z() + other.z());
    }

    Vector3D operator-(const Vector3D& other) const {
        return Vector3D(x() - other.x(), y() - other.y(), z() - other.z());
    }

    Vector3D operator*(FloatType scalar) const {
        return Vector3D(x() * scalar, y() * scalar, z() * scalar);
    }

    Vector3D operator/(FloatType scalar) const {
        return Vector3D(x() / scalar, y() / scalar, z() / scalar);
    }

    Vector3D& operator+=(const Vector3D& other) {
        x() += other.x(); y() += other.y(); z() += other.z();
        return *this;
    }

    Vector3D& operator-=(const Vector3D& other) {
        x() -= other.x(); y() -= other.y(); z() -= other.z();
        return *this;
    }

    Vector3D operator-() const {
        return Vector3D(- x(), - y(), - z());
    }

    Vector3D& operator*=(FloatType scalar) {
        x() *= scalar; y() *= scalar; z() *= scalar;
        return *this;
    }

    Vector3D& operator/=(FloatType scalar) {
        x() /= scalar; y() /= scalar; z() /= scalar;
        return *this;
    }

    bool operator==(const Vector3D& other) const {
        return x() == other.x() && y() == other.y() && z() == other.z();
    }

    // Magnitude and normalization
    FloatType magnitude() const {
        return std::sqrt(x() * x() + y() * y() + z() * z());
    }

    Vector3D normalized() const {
        FloatType mag = magnitude();
        return (mag == 0) ? Vector3D(0, 0, 0) : *this / mag;
    }

    Vector3D cross(const Vector3D& other) const {
        return Vector3D(
            y() * other.z() - z() * other.y(),
            z() * other.x() - x() * other.z(),
            x() * other.y() - y() * other.x()
        );
    }

    FloatType dot(const Vector3D& other) const {
        FloatType result = x() * other.x() + y() * other.y() + z() * other.z();
        return result;
    }

    // print the content
    friend std::ostream& operator<<(std::ostream& os, const Vector3D& v) {
    return os << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    }
};


enum class BFM_Model {
    NONE = 0,
    HALL = 1,
    HALL_THOLLET = 2,
    LIFT_DRAG = 3
};


enum class BoundaryIndices {
    I_START = 0,
    I_END = 1,
    J_START = 2,
    J_END = 3,
    K_START = 4,
    K_END = 5
};


enum class KindSolver {
    EULER = 0,
    RANS = 1
};


enum class Topology{
    TWO_DIMENSIONAL = 0,
    THREE_DIMENSIONAL = 1,
    AXISYMMETRIC = 2
};


enum class BoundaryType {
    INLET = 0,
    OUTLET = 1,
    OUTLET_RADIAL_EQUILIBRIUM = 2,
    WALL = 3,
    EMPTY = 4,
    WEDGE = 5
};


enum class TimeIntegration {
    RUNGE_KUTTA_4 = 0,
    RUNGE_KUTTA_3 = 1
};


enum class TimeStepMethod {
    LOCAL = 0,
    GLOBAL = 1
};


enum class ConvectionScheme {
    JST = 0
};


// 2D matrix that can include floats or vectors inside
template<typename T>
class Matrix2D {
public:
    // Default constructor
    Matrix2D() : _ni(0), _nj(0) {}

    // Constructor with dimensions
    Matrix2D(size_t ni, size_t nj) : _ni(ni), _nj(nj), _data(ni * nj) {}

    // Resize the matrix to new dimensions
    void resize(size_t ni, size_t nj) {
        _ni = ni;
        _nj = nj;
        _data.resize(ni * nj);
    }

    // Overloaded () operator for access with bounds checking
    T& operator()(size_t i, size_t j) {
        check_bounds(i, j);
        return _data[i * _nj + j];
    }

    const T& operator()(size_t i, size_t j) const {
        check_bounds(i, j);
        return _data[i * _nj + j];
    }

    size_t sizeI() const { return _ni; }
    size_t sizeJ() const { return _nj; }

private:
    size_t _ni, _nj;
    std::vector<T> _data;

    void check_bounds(size_t i, size_t j) const {
        if (i >= _ni || j >= _nj) {
            throw std::out_of_range("Matrix2D::operator() index out of range");
        }
    }
};


// 3D matrix that can include floats or vectors inside
template <typename T>
class Matrix3D {
public:
    // Default constructor
    Matrix3D() : _ni(0), _nj(0), _nk(0) {}

    // Constructor with dimensions and optional initial data
    Matrix3D(size_t ni, size_t nj, size_t nk, const std::vector<T>& initial_data = {})
        : _ni(ni), _nj(nj), _nk(nk) {
        if (!initial_data.empty()) {
            if (initial_data.size() != ni * nj * nk) {
                throw std::invalid_argument("Initial data size does not match matrix dimensions");
            }
            _data = initial_data;
        } else {
            _data.resize(ni * nj * nk);
        }
    }

    // Resize the matrix
    void resize(size_t ni, size_t nj, size_t nk) {
        _ni = ni;
        _nj = nj;
        _nk = nk;
        _data.resize(ni * nj * nk);
    }

    // Access and modify using operator() with bounds checking
    T& operator()(size_t i, size_t j, size_t k) {
        check_bounds(i, j, k);
        return _data[i * _nj * _nk + j * _nk + k];
    }

    const T& operator()(size_t i, size_t j, size_t k) const {
        check_bounds(i, j, k);
        return _data[i * _nj * _nk + j * _nk + k];
    }

    // Size getters
    size_t sizeI() const { return _ni; }
    size_t sizeJ() const { return _nj; }
    size_t sizeK() const { return _nk; }

    // Compute L2 norm of all components of the matrix
    T norm() const {
        T sum = 0;
        for (size_t i = 0; i < _ni; ++i) {
            for (size_t j = 0; j < _nj; ++j) {
                for (size_t k = 0; k < _nk; ++k) {
                    sum += std::pow(operator()(i, j, k), 2);  // Square each element
                }
            }
        }
        return std::sqrt(sum);  // Return the square root of the sum
    }

    // In-place addition
    Matrix3D& operator+=(const Matrix3D& other) {
        if (_ni != other._ni || _nj != other._nj || _nk != other._nk) {
            throw std::invalid_argument("Matrix3D::operator+= dimension mismatch");
        }
        for (size_t idx = 0; idx < _data.size(); ++idx) {
            _data[idx] += other._data[idx];
        }
        return *this;
    }

    // In-place multiplication
    Matrix3D& operator*=(const Matrix3D& other) {
        if (_ni != other._ni || _nj != other._nj || _nk != other._nk) {
            throw std::invalid_argument("Matrix3D::operator+= dimension mismatch");
        }
        for (size_t idx = 0; idx < _data.size(); ++idx) {
            _data[idx] *= other._data[idx];
        }
        return *this;
    }

    Matrix3D& operator*=(T scalar) {
        for (size_t idx = 0; idx < _data.size(); ++idx) {
            _data[idx] *= scalar;
        }
        return *this;
    }

    // In-place subtraction
    Matrix3D& operator-=(const Matrix3D& other) {
        if (_ni != other._ni || _nj != other._nj || _nk != other._nk) {
            throw std::invalid_argument("Matrix3D::operator-= dimension mismatch");
        }
        for (size_t idx = 0; idx < _data.size(); ++idx) {
            _data[idx] -= other._data[idx];
        }
        return *this;
    }

    std::vector<T> getData() const {return _data;}

private:
    size_t _ni, _nj, _nk;
    std::vector<T> _data;

    void check_bounds(size_t i, size_t j, size_t k) const {
        if (i >= _ni || j >= _nj || k >= _nk) {
            throw std::out_of_range("Matrix3D::operator() index out of range");
        }
    }
};


// Wrapper of std::array<FloatType, 5> that is used for primitive, conservative variables. Overloaded methods
class StateVector {
    private:
        static constexpr std::size_t Size = 5; // for the moment always 5, later it could grow
        std::array<FloatType, Size> _data;
    
    public:
        // Default constructor
        StateVector() = default;
    
        // Constructor from std::array
        explicit StateVector(const std::array<FloatType, Size>& arr) : _data(arr) {}
    
        // Access operators
        FloatType& operator[](std::size_t idx) { return _data[idx]; }

        const FloatType& operator[](std::size_t idx) const { return _data[idx]; }
    
        // Element-wise addition
        StateVector operator+(const StateVector& other) const {
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i)
                result[i] = _data[i] + other[i];
            return result;
        }
    
        // Element-wise subtraction
        StateVector operator-(const StateVector& other) const {
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i)
                result[i] = _data[i] - other[i];
            return result;
        }
    
        // Element-wise multiplication
        StateVector operator*(const StateVector& other) const {
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i)
                result[i] = _data[i] * other[i];
            return result;
        }
    
        // Element-wise division
        StateVector operator/(const StateVector& other) const {
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i) {
                if (other[i] == 0.0)
                    throw std::runtime_error("Division by zero in StateVector");
                result[i] = _data[i] / other[i];
            }
            return result;
        }
    
        // Scalar operations
        StateVector operator*(FloatType scalar) const {
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i)
                result[i] = _data[i] * scalar;
            return result;
        }
    
        StateVector operator/(FloatType scalar) const {
            if (scalar == 0.0)
                throw std::runtime_error("Division by zero in scalar operation");
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i)
                result[i] = _data[i] / scalar;
            return result;
        }
    
        StateVector operator+(FloatType scalar) const {
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i)
                result[i] = _data[i] + scalar;
            return result;
        }
    
        StateVector operator-(FloatType scalar) const {
            StateVector result;
            for (std::size_t i = 0; i < Size; ++i)
                result[i] = _data[i] - scalar;
            return result;
        }
    
        // Assignment operators
        StateVector& operator+=(const StateVector& other) {
            for (std::size_t i = 0; i < Size; ++i)
                _data[i] += other[i];
            return *this;
        }
    
        StateVector& operator-=(const StateVector& other) {
            for (std::size_t i = 0; i < Size; ++i)
                _data[i] -= other[i];
            return *this;
        }
    
        // Convenience conversion to std::array
        operator std::array<FloatType, Size>() const {
            return _data;
        }
    
        // dot product
        FloatType dot(const StateVector& other) const {
            FloatType result = 0.0;
            for (std::size_t i = 0; i < Size; ++i)
                result += _data[i] * other[i];
            return result;
        }
    
        // Optional: L2 norm
        FloatType norm() const {
            return std::sqrt(this->dot(*this));
        }
    };



struct FlowSolution {
    Matrix3D<FloatType> _rho;
    Matrix3D<FloatType> _rhoU;
    Matrix3D<FloatType> _rhoV;
    Matrix3D<FloatType> _rhoW;
    Matrix3D<FloatType> _rhoE;

    // Constructor
    FlowSolution(size_t ni, size_t nj, size_t nk) {
        resize(ni, nj, nk);
    }

    FlowSolution() = default;

    FloatType norm(size_t iVar) const {
        switch (iVar)
        {
        case 0:
            return _rho.norm();
        case 1:
            return _rhoU.norm();
        case 2:
            return _rhoV.norm();
        case 3:
            return _rhoW.norm();
        case 4:
            return _rhoE.norm();        
        default:
            throw std::out_of_range("Index out of range for norm calculation");
        }
    }

    void add(size_t i, size_t j, size_t k, const StateVector& delta) {
        _rho(i,j,k)   += delta[0];
        _rhoU(i,j,k)  += delta[1];
        _rhoV(i,j,k)  += delta[2];
        _rhoW(i,j,k)  += delta[3];
        _rhoE(i,j,k)  += delta[4];
    }

    void subtract(size_t i, size_t j, size_t k, const StateVector& delta) {
        _rho(i,j,k)   -= delta[0];
        _rhoU(i,j,k)  -= delta[1];
        _rhoV(i,j,k)  -= delta[2];
        _rhoW(i,j,k)  -= delta[3];
        _rhoE(i,j,k)  -= delta[4];
    }

    StateVector at(size_t i, size_t j, size_t k) const {
        return StateVector({_rho(i,j,k), _rhoU(i,j,k), _rhoV(i,j,k), _rhoW(i,j,k), _rhoE(i,j,k)});
    }

    void set(size_t i, size_t j, size_t k, const StateVector& vals) {
        _rho(i,j,k) = vals[0];
        _rhoU(i,j,k) = vals[1];
        _rhoV(i,j,k) = vals[2];
        _rhoW(i,j,k) = vals[3];
        _rhoE(i,j,k) = vals[4];
    }

    void set(size_t i, size_t j, size_t k, size_t eq, FloatType val) {
        switch (eq) {
            case 0:
                _rho(i, j, k) = val;
                break;
            case 1:
                _rhoU(i, j, k) = val;
                break;
            case 2:
                _rhoV(i, j, k) = val;
                break;
            case 3:
                _rhoW(i, j, k) = val;
                break;
            case 4:
                _rhoE(i, j, k) = val;
                break;
            default:
                throw std::out_of_range("Invalid equation index: " + std::to_string(eq));
        }
    }

    void resize(size_t i, size_t j, size_t k){
        _rho.resize(i,j,k);
        _rhoU.resize(i,j,k);
        _rhoV.resize(i,j,k);
        _rhoW.resize(i,j,k);
        _rhoE.resize(i,j,k);
    }

    // Addition: FlowSolution + FlowSolution
    FlowSolution operator+(const FlowSolution& other) const {
        FlowSolution result = *this;
        result += other;
        return result;
    }

    // In-place addition: FlowSolution += FlowSolution
    FlowSolution& operator+=(const FlowSolution& other) {
        _rho  += other._rho;
        _rhoU += other._rhoU;
        _rhoV += other._rhoV;
        _rhoW += other._rhoW;
        _rhoE += other._rhoE;
        return *this;
    }

    FlowSolution& operator*(const Matrix3D<FloatType>& other) {
        _rho  *= other;
        _rhoU += other;
        _rhoV += other;
        _rhoW += other;
        _rhoE += other;
        return *this;
    }

    // Subtraction: FlowSolution - FlowSolution
    FlowSolution operator-(const FlowSolution& other) const {
        FlowSolution result = *this;
        result -= other;
        return result;
    }

    // In-place subtraction: FlowSolution -= FlowSolution
    FlowSolution& operator-=(const FlowSolution& other) {
        _rho  -= other._rho;
        _rhoU -= other._rhoU;
        _rhoV -= other._rhoV;
        _rhoW -= other._rhoW;
        _rhoE -= other._rhoE;
        return *this;
    }

    FlowSolution& operator*=(FloatType scalar) {
        _rho  *= scalar;
        _rhoU *= scalar;
        _rhoV *= scalar;
        _rhoW *= scalar;
        _rhoE *= scalar;
        return *this;
    }

    // Multiplication: FlowSolution * scalar
    FlowSolution operator*(FloatType scalar) const {
        FlowSolution result = *this;
        result *= scalar;   // GOOD: calls the in-place multiplication, no recursion
        return result;
    }
};

enum class FluxDirection {
    I=0,
    J=1,
    K=2,
};


struct Statistics {
    FloatType mean = 0;
    FloatType min = 0;
    FloatType max = 0;
    FloatType stdDev = 0;
    size_t size = 0;

    Statistics() = default; // allow default construction

    // Method that computes statistics
    Statistics(const Matrix3D<FloatType>& matrix) {
        if (matrix.sizeI() * matrix.sizeJ() * matrix.sizeK() == 0) {
            mean = min = max = stdDev = 0;
            size = 0;
            return;
        }

        std::vector<FloatType> data = matrix.getData();
        auto [minIt, maxIt] = std::minmax_element(data.begin(), data.end());
        min = *minIt;
        max = *maxIt;

        FloatType sum = std::accumulate(data.begin(), data.end(), static_cast<FloatType>(0));
        mean = sum / static_cast<FloatType>(data.size());

        FloatType sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), static_cast<FloatType>(0));
        stdDev = std::sqrt(sq_sum / static_cast<FloatType>(data.size()) - mean * mean);

        size = data.size();
    }
};