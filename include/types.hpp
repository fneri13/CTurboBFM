#pragma once
#include <iostream>

using FloatType = double;


// This class represents a 3D vector with basic operations
class Vector3D {
private:
    std::array<FloatType, 3> _data;

public:
    // Constructors
    Vector3D() : _data{0.0, 0.0, 0.0} {}
    Vector3D(FloatType x, FloatType y, FloatType z) : _data{x, y, z} {}

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

    // Add this inside the Vector3D class (outside public/private blocks is fine)
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

private:
    size_t _ni, _nj, _nk;
    std::vector<T> _data;

    void check_bounds(size_t i, size_t j, size_t k) const {
        if (i >= _ni || j >= _nj || k >= _nk) {
            throw std::out_of_range("Matrix3D::operator() index out of range");
        }
    }
};


