#pragma once

using FloatType = double;

struct Point3D {
    FloatType x, y, z;
};

enum class BFM_Model {
    NONE = 0,
    HALL = 1,
    HALL_THOLLET = 2,
    LIFT_DRAG = 3
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
class Matrix1D {
public:
    Matrix1D(size_t n)
        : _n(n), _data(n) {}

    // Access and modify with bounds checking
    T& at(size_t i) {
        check_bounds(i);
        return _data[i];
    }

    const T& at(size_t i) const {
        check_bounds(i);
        return _data[i];
    }

    size_t size() const { return _n; }

private:
    size_t _n;
    std::vector<T> _data;

    void check_bounds(size_t i) const {
        if (i >= _n) {
            throw std::out_of_range("Matrix1D::at() index out of range");
        }
    }
};

template <typename T>
class Matrix2D {
public:
    // Default constructor (initially no size)
    Matrix2D() : _ni(0), _nj(0) {}

    // Constructor with known dimensions
    Matrix2D(size_t ni, size_t nj, const std::vector<T>& initial_data = {})
        : _ni(ni), _nj(nj) {
        if (!initial_data.empty()) {
            _data = initial_data;
        } else {
            _data.resize(ni * nj);
        }
    }

    // Resize the matrix to new dimensions
    void resize(size_t ni, size_t nj) {
        _ni = ni;
        _nj = nj;
        _data.resize(ni * nj);
    }

    // Access and modify with bounds checking
    T& at(size_t i, size_t j) {
        check_bounds(i, j);
        return _data[i * _nj + j];
    }

    const T& at(size_t i, size_t j) const {
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
            throw std::out_of_range("Matrix2D::at() index out of range");
        }
    }
};


template <typename T>
class Matrix3D {
public:
    // Default constructor (initially no size)
    Matrix3D() : _ni(0), _nj(0), _nk(0) {}

    // Constructor with known dimensions
    Matrix3D(size_t ni, size_t nj, size_t nk, const std::vector<T>& initial_data = {})
        : _ni(ni), _nj(nj), _nk(nk) {
        if (!initial_data.empty()) {
            _data = initial_data;
        } else {
            _data.resize(ni * nj * nk);
        }
    }

    // Resize the matrix to new dimensions
    void resize(size_t ni, size_t nj, size_t nk) {
        _ni = ni;
        _nj = nj;
        _nk = nk;
        _data.resize(ni * nj * nk);
    }

    // Access and modify with bounds checking
    T& at(size_t i, size_t j, size_t k) {
        check_bounds(i, j, k);
        return _data[i * _nj * _nk + j * _nk + k];
    }

    const T& at(size_t i, size_t j, size_t k) const {
        check_bounds(i, j, k);
        return _data[i * _nj * _nk + j * _nk + k];
    }

    size_t sizeI() const { return _ni; }
    size_t sizeJ() const { return _nj; }
    size_t sizeK() const { return _nk; }

private:
    size_t _ni, _nj, _nk;
    std::vector<T> _data;

    void check_bounds(size_t i, size_t j, size_t k) const {
        if (i >= _ni || j >= _nj || k >= _nk) {
            throw std::out_of_range("Matrix3D::at() index out of range");
        }
    }
};



template <typename T>
class Matrix4D {
public:
    // Default constructor (initially no size)
    Matrix4D() : _ni(0), _nj(0), _nk(0), _nl(0) {}

    // Constructor with known dimensions
    Matrix4D(size_t ni, size_t nj, size_t nk, size_t nl, const std::vector<T>& initial_data = {})
        : _ni(ni), _nj(nj), _nk(nk), _nl(nl) {
        if (!initial_data.empty()) {
            _data = initial_data;
        } else {
            _data.resize(ni * nj * nk * nl);
        }
    }

    // Resize the matrix to new dimensions
    void resize(size_t ni, size_t nj, size_t nk, size_t nl) {
        _ni = ni;
        _nj = nj;
        _nk = nk;
        _nl = nl;
        _data.resize(ni * nj * nk * nl);
    }

    // Access and modify with bounds checking
    T& at(size_t i, size_t j, size_t k, size_t l) {
        check_bounds(i, j, k, l);
        return _data[i * _nj * _nk * _nl + j * _nk * _nl + k * _nl + l];
    }

    const T& at(size_t i, size_t j, size_t k, size_t l) const {
        check_bounds(i, j, k, l);
        return _data[i * _nj * _nk * _nl + j * _nk * _nl + k * _nl + l];
    }

    size_t sizeI() const { return _ni; }
    size_t sizeJ() const { return _nj; }
    size_t sizeK() const { return _nk; }
    size_t sizeL() const { return _nl; }

private:
    size_t _ni, _nj, _nk, _nl;
    std::vector<T> _data;

    void check_bounds(size_t i, size_t j, size_t k, size_t l) const {
        if (i >= _ni || j >= _nj || k >= _nk || l >= _nl) {
            throw std::out_of_range("Matrix4D::at() index out of range");
        }
    }
};

