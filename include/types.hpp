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