#ifndef POTENTIALS_AND_FORCES_HPP
#define POTENTIALS_AND_FORCES_HPP

#include "Vector2D.hpp"
#include <cmath>

namespace potentials {

constexpr double R_M = 1.737051787;

// phi*(x) — dimensionless modified LJ potential, where x = r/sigma.
inline double lj_potential(double x) {
    if (x <= 0.85) {
        return -39170.516264*x*x*x + 102795.709007*x*x - 90115.812549*x + 26401.649473;
    }
    if (x <= 1.24445506) {
        double const x6  = x*x*x * x*x*x;
        double const x12 = x6*x6;
        return 4.0 * (1.0/x12 - 1.0/x6);
    }
    if (x <= R_M) {
        double u = R_M - x;
        return 3.2920028*u*u*u - 4.86489008*u*u;
    }
    return 0.0;
}

// -d(phi*)/dx — positive means repulsive.
inline double lj_force_magnitude_signed(double x) {
    if (x <= 0.85) {
        return -(3.0*(-39170.516264)*x*x + 2.0*102795.709007*x - 90115.812549);
    }
    if (x <= 1.24445506) {
        double x7  = x*x*x * x*x*x * x;
        double x13 = x7 * x*x*x * x*x*x;
        return 48.0/x13 - 24.0/x7;
    }
    if (x <= R_M) {
        double u = R_M - x;
        return 3.0*3.2920028*u*u + 2.0*(-4.86489008)*u;
    }
    return 0.0;
}

// Force on particle i due to particle j, given r_ij = r_i - r_j
inline Vector2D pairwise_force(const Vector2D& r_ij) {
    double x = r_ij.norm();
    if (x < 1e-12 || x > R_M) return Vector2D(0.0, 0.0);
    double f = lj_force_magnitude_signed(x);
    return r_ij * (f / x);
}

} // namespace potentials

#endif // POTENTIALS_AND_FORCES_HPP
