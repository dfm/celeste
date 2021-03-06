#ifndef _CELESTE_BODY_H_
#define _CELESTE_BODY_H_

#include <cmath>
#include <cfloat>

#include "celeste/parameter.h"
#include "celeste/central_body.h"

namespace celeste {

#define G_GRAV 2945.462538537765

#define INVALID_IMPACT_PARAMETER 1

class Body : public ParameterizedModel {
friend class System;
public:
    // The parameters are public to enable freezing and thawing.
    Parameter log_mass, log_radius, ref_time,
              sqrt_a_sin_ix_cos_iy, sqrt_a_cos_ix_cos_iy, sqrt_a_sin_iy,
              sqrt_e_sin_omega, sqrt_e_cos_omega;

    Body(
        double log_mass = 0.0,          // Solar masses
        double log_radius = 0.0,        // Solar radii
        double ref_time = 0.0,          // days

        double log_semimajor = 0.0,     // Solar radii
        double ix = 0.0,                // radians
        double iy = 0.0,                // radians

        double eccen = 0.0,             // dimensionless
        double omega = 0.0              // radians
    )
        : log_mass(log_mass),
          log_radius(log_radius),
          ref_time(ref_time)
    {
          set_semimajor_terms(log_semimajor, ix, iy);
          set_eccen_terms(eccen, omega);
    };

    size_t full_size () const { return 8; };
    GET_PARAM_FUNCTION(
        GET_PARAM(log_mass)
        GET_PARAM(log_radius)
        GET_PARAM(ref_time)
        GET_PARAM(sqrt_a_sin_ix_cos_iy)
        GET_PARAM(sqrt_a_cos_ix_cos_iy)
        GET_PARAM(sqrt_a_sin_iy)
        GET_PARAM(sqrt_e_sin_omega)
        GET_PARAM(sqrt_e_cos_omega)
    )
    SET_PARAM_FUNCTION(
        SET_PARAM(log_mass)
        SET_PARAM(log_radius)
        SET_PARAM(ref_time)
        SET_PARAM(sqrt_a_sin_ix_cos_iy)
        SET_PARAM(sqrt_a_cos_ix_cos_iy)
        SET_PARAM(sqrt_a_sin_iy)
        SET_PARAM(sqrt_e_sin_omega)
        SET_PARAM(sqrt_e_cos_omega)
    )

    // Reparameterizations.
    double get_semimajor () const {
        return sqrt_a_sin_ix_cos_iy * sqrt_a_sin_ix_cos_iy
             + sqrt_a_cos_ix_cos_iy * sqrt_a_cos_ix_cos_iy
             + sqrt_a_sin_iy        * sqrt_a_sin_iy;
    };

    double get_eccen () const {
        return sqrt_e_sin_omega * sqrt_e_sin_omega
             + sqrt_e_cos_omega * sqrt_e_cos_omega;
    };

    double get_omega () const {
        return atan2(sqrt_e_sin_omega, sqrt_e_cos_omega);
    };

    double get_ix () const {
        return atan2(sqrt_a_sin_ix_cos_iy, sqrt_a_cos_ix_cos_iy);
    };
    double get_iy () const {
        return atan2(sqrt(sqrt_a_sin_ix_cos_iy * sqrt_a_sin_ix_cos_iy
                        + sqrt_a_cos_ix_cos_iy * sqrt_a_cos_ix_cos_iy),
                     sqrt_a_sin_iy);
    };

    double get_inclination () const {
        return 0.5*M_PI + get_ix();
    };
    void set_inclination (double incl) {
        set_semimajor_terms(log(get_semimajor()), incl - 0.5*M_PI, get_iy());
    };

    double get_period (const CentralBody& central) const {
        double a = get_semimajor(),
               mu = exp(central.log_mass)+exp(log_mass);
        return 2.0*M_PI*sqrt(a*a*a/G_GRAV/mu);
    };
    void set_period (const CentralBody& central, double period) {
        double mu = exp(central.log_mass)+exp(log_mass),
               log_a = log(G_GRAV*period*period*mu/(4*M_PI*M_PI)) / 3.0;
        set_semimajor_terms(log_a);
    };

    double get_impact (const CentralBody& central) const {
        double rstar = exp(central.log_radius),
               incl = get_inclination(),
               e = get_eccen(),
               factor = 1.0;
        if (e > DBL_EPSILON) factor = (1.0-e*e)/(1.0+e*sin(get_omega()));
        return get_semimajor() * cos(incl) / rstar * factor;
    };
    void set_impact (const CentralBody& central, double impact) {
        double e = get_eccen(),
               arg = impact * exp(central.log_radius) / get_semimajor();
        if (e > DBL_EPSILON) arg *= (1.0+e*sin(get_omega()))/(1.0-e*e);
        if (arg > 1.0 || arg < -1.0) throw INVALID_IMPACT_PARAMETER;
        set_inclination(acos(arg));
    };

private:
    void set_semimajor_terms (double log_a, double ix, double iy) {
        double sqrt_a = exp(0.5 * log_a),
               cos_iy = cos(iy);
        sqrt_a_sin_ix_cos_iy = sqrt_a * sin(ix) * cos_iy;
        sqrt_a_cos_ix_cos_iy = sqrt_a * cos(ix) * cos_iy;
        sqrt_a_sin_iy        = sqrt_a * sin(iy);
    };

    void set_semimajor_terms (double log_a) {
        double f = exp(0.5 * log_a) / sqrt(get_semimajor());
        sqrt_a_sin_ix_cos_iy = f * sqrt_a_sin_ix_cos_iy;
        sqrt_a_cos_ix_cos_iy = f * sqrt_a_cos_ix_cos_iy;
        sqrt_a_sin_iy        = f * sqrt_a_sin_iy;
    };

    void set_eccen_terms (double eccen, double omega) {
        double sqrt_e = sqrt(eccen);
        sqrt_e_sin_omega = sqrt_e * sin(omega);
        sqrt_e_cos_omega = sqrt_e * cos(omega);
    }

};  // class Body

};  // namespace celeste

#endif  // _CELESTE_BODY_H_
