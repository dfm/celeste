#ifndef _CELESTE_BODY_H_
#define _CELESTE_BODY_H_

#include <cmath>
#include <cstddef>
#include "celeste/central_body.h"

namespace celeste {

#define G_GRAV 2945.462538537765

class Body {
public:
    Body(
        double log_mass = 0.0,
        double log_radius = 0.0,
        double ref_time = 0.0,

        double log_semimajor = 0.0,
        double ix = 0.0,
        double iy = 0.0,

        double eccen = 0.0,
        double omega = 0.0,

        const CentralBody* central = NULL
    )
        : log_mass_(log_mass),
          log_radius_(log_radius),
          ref_time_(ref_time),
          central_(central)
    {
          set_semimajor_terms(log_semimajor, ix, iy);
          set_eccen_terms(eccen, omega);
    };

    std::size_t size () const { return 8; };

    void set_central (const CentralBody* central) { central_ = central; };

    double get_semimajor () const {
        return sqrt_a_sin_ix_cos_iy_ * sqrt_a_sin_ix_cos_iy_
             + sqrt_a_cos_ix_cos_iy_ * sqrt_a_cos_ix_cos_iy_
             + sqrt_a_sin_iy_        * sqrt_a_sin_iy_;
    };

    double get_period () const {
        double mstar = exp(central_.log_mass),
               a = get_semimajor();
        return 2 * np.pi * np.sqrt(a * a * a / _G / (mstar + self.mass))
    }

    void get_parameter_vector (double* params) const {
        params[0] = log_mass_;
        params[1] = log_radius_;
        params[2] = ref_time_;
        params[3] = sqrt_a_sin_ix_cos_iy_;
        params[4] = sqrt_a_cos_ix_cos_iy_;
        params[5] = sqrt_a_sin_iy_;
        params[6] = sqrt_e_sin_omega_;
        params[7] = sqrt_e_cos_omega_;
    };

    void set_parameter_vector (const double* params) {
        log_mass_             = params[0];
        log_radius_           = params[1];
        ref_time_             = params[2];
        sqrt_a_sin_ix_cos_iy_ = params[3];
        sqrt_a_cos_ix_cos_iy_ = params[4];
        sqrt_a_sin_iy_        = params[5];
        sqrt_e_sin_omega_     = params[6];
        sqrt_e_cos_omega_     = params[7];
    };

private:

    void set_semimajor_terms (double log_a, double ix, double iy) {
        double sqrt_a = exp(0.5 * log_a),
               cos_iy = cos(iy);
        sqrt_a_sin_ix_cos_iy_ = sqrt_a * sin(ix) * cos_iy;
        sqrt_a_cos_ix_cos_iy_ = sqrt_a * cos(ix) * cos_iy;
        sqrt_a_sin_iy_        = sqrt_a * sin(iy);
    }

    void set_eccen_terms (double eccen, double omega) {
        double sqrt_e = sqrt(eccen);
        sqrt_e_sin_omega_ = sqrt_e * sin(omega);
        sqrt_e_cos_omega_ = sqrt_e * cos(omega);
    }

    const CentralBody* central_;
    double log_mass_, log_radius_, ref_time_,
           sqrt_a_sin_ix_cos_iy_, sqrt_a_cos_ix_cos_iy_, sqrt_a_sin_iy_,
           sqrt_e_sin_omega_, sqrt_e_cos_omega_;

};  // class Body

};  // namespace celeste

#endif  // _CELESTE_BODY_H_
