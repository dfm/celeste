#ifndef _CELESTE_CENTRAL_BODY_H_
#define _CELESTE_CENTRAL_BODY_H_

#include "celeste/parameter.h"

namespace celeste {

class CentralBody {
public:
    // The parameters are public to enable freezing and thawing.
    Parameter<double> log_mass, log_radius;

    CentralBody(
        double log_mass0 = 0.0,      // Solar masses
        double log_radius0 = 0.0    // Solar radii
    )
        : log_mass(log_mass0),
          log_radius(log_radius0) {};

    // Get the size of the unfrozen parameter vector.
    size_t size () const {
        size_t count = 0;
        COUNT_PARAM(log_mass)
        COUNT_PARAM(log_radius)
        return count;
    };

    // Parameter vector interface.
    void get_parameter_vector (double* params) const {
        size_t count = 0;
        GET_PARAM(log_mass)
        GET_PARAM(log_radius)
    };
    void set_parameter_vector (const double* params) {
        size_t count = 0;
        SET_PARAM(log_mass)
        SET_PARAM(log_radius)
    };

};  // class CentralBody

};  // namespace celeste

#endif  // _CELESTE_CENTRAL_BODY_H_
