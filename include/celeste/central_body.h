#ifndef _CELESTE_CENTRAL_BODY_H_
#define _CELESTE_CENTRAL_BODY_H_

#include <vector>

#include "celeste/parameter.h"

namespace celeste {

class CentralBody : public ParameterizedModel {
friend class System;
public:
    // The parameters are public to enable freezing and thawing.
    Parameter log_mass, log_radius;

    CentralBody(
        double log_mass0 = 0.0,      // Solar masses
        double log_radius0 = 0.0    // Solar radii
    )
        : log_mass(log_mass0),
          log_radius(log_radius0) {};

    size_t full_size () const { return 2; };
    GET_PARAM_FUNCTION(
        GET_PARAM(log_mass)
        GET_PARAM(log_radius)
    )
    SET_PARAM_FUNCTION(
        SET_PARAM(log_mass)
        SET_PARAM(log_radius)
    )

};  // class CentralBody

};  // namespace celeste

#endif  // _CELESTE_CENTRAL_BODY_H_
