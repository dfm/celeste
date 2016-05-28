#ifndef _CELESTE_CENTRAL_BODY_H_
#define _CELESTE_CENTRAL_BODY_H_

#include <cstddef>

namespace celeste {

class CentralBody {
friend class Body;
public:
    CentralBody(double log_mass = 0.0, double log_radius = 0.0)
        : log_mass_(log_mass), log_radius_(log_radius) {};

    std::size_t size () const { return 2; };

    void get_parameter_vector (double* params) const {
        params[0] = log_mass_;
        params[1] = log_radius_;
    };

    void set_parameter_vector (const double* params) {
        log_mass_   = params[0];
        log_radius_ = params[1];
    };

private:
    double log_mass_, log_radius_;

};  // class CentralBody

};  // namespace celeste

#endif  // _CELESTE_CENTRAL_BODY_H_
