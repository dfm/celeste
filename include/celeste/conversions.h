#ifndef _CELESTE_CENTRAL_BODY_H_
#define _CELESTE_CENTRAL_BODY_H_

#include <cmath>

namespace celeste {

double get_period (const CentralBody& central, const Body& body)
{
    double a = body.get_semimajor(),
           mu = exp(central.log_mass_)+exp(body.log_mass_);
    return 2.0*M_PI*sqrt(a*a*a/G_GRAV/mu);
}

double set_period (const CentralBody& central, Body& body, double period)
{
    double mu = exp(central.log_mass_)+exp(body.log_mass_),
           a = pow(G_GRAV*period*period*mu/(4*M_PI*M_PI), 1./3.);
    body.set_semimajor_terms(a, body.ix, body.iy);
}

};  // namespace celeste

#endif  // _CELESTE_CENTRAL_BODY_H_
