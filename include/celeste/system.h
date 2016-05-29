#ifndef _CELESTE_SYSTEM_H_
#define _CELESTE_SYSTEM_H_

#include <vector>

#include <Eigen/Core>
#include <unsupported/Eigen/AutoDiff>

#include "celeste/body.h"
#include "celeste/central_body.h"

namespace celeste {

typedef Eigen::AutoDiffScalar<Eigen::VectorXd> DiffType;
typedef Eigen::Matrix<DiffType, Eigen::Dynamic, 1> GradVector;

class System {
public:
    System () {};
    System (CentralBody central) : central_(central) {};
    void add_body (Body body) {
        bodies_.push_back(body);
    };

    size_t full_size () const {
        size_t count = central_.full_size();
        for (size_t i = 0; i < bodies_.size(); ++i) count += bodies_[i].full_size();
        return count;
    };
    size_t size () const {
        size_t count = central_.size();
        for (size_t i = 0; i < bodies_.size(); ++i) count += bodies_[i].size();
        return count;
    };

    CentralBody& central () { return central_; };
    Body& body (size_t index = 0) { return bodies_[index]; };

    Eigen::VectorXd get_parameter_vector () const {
        Eigen::VectorXd params(size());
        size_t count = central_.get_parameter_vector(&(params[0]));
        for (size_t i = 0; i < bodies_.size(); ++i)
            count += bodies_[i].get_parameter_vector(&(params[count]));
        return params;
    };

    void set_parameter_vector (const Eigen::VectorXd& params) {
        size_t count = central_.set_parameter_vector(&(params[0]));
        for (size_t i = 0; i < bodies_.size(); ++i)
            count += bodies_[i].set_parameter_vector(&(params[count]));
    };

    Eigen::VectorXd get_full_vector () const {
        Eigen::VectorXd params(full_size());
        size_t count = central_.get_full_vector(&(params[0]));
        for (size_t i = 0; i < bodies_.size(); ++i)
            count += bodies_[i].get_full_vector(&(params[count]));
        return params;
    };

    GradVector get_full_gradient_vector () const {
        size_t nder = size(), derind = 0, count = 0;
        GradVector params(full_size());
        for (size_t k = 0; k < central_.full_size(); ++k) {
            const Parameter& par = central_.parameter(k);
            if (par.is_frozen()) params[count++] = DiffType(par.get_value());
            else params[count++] = DiffType(par.get_value(), nder, derind++);
        }
        for (size_t i = 0; i < bodies_.size(); ++i) {
            const Body& body = bodies_[i];
            for (size_t k = 0; k < body.full_size(); ++k) {
                const Parameter& par = body.parameter(k);
                if (par.is_frozen()) params[count++] = DiffType(par.get_value());
                else params[count++] = DiffType(par.get_value(), nder, derind++);
            }
        }
        return params;
    };

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, 1> kepler_reparameterize (const Eigen::Matrix<T, Eigen::Dynamic, 1> params) const {
        Eigen::Matrix<T, Eigen::Dynamic, 1> result(2 + 11 * bodies_.size());

        result(0) = exp(params(0));  // Central mass
        result(1) = exp(params(1));  // Central radius
        T m_central = result(1);

        size_t i, j, n;
        for (i = 0, j = 2, n = 2; i < bodies_.size(); ++i, j += 8, n += 11) {
            // Access the parameters.
            T mass               = exp(params[j]),
              t0                 = params[j + 2],
              sqrt_a_cosiy_cosix = params[j + 3],
              sqrt_a_cosiy_sinix = params[j + 4],
              sqrt_a_siniy       = params[j + 5],
              sqrt_e_sinw        = params[j + 6],
              sqrt_e_cosw        = params[j + 7],
              tref,
              e, sqrt_e,
              sin_omega, cos_omega,
              a, sqrt_a,
              two_pi_over_period,
              E0;

            // Compute the period.
            a = sqrt_a_cosiy_cosix * sqrt_a_cosiy_cosix
              + sqrt_a_cosiy_sinix * sqrt_a_cosiy_sinix
              + sqrt_a_siniy       * sqrt_a_siniy;
            sqrt_a = sqrt(a);
            two_pi_over_period = sqrt((m_central + mass) * G_GRAV / (a*a*a));

            // Compute the eccentricity and special case circular orbits.
            e = sqrt_e_cosw * sqrt_e_cosw + sqrt_e_sinw * sqrt_e_sinw;
            if (e > DBL_EPSILON) {
                sqrt_e = sqrt(e);

                // omega
                sin_omega = sqrt_e_sinw / sqrt_e;
                cos_omega = sqrt_e_cosw / sqrt_e;

                // Compute the reference time using the time of transit.
                //  -> f = 0.5*pi - w
                //    -> sin(f) = cos(w)
                //    -> cos(f) = sin(w)
                //  -> tan(0.5*f) = sin(f) / (1 + cos(f))
                //                = cos(w) / (1 + sin(w))
                E0 = 2.0 * atan2(sqrt(1.0 - e) * cos_omega,
                                 sqrt(1.0 + e) * (1.0 + sin_omega));
                tref = t0 - (E0 - e * sin(E0)) / two_pi_over_period;
            } else {
                sin_omega = T(0.0);
                cos_omega = T(1.0);
                tref = t0 - T(0.5 * M_PI) / two_pi_over_period;
            }

            result[n]    = mass;                         // Mass
            result[n+1]  = exp(params[j+1]);             // Radius
            result[n+2]  = tref;                         // Reference time
            result[n+3]  = e;                            // Eccentricity
            result[n+4]  = cos_omega;                    // ...
            result[n+5]  = sin_omega;                    // ...
            result[n+6]  = a;                            // Semi-major axis
            result[n+7]  = sqrt_a_cosiy_cosix / sqrt_a;  // Inclination
            result[n+8]  = sqrt_a_cosiy_sinix / sqrt_a;  // ...
            result[n+9]  = sqrt_a_siniy       / sqrt_a;  // ...
            result[n+10] = two_pi_over_period;           // Time factor
        }

        return result;
    };

private:
    CentralBody central_;
    std::vector<Body> bodies_;

};  // class System

};  // namespace celeste

#endif  // _CELESTE_SYSTEM_H_
