#ifndef _CELESTE_PARAMETER_H_
#define _CELESTE_PARAMETER_H_

namespace celeste {

#define INVALID_PARAMETER_INDEX -count

#define GET_PARAM_FUNCTION(BODY)                      \
    const Parameter& parameter (size_t index) const { \
        size_t count = 0;                             \
        BODY                                          \
        throw INVALID_PARAMETER_INDEX;                \
    };
#define GET_PARAM(NAME)  if (count++ == index) return NAME;

#define SET_PARAM_FUNCTION(BODY)                  \
    void parameter (size_t index, double value) { \
        size_t count = 0;                         \
        BODY                                      \
        throw INVALID_PARAMETER_INDEX;            \
    };
#define SET_PARAM(NAME)   if (count++ == index) { NAME = value; return; }


class Parameter {
public:
    Parameter (double value = 0.0, bool frozen = true) : value_(value), frozen_(frozen) {};

    void thaw () { frozen_ = false; };
    void freeze () { frozen_ = true; };
    bool is_frozen () const { return frozen_; };

    operator double() const { return value_; };

    void operator= (const double& value) { set_value(value); };
    void set_value (const double& value) { value_ = double(value); };
    double get_value () const { return value_; };

private:
    bool frozen_;
    double value_;

};  // class Parameter

class ParameterizedModel {
public:
    virtual size_t full_size () const = 0;
    size_t size () const {
        size_t count = 0, n = full_size();
        for (size_t i = 0; i < n; ++i)
            if (!parameter(i).is_frozen()) count++;
        return count;
    };

    virtual const Parameter& parameter (size_t index) const = 0;
    virtual void parameter (size_t index, double value) = 0;

    size_t get_parameter_vector (double* params) const {
        size_t count = 0, n = full_size();
        for (size_t i = 0; i < n; ++i) {
            const Parameter& par = parameter(i);
            if (!par.is_frozen()) params[count++] = par;
        }
        return count;
    };
    size_t set_parameter_vector (const double* params) {
        size_t count = 0, n = full_size();
        for (size_t i = 0; i < n; ++i)
            if (!parameter(i).is_frozen()) parameter(i, params[count++]);
        return count;
    };

    size_t get_full_vector (double* params) const {
        size_t count = 0, n = full_size();
        for (size_t i = 0; i < n; ++i) params[count++] = parameter(i);
        return count;
    };

};  // class ParameterizedModel

};  // namespace celeste

#endif  // _CELESTE_PARAMETER_H_
