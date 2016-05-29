#ifndef _CELESTE_PARAMETER_H_
#define _CELESTE_PARAMETER_H_

namespace celeste {

#define COUNT_PARAM(NAME) if (!NAME.is_frozen()) count++;
#define GET_PARAM(NAME)   if (!NAME.is_frozen()) params[count++] = NAME;
#define SET_PARAM(NAME)   if (!NAME.is_frozen()) NAME = params[count++];

template <typename T>
class Parameter {
public:
    Parameter (T value = T(0.0), bool frozen = true) : value_(value), frozen_(frozen) {};

    void thaw () { frozen_ = false; };
    void freeze () { frozen_ = true; };
    bool is_frozen () const { return frozen_; };

    operator T() const { return value_; };

    void operator= (const T& value) { set_value(value); };
    void set_value (const T& value) { value_ = T(value); };
    T get_value () const { return value_; };

private:
    bool frozen_;
    T value_;

};  // class Parameter

};  // namespace celeste

#endif  // _CELESTE_PARAMETER_H_
