#ifndef _CELESTE_SYSTEM_H_
#define _CELESTE_SYSTEM_H_

#include <cmath>
#include <cfloat>
#include <vector>

#include "celeste/central_body.h"
#include "celeste/body.h"

namespace celeste {

class System {
public:
    System (CentralBody& central) : central_body_(central) {};
    void add_body (Body body) {
        body.set_central(&central_body_);
        orbiting_bodies_.push_back(body);
    };

private:
    CentralBody central_body_;
    std::vector<Body> orbiting_bodies_;

};  // class System

};  // namespace celeste

#endif  // _CELESTE_SYSTEM_H_
