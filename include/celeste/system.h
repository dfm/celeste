#ifndef _CELESTE_SYSTEM_H_
#define _CELESTE_SYSTEM_H_

#include <vector>

#include "celeste/central_body.h"
#include "celeste/body.h"

namespace celeste {

class System {
public:
    System () {};
    System (CentralBody central) : central_(central) {};
    void add_body (Body body) {
        bodies_.push_back(body);
    };

    size_t size () const {
        size_t count = central_.size();
        for (size_t i = 0; i < bodies_.size(); ++i) count += bodies_[i].size();
        return count;
    };

    CentralBody& central () { return central_; };
    Body& body (size_t index = 0) { return bodies_[index]; };

private:
    CentralBody central_;
    std::vector<Body> bodies_;

};  // class System

};  // namespace celeste

#endif  // _CELESTE_SYSTEM_H_
