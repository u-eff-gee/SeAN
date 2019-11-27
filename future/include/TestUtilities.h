#pragma once

#include <cmath>

using std::fabs;

class TestUtilities{

public:
    TestUtilities() = default;
    ~TestUtilities() = default;

    template<typename T, typename U>
    inline void is_equal(T t, U u) const {
        if(t != u)
            throw 0;
    }

    inline void is_close_absolute(double x, double y, double epsilon) const {
        if(fabs(x-y) > epsilon)
            throw 0;
    }

    inline void is_close_relative(double x, double y, double epsilon_relative) const {
        if(fabs((x-y)/(0.5*(x+y))) > epsilon_relative)
            throw 0;
    }
};