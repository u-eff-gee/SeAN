#pragma once

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;
using std::fabs;

class TestUtilities{

public:
    TestUtilities() = default;
    ~TestUtilities() = default;

    template<typename T, typename U>
    inline void is_equal(T t, U u) const {
        if(t != u)
            throw;
    }

    inline void is_close_absolute(double x, double y, double epsilon) const {
        if(fabs(x-y) > epsilon)
            throw;
    }

    inline void is_close_relative(double x, double y, double epsilon_relative) const {
        if(fabs((x-y)/(0.5*(x+y))) > epsilon_relative){
            cout << __FILE__ << ": " << __FUNCTION__ << "( " << x << ", " << y << ", " << epsilon_relative << " )" << endl;
            throw;
        }
    }

    static constexpr double num_tol_rel = 1e-5;
};