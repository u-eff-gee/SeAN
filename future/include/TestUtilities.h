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
    inline void is_equal(const T t, const U u) const {
        if(t != u){
            cout << __FILE__ << ": " << __FUNCTION__ << "( " << t << ", " << u << " )" << endl;
            throw;
        }
    }

    inline void is_close_relative(const double x, const double y, const double epsilon_relative) const {
        if(fabs((x-y)/(0.5*(x+y))) > epsilon_relative){
            cout << __FILE__ << ": " << __FUNCTION__ << "( " << x << ", " << y << ", " << epsilon_relative << " )" << endl;
            throw;
        }
    }

    inline void larger_than(const double x, const double y) const {
        if(y > x){
            cout << __FILE__ << ": " << __FUNCTION__ << "( " << x << ", " << y << " )" << endl;
            throw;            
        }
    }

    static constexpr double num_tol_rel = 1e-5;
};