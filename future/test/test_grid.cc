#include "Grid.h"
#include "TestUtilities.h"

int main(){

    Grid grid;
    TestUtilities test;

    double range = 5e4; // An extremely large range needs to be chosen for the required precision.
    test.test_closeness(grid.breit_wigner_coverage(0., 1., -range, range), 
        1., test.num_tol_rel);
    vector<double> x_values = grid.breit_wigner(0., 1., -range, range, 5);
    test.test_equality<size_t, size_t>(x_values.size(), 5);
    test.test_equality<double, double>(x_values[0], -range);
    test.test_closeness(x_values[1], -0.5, test.num_tol_rel);
    test.test_closeness(x_values[2], 0., test.num_tol_rel);
    test.test_closeness(x_values[3], 0.5, test.num_tol_rel);
    test.test_equality<double, double>(x_values[4],  range);

    range = 5.;
     test.test_closeness(grid.normal_coverage(0., 1., -range, range), 
        1., test.num_tol_rel);
    x_values = grid.normal(0., 1., -range, range, 5);
    test.test_equality<size_t, size_t>(x_values.size(), 5);
    test.test_equality<double, double>(x_values[0], -range);
    test.test_closeness(x_values[1], -0.6744897501960817, test.num_tol_rel);
    test.test_closeness(x_values[2], 0., test.num_tol_rel);
    test.test_closeness(x_values[3], 0.6744897501960817, test.num_tol_rel);
    test.test_equality<double, double>(x_values[4],  range);

    // Test strictly_increasing function which makes each value unique and sorts the array.
    vector<double> x{1., 0., 4., 4., 4., 4., 2., 2., 3.};
    grid.strictly_increasing(x);
    test.test_equality<double, double>(x[0], 0.  );
    test.test_equality<double, double>(x[1], 1.  );
    test.test_equality<double, double>(x[2], 2.  );
    test.test_equality<double, double>(x[3], 2.5 );
    test.test_equality<double, double>(x[4], 3.  );
    test.test_equality<double, double>(x[5], 3.25);
    test.test_equality<double, double>(x[6], 3.5 );
    test.test_equality<double, double>(x[7], 3.75);
    test.test_equality<double, double>(x[8], 4.  );
}