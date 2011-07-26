/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "tools.hpp"
#include "amplitude.hpp"
#include "solver.hpp"
#include <iostream>
#include <gsl/gsl_errno.h>

int main(int argc, char* argv[])
{
    AmplitudeR N;
    Solver s(&N);

    gsl_set_error_handler(&ErrHandler);

    s.Solve(2);




    return 0;
}
