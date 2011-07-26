/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _SOLVER_H
#define _SOLVER_H

#include "amplitude.hpp"

class Solver
{
    public:
        Solver(AmplitudeR* N_);
        void Solve(REAL maxy);

        REAL RapidityDerivative(REAL y, REAL lnr01, REAL lnb01, REAL thetab, const REAL* data);

        REAL Kernel(REAL r01, REAL r02, REAL r12, REAL y=0,
            REAL b01=0, REAL thetab=0, REAL theta2=0 );
        REAL InterpolateN(REAL lnr, REAL lnb, REAL thetab, const REAL* data);

    private:
        AmplitudeR* N;

};

const int INTERPOLATION_POINTS = 6;

#endif
