/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

extern "C"
{
    #include "fourier/fourier.h"
}
#include "ugd.hpp"
#include <gsl/gsl_sf_bessel.h>

UGD::UGD(AmplitudeLib* N_)
{
    N=N_;

    // Some initialisation stuff -
    set_fpu_state();
    init_workspace_fourier(1500);   // number of bessel zeroes, max 2000
}


/*
 * Evaluate unintegrated gluon density at point k,y
 * \psi(k,y) = C_F/\alpha_s(k) 1/(2\pi)^3 \int d^2 r exp(-ik.r)
 *              * \nabla_r^2 (2N(r,y) - N(r,y)^2 )
 * \nabla_r^2 = 1/r \partial_r + \partial_r^2
 */

struct HankelHelper
{
    REAL y;
    REAL k;
    AmplitudeLib* N;
};

REAL Hankelhelperf(REAL r, void* p)
{
    HankelHelper* par = (HankelHelper*)p;
    if (r<par->N->MinR() or r>par->N->MaxR()) return 0;
    REAL result = r*( 1.0/r*par->N->N(r, par->y, 1) + par->N->N(r, par->y, 2)
                    - 1.0/r*par->N->N(r, par->y, 1, true) - par->N->N(r, par->y, 2, true)
     );
    return result;
}
REAL UGD::Evaluate(REAL k, REAL y)
{
    HankelHelper par;
    par.y=y; par.N=N; par.k=k;
    REAL result = fourier_j0(k,Hankelhelperf,&par);
    
    REAL Cf = (SQR(Nc)-1)/(2.0*Nc);
    result *= Cf/Alpha_s(SQR(k));  //Todo: scaling?
    result /= SQR(2.0*M_PI);
    return result;
}
