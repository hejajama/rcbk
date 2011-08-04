/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _AMPLITUDELIB_H
#define _AMPLITUDELIB_H

#include "../src/config.hpp"
#include "../src/interpolation.hpp"
#include <vector>

class AmplitudeLib
{
    public:
        AmplitudeLib(std::string datafile);
        ~AmplitudeLib();

        // der w.r.t r der times.
        REAL N(REAL r, REAL y, int der=0);
        
        // Amplitude in k-space 
        REAL N_k(REAL kt, REAL y);

        // Unintegrated gluon density
        REAL UGD(REAL k, REAL y);

        // k_T factorization: d\sigma / (dyd^2p_T)
        // = const * 1/p_T^2 \int d^2 k_T/4 \alphas_(Q) \psi(|p_t+k_T|/2,x1)
        //   * \psi(|p_t-k_T|/2, x2)
        REAL dSigmadyd2pt(REAL pt, REAL x1, REAL x2);
        REAL dSigmady(REAL y, REAL sqrts);

        // d ln N / d ln r^2
        REAL LogLogDerivative(REAL r, REAL y);

        // Saturation scale N(r, y) = Ns
        REAL SaturationScale(REAL y, REAL Ns);

        void InitializeInterpoaltion(REAL y);
    
        int YVals();
        int RPoints();
        REAL MinR();
        REAL MaxR();
        REAL MaxY();
        
        
    private:
        std::vector< std::vector<REAL> > n;
        std::vector<REAL> yvals;
        std::vector<REAL> lnrvals;
        std::vector<REAL> rvals;
        Interpolator *interpolator;

        REAL interpolator_y;
        REAL* tmprarray;
        REAL* tmpnarray;

        REAL minr;
        REAL rmultiplier;
        int rpoints;
};

const int INTERPOLATION_POINTS = 8;
const REAL UGD_IR_CUTOFF=0.3;   // ugd(k<UGD_IR_CUTOFF)=0     BAD?????

#endif
