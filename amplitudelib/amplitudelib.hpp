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

        // der w.r.t r der times. if sqr, then calculate N^2, not N
        REAL N(REAL r, REAL y, int der=0, bool sqr=false);
        
        // Amplitude in k-space 
        REAL N_k(REAL kt, REAL y);

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
#endif
