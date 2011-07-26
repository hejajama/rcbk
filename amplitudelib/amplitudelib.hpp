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

        REAL N(REAL r, REAL y);
    
        int YVals();
        int RPoints();
        REAL MinR();
        REAL MaxR();
        
        
    private:
        std::vector< std::vector<REAL> > n;
        std::vector<REAL> yvals;
        std::vector<REAL> lnrvals;

        REAL minr;
        REAL rmultiplier;
        int rpoints;
};

const int INTERPOLATION_POINTS = 8;
#endif
