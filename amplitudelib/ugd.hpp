/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _UGD_H
#define _UGD_H

#include "../src/config.hpp"
#include "amplitudelib.hpp"

/*
 * Unintegrated gluon distribution
 */
 
class UGD
{
    public:
        UGD(AmplitudeLib* N_);

        REAL Evaluate(REAL k, REAL y);

    private:
        AmplitudeLib* N;
        
};

#endif
