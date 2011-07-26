/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _RAMPLITUDE_H
#define _RAMPLITUDE_H

#include "config.hpp"
#include <vector>
#include <cmath>

class AmplitudeR
{
public:
    AmplitudeR();
    void Intialize();

    // Returns index of the added rapidity
    int AddRapidity(REAL y);

    REAL InitialCondition(REAL r, REAL b);

    void AddDataPoint(int yind, int rind, int bind, int thetaind, REAL value);

    REAL RVal(int rind);
    REAL LogRVal(int rind);
    REAL BVal(int bind);
    REAL LogBVal(int bind);
    REAL ThetaVal(int thetaind);
    int RPoints();
    int YPoints();
    int BPoints();
    int ThetaPoints();
    REAL MinR();
    REAL MaxR();
    REAL MaxLnR();
    REAL MinLnR();
    REAL RMultiplier();

    bool ImpactParameter(); // return bdep

    std::vector<REAL>& LogRVals();
    std::vector<REAL>& LogBVals();
    std::vector<REAL>& ThetaVals();


    // amplitude[yind][rind][bind][thetaind]
    std::vector < std::vector< std::vector< std::vector<REAL> > > > n;
private:
    std::vector<REAL> logrvals;
    std::vector<REAL> yvals;
    std::vector<REAL> logbvals;
    std::vector<REAL> thetavals;
    bool bdep;      // do we take into account impact parameter dependency
};

const REAL MINLN_N = -999;

#endif
