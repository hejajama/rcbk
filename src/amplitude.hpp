/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _RAMPLITUDE_H
#define _RAMPLITUDE_H

#include <tools/config.hpp>
#include <vector>
#include <cmath>

enum InitialConditionR
{
    GBW,    // ref 0902.1112, no anomalous dimension
    MV,     // ref 0902.1112, anomalous dimension
    MV1,    // ref 0902.1112, same as MV but w.o. anomalous dimension
            //    no fitted parameters
    AN06    // ref e.g. 0704.012, 1_exp(-(rQ_s)^(2\gamma)/4)
};

// AmplitudeR::Initialize() must be called before this class is used, but first
// one needs to set up this (e.g. call SetInitialCondition etc).

class AmplitudeR
{
public:
    AmplitudeR();
    void Initialize();

    // Return tabulated value of amplitude
    REAL Ntable(int yind, int rind, int bind=0, int thetaind=0);



    // Returns index of the added rapidity
    int AddRapidity(REAL y);

    REAL InitialCondition(REAL r, REAL b);

    void AddDataPoint(int yind, int rind, int bind, int thetaind, REAL value);

    REAL RVal(int rind);
    REAL LogRVal(int rind);
    REAL BVal(int bind);
    REAL LogBVal(int bind);
    REAL ThetaVal(int thetaind);
    REAL YVal(int yind);
    int RPoints();
    int YPoints();
    int BPoints();
    int ThetaPoints();
    REAL MinR();
    REAL MaxR();
    REAL MaxLnR();
    REAL MinLnR();
    REAL RMultiplier();

    void SetMinR(REAL minr_);
    
    // Initial condition dependent \\alpha_s
    REAL Alpha_s_ic(REAL rsqr, REAL scaling=1.0);

    bool ImpactParameter(); // return bdep

    std::vector<REAL>& LogRVals();
    std::vector<REAL>& LogBVals();
    std::vector<REAL>& ThetaVals();
    
    InitialConditionR InitialCondition();
    void SetInitialCondition(InitialConditionR ic_);
    
    void SetAlphasScaling(REAL scaling);

    REAL InitialSaturationScaleSqr();


    // amplitude[yind][rind][bind][thetaind]
    std::vector < std::vector< std::vector< std::vector<REAL> > > > n;
private:
    std::vector<REAL> logrvals;
    std::vector<REAL> yvals;
    std::vector<REAL> logbvals;
    std::vector<REAL> rvals;
    std::vector<REAL> thetavals;
    bool bdep;      // do we take into account impact parameter dependency
    InitialConditionR ic;
    
    // Parameters which may depend on the IC
    REAL lambdaqcd2;
    REAL alphas_scaling;
    REAL Cfactorsqr;   // \alpha_s \sim 1/log(C^2/(r^2\lambdaqcd^2))
    REAL maxalphas;
    REAL minr;
    REAL Q_s0sqr;    // Initial saturation scale sqr
};

const REAL MINLN_N = -999;

#endif
