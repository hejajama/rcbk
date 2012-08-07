/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _RAMPLITUDE_H
#define _RAMPLITUDE_H

#include <tools/config.hpp>
#include <vector>
#include <cmath>
#include "ic.hpp"

enum InitialConditionR
{
    GBW,    // ref 0902.1112, no anomalous dimension
    MVic,     // ref 0902.1112, anomalous dimension
    MV_Au,	// Glauber-type expenentiated MVgamma
    MV1,    // ref 0902.1112, same as MV but w.o. anomalous dimension
            //    no fitted parameters
    AN06,    // ref e.g. 0704.012, 1_exp(-(rQ_s)^(2\gamma)/4)
    MV1_dAu, // ref 1001.1378, RHIC dAu data
    MV1_OSC // ref 0708.0231, infrared-oscillation-cancelling factor
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
    
    InitialCondition* GetInitialCondition();
    void SetInitialCondition(InitialCondition* ic_);
    
    void SetAlphasScaling(REAL scaling);
    double GetAlphasScaling();
    void SetLambdaQcd(double lambda);
    double GetLambdaQcd();	// in GeV


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
    
    InitialCondition *initial_condition;
    
    // Parameters which may depend on the IC
    REAL Csqr;   // \alpha_s \sim 1/log(4 C^2/(r^2\lambdaqcd^2))
    REAL maxalphas;
    REAL minr;
    double lambdaqcd;
};

const REAL MINLN_N = -999;

#endif
