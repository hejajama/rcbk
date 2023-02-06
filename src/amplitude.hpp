/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _RAMPLITUDE_H
#define _RAMPLITUDE_H

#include <vector>
#include <cmath>
#include "ic.hpp"
#include "config.hpp"

// AmplitudeR::Initialize() must be called before this class is used, but first
// one needs to set up this (e.g. call SetInitialCondition etc).

// Whether we have fixed number of flavours (nf=Nf=3 defined in config.hpp) or
// we include heavy quarks and use the fact that number of active flavors is changing
// see 1012.4408 Sec. 2.2
enum AlphasFlavours
{
    LIGHTQ,
    HEAVYQ    // Include also heavy flavours (c and b)
};

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
    std::string Alpha_s_str();

    bool ImpactParameter(); // return bdep

    std::vector<REAL>& LogRVals();
    std::vector<REAL>& LogBVals();
    std::vector<REAL>& ThetaVals();
    
    InitialCondition* GetInitialCondition();
    void SetInitialCondition(InitialCondition* ic_);
    
    void SetAlphasScaling(REAL scaling);
    double GetAlphasScaling();

    void SetAlphasFlavours(AlphasFlavours f);
    AlphasFlavours GetAlphasFlavours();

    // Lambda_QCD values in case of fixed number of flavors
    void SetLambdaQcd(double lambda);
    double GetLambdaQcd();	// in GeV
    
    void SetAlphasFreeze(REAL c);
    double GetAlphasFreeze();

    void SetRPoints(int rp);


    // amplitude[yind][rind][bind][thetaind]
    std::vector < std::vector< std::vector< std::vector<REAL> > > > n;
private:
    std::vector<REAL> logrvals;
    std::vector<REAL> yvals;
    std::vector<REAL> logbvals;
    std::vector<REAL> rvals;
    std::vector<REAL> thetavals;
    bool bdep;      // do we take into account impact parameter dependency
    
    InitialCondition *initial_condition;
    
    // Parameters which may depend on the IC
    REAL Csqr;   // \alpha_s \sim 1/log(4 C^2/(r^2\lambdaqcd^2))
    REAL maxalphas;
    REAL alphas_freeze_c;	// describes how smoothly alphas is cutoff at maxalphas, if =0, then sharp cutoff at maxalphas
    REAL minr;
    double lambdaqcd;
    
    int RPOINTS;

    AlphasFlavours alphas_flavours;
};

const REAL MINLN_N = -999;
/// Accuracy settings
const double MAXR = 100; // orig: 50


#endif
