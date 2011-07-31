/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
#include <cmath>
#include <iostream>
using std::cout; using std::endl;

AmplitudeR::AmplitudeR()
{
    bdep=false;
    alphas_scaling=1.0;
    SetInitialCondition(IPSAT);
    
}

void AmplitudeR::Initialize()
{
    if (ImpactParameter())
    {
        for (int thetaind=0; thetaind<ThetaPoints(); thetaind++)
        {
            REAL tmptheta = thetaind*2.0*M_PI / static_cast<REAL>(ThetaPoints()-1.0);
            thetavals.push_back( tmptheta );
        }
    }
    else
    {
        thetavals.push_back(0);
        logbvals.push_back(MINLN_N);
    }

    for (int rind=0; rind < RPoints(); rind++)
    {
        logrvals.push_back(std::log(MinR()
                * std::pow(RMultiplier(), rind) ) );
        if (ImpactParameter()) logbvals.push_back(logrvals[rind]);
    }
    if (logrvals[logrvals.size()-1]<0)
    {
        cerr << "MaxR must be at least " << std::exp(logrvals[logrvals.size()-1])
            << endl;
        return;
    }

    AddRapidity(0);

    for (int rind=0; rind<RPoints(); rind++)
    {
        for (int bind=0; bind<BPoints(); bind++)
        {
            for (int thetaind=0; thetaind < ThetaPoints(); thetaind++)
            {
                n[0][rind][bind][thetaind]
                    = InitialCondition(RVal(rind), BVal(bind) );
            }
        }
    }

}


/*
 * Add new rapidity value for tables
 * Returns index of the added rapidity, -1 if error occurs
 */
int AmplitudeR::AddRapidity(REAL y)
{
    if (yvals.size()>0)
    {
        if (y <= yvals[yvals.size()-1])
        {
            cerr << "Trying to add rapidity value " << y << " but currently "
            << "the largest rapidity tabulated is " << yvals[yvals.size()-1] << endl;
            return -1;
        }
    }

    yvals.push_back(y);

    // n[yind][rind][bind][thetaind]
    std::vector< std::vector< std::vector< REAL> > > tmprvec;
    for (int rind=0; rind<RPoints(); rind++)
    {
        std::vector< std::vector<REAL> > tmpbvec;
        for (int bind=0; bind<BPoints(); bind++)
        {
            std::vector<REAL> tmpthetavec;
            for (int thetaind=0; thetaind < ThetaPoints(); thetaind++)
            {
                tmpthetavec.push_back(0);
            }
            tmpbvec.push_back(tmpthetavec);
        }
        tmprvec.push_back(tmpbvec);
    }
    n.push_back(tmprvec);

    return yvals.size()-1;
}


REAL AmplitudeR::InitialCondition(REAL r, REAL b)
{
    if (ic == IPSAT)
    {
        const REAL Q_s0sqr = 0.24; // Fitted to HERA data at arXiv:0902.1112
        return 1.0 - std::exp(-SQR(r)*Q_s0sqr / 4.0 * std::exp( -SQR(b)/2 ) );
    }
    if (ic == AN06)
    {
        const REAL Q_s0sqr = 1.0;   // following arXiv:0704.0612, not fitted
        const REAL gamma = 0.6;
        return 1.0 - std::exp( -std::pow( SQR(r)*Q_s0sqr, gamma )/4.0 );
    }
    cerr << "Unkown initial condition set! " << LINEINFO << endl;
    return 0;
}


void AmplitudeR::AddDataPoint(int yind, int rind, int bind,
    int thetaind, REAL value)
{
    n[yind][rind][bind][thetaind] = value;
}

/*
 * Return tabulated value of amplitude
 */
REAL AmplitudeR::Ntable(int yind, int rind, int bind, int thetaind)
{
    return n[yind][rind][bind][thetaind];
}

/*
 * Initial condition dependent strong coupling constant
 * \lambda_{QCD}^2 and the scaling factor (whether or not there is factor
 * 4 in the ln[1/(r^2\lambdaQCD^2)] term) may depend on the IC
 * if alphas_scaling is given, it overrides saved alphas_scaling
 * 
 * SetInitialCondition sets appropriate values to required variables
 */
REAL AmplitudeR::Alpha_s_ic(REAL rsqr, REAL scaling)
{
    REAL scalefactor=0;
    if (std::abs(scaling-1.0)>0.0001)   // Don't use stored value
        scalefactor = Cfactorsqr*scaling;
    else
        scalefactor = Cfactorsqr*alphas_scaling;
    
    if (scalefactor/(rsqr*lambdaqcd2) < 1.0) return maxalphas;
    
    REAL alpha = 12.0*M_PI/( (33.0-2.0*Nf)*std::log(scalefactor/ (rsqr*lambdaqcd2) ) );
    if (alpha>maxalphas)
        return maxalphas;
    
    return alpha;

}
int AmplitudeR::RPoints()
{
    return 400;
    //return 90;
}

int AmplitudeR::YPoints()
{
    return yvals.size();
}

int AmplitudeR::BPoints()
{
    if (!ImpactParameter()) return 1;
    return RPoints();
}

int AmplitudeR::ThetaPoints()
{
    if (!ImpactParameter()) return 1;
    return 10;
}

REAL AmplitudeR::MinR()
{
    //return 1e-7;
    return 1e-9;
    //return 1e-5;
}

REAL AmplitudeR::RMultiplier()
{
    //return 1.08;
    return std::pow(50.0/MinR(), 1.0/(RPoints()-1));
}

REAL AmplitudeR::MaxR()
{
    return MinR()*std::pow(RMultiplier(), RPoints()-1);
}

REAL AmplitudeR::MaxLnR()
{
    return logrvals[logrvals.size()-1];
}

REAL AmplitudeR::MinLnR()
{
    return logrvals[0];
}

REAL AmplitudeR::RVal(int rind)
{
    return std::exp( logrvals[rind] );
}

REAL AmplitudeR::LogRVal(int rind)
{
    return logrvals[rind];
}

REAL AmplitudeR::BVal(int bind)
{
    return std::exp( logbvals[bind] );
}

REAL AmplitudeR::LogBVal(int rind)
{
    return logbvals[rind];
}

REAL AmplitudeR::ThetaVal(int thetaind)
{
    return thetavals[thetaind];
}

REAL AmplitudeR::YVal(int yind)
{
    return yvals[yind];
}

std::vector<REAL>& AmplitudeR::LogRVals()
{
    return logrvals;
}

std::vector<REAL>& AmplitudeR::LogBVals()
{
    return logbvals;
}

std::vector<REAL>& AmplitudeR::ThetaVals()
{
    return thetavals;
}

bool AmplitudeR::ImpactParameter()
{
    return bdep;
}

void AmplitudeR::SetInitialCondition(InitialConditionR ic_)
{
    ic=ic_;
    switch(ic)
    {
        case IPSAT:
            // Values from the fit to the HERA data 0902.1112
            lambdaqcd2=0.241*0.241;
            Cfactorsqr=4.0;
            maxalphas=0.7;
            break;
        case AN06:
            // ref. 0704.0612
            lambdaqcd2 = 0.2*0.2;
            Cfactorsqr = 1.0;
            maxalphas=0.5;
            break;
    };

}

InitialConditionR AmplitudeR::InitialCondition()
{
    return ic;
}

void AmplitudeR::SetAlphasScaling(REAL scaling)
{
    alphas_scaling=scaling;
}


