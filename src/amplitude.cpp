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
    Intialize();
}

void AmplitudeR::Intialize()
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
    const REAL Q_s0sqr = 0.24; // Fitted to HERA data at arXiv:0902.1112
    return 1.0 - exp(-SQR(r)*Q_s0sqr / 4.0 * std::exp( -SQR(b)/2 ) );
    //return 1.0 - std::exp(-10*SQR(r)*std::exp( -SQR(b)/2 ) );
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

int AmplitudeR::RPoints()
{
    return 155;
    //return 95;
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
    return 1e-5;
    //return 1e-6;
}

REAL AmplitudeR::RMultiplier()
{
    return 1.1;
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
