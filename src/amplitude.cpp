/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
#include <tools/tools.hpp>
#include <cmath>
#include <iostream>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

using std::cout; using std::endl;

AmplitudeR::AmplitudeR()
{
    bdep=false;
    alphas_scaling=1.0;
    minr=1e-9;
    SetInitialCondition(GBW);
    Q_s0sqr=0;
    
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
        rvals.push_back(std::exp(logrvals[rind]));
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

double Inthelperf_mvinfra(double k, void* p);
struct Inthelper_mvinfra{ double lambdaqcd2; double qsqr; double r; };
REAL AmplitudeR::InitialCondition(REAL r, REAL b)
{
	const REAL e = 2.7182818;
    if (ic == GBW)
    {
        if (r<3e-6) return SQR(r)*Q_s0sqr / 4.0 * std::exp( -SQR(b)/2 );
        return 1.0 - std::exp(-SQR(r)*Q_s0sqr / 4.0 * std::exp( -SQR(b)/2 ) );
    }
    if (ic == MV)
    {   // same ref as for GBW
        const REAL anomalous_dimension = 1.13;
        if (r < 2e-6)   ///NOTE: factor 1/4 "correctly", not as in ref.
            return std::pow(SQR(r)*Q_s0sqr, anomalous_dimension)/4.0
            * std::log( 1.0/(r*std::sqrt(lambdaqcd2)) + e);
        return 1.0 - std::exp(-std::pow(SQR(r)*Q_s0sqr, anomalous_dimension)/4.0
            * std::log( 1.0/(r*std::sqrt(lambdaqcd2)) + e) );
    }
    if (ic == MV_Au)
    {
		// Exponentaited proton  => nucleus
		// N(r,b) = 1 - exp(-\sigma_0 T_A(b) N(r) ), where
		// N(r) is the MV model (or any)
		const double anomalous_dimension = 1.13;
		
		/// These coefficients depend on A and dipole model!
		const double A=197;
		const double b = 0.172*FMGEV;	 	// Impact parameter, 10% most central
		const double sigma_pp = 32.77 * 0.1 * FMGEV*FMGEV;		// convert mb in 1/GeV^2
		double np = 1.0 - std::exp(-std::pow(SQR(r)*Q_s0sqr, anomalous_dimension)/4.0
            * std::log( 1.0/(r*std::sqrt(lambdaqcd2)) + e) );
            
         
        return 1.0 - std::exp( -sigma_pp * A * T_A(b, A) * np);
        
		
		
	}
    if (ic == MV1 or ic == MV1_dAu)
    {
        // Same as previoius but w.o. anomalous dimension
        if (r < 2e-6)
            return SQR(r)*Q_s0sqr/4.0
                * std::log( 1.0/(r*std::sqrt(lambdaqcd2)) + e);
        return 1.0 - std::exp(-SQR(r)*Q_s0sqr/4.0
            * std::log( 1.0/(r*std::sqrt(lambdaqcd2)) + e) );
    }
    if (ic == AN06)
    {
        const REAL gamma = 0.6;
        if (r<1e-10) return std::pow(SQR(r)*Q_s0sqr,gamma)/4.0;
        return 1.0 - std::exp( -std::pow( SQR(r)*Q_s0sqr, gamma )/4.0 );
    }
    if (ic == MV1_OSC)
    {
        Inthelper_mvinfra helper;
        helper.qsqr = Q_s0sqr; helper.lambdaqcd2 = lambdaqcd2;
        helper.r=r;
        gsl_function f;
        f.params=&helper;
        f.function = Inthelperf_mvinfra;
        REAL result, abserr;
        gsl_integration_workspace *workspace 
          = gsl_integration_workspace_alloc(500);
        int status = gsl_integration_qag(&f, 0.0001, 100,
            0, 0.01, 100, GSL_INTEG_GAUSS51, workspace, &result, &abserr);
        gsl_integration_workspace_free(workspace);
        if (status)
        {
            cerr << "icint failed at " << LINEINFO << " intresult: " << result
                    << " relerr " << std::abs(abserr/result) << " r: " << r << endl;
        }
        if (result < 1e-6)
            return result;
        else
            return 1.0 - std::exp(-result);
    }
    cerr << "Unkown initial condition set! " << LINEINFO << endl;
    return 0;
}

REAL Inthelperf_mvinfra(double k, void* p)
{
    Inthelper_mvinfra * par = (Inthelper_mvinfra*) p;
    REAL result=0;

    result = SQR(par->r)*par->qsqr
        *(1.0-gsl_sf_bessel_J0(SQR(par->r)*par->lambdaqcd2))/(M_PI*k*k*k)
        * std::atan(SQR(k)/(SQR(par->r)*par->lambdaqcd2))*2/(M_PI);
    return result;
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
 * if scaling is given, it overrides saved alphas_scaling
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
    return minr;
    //return 2e-6;  // kw mv
    //return 8e-7;  // kw gbw
    //return 5e-7;  // balitsky mv
    //return 5e-8;    // balitsky gbw
    //return 1e-8;
    //return 1e-9;
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
    return rvals[rind];
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
        case GBW:
            Q_s0sqr = 0.24;
            lambdaqcd2=0.241*0.241;
            Cfactorsqr=4.0;
            maxalphas=0.7;
            x0=0.01;
            break;
        case MV_Au:	 // Exponentiated proton
        case MV:
            Q_s0sqr = 0.15;
            lambdaqcd2=0.241*0.241;
            Cfactorsqr=4.0;
            maxalphas=0.7;
            x0=0.01;
            break;
        case MV1:
            // Values from the fit to the HERA data 0902.1112
            Q_s0sqr = 0.2;
            lambdaqcd2=0.241*0.241;
            Cfactorsqr=4.0;
            maxalphas=0.7;
            x0=0.01;
            break;
        case AN06:
            // ref. 0704.0612
            Q_s0sqr = 1.0;
            lambdaqcd2 = 0.2*0.2;
            Cfactorsqr = 1.0;
            maxalphas=0.5;
            x0=0.01;
            break;
        case MV1_dAu:   // ref 1001.1378 , 1009.3215
            Q_s0sqr = 0.6 ; //0.4;  // 0.6: Most central
            lambdaqcd2 = 0.241*0.241;
            x0=0.02;
            maxalphas=0.7;
            Cfactorsqr=4.0;
            break;
        case MV1_OSC:
            Q_s0sqr = 2;
            lambdaqcd2 = 0.241*0.241;
            x0=0.01;
            maxalphas=0.7;
            Cfactorsqr = 0.4;
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

/*
 * Set smallest r for the grid
 * Must be set before Initialize() is called!
 */
void AmplitudeR::SetMinR(REAL minr_)
{
    if (rvals.size()!=0)
    {
        cerr << "SetMinR() called after AmplitudeR::Initialize() is called, "
        << "can't do anything... " << LINEINFO << endl;
        return;
    }
    minr=minr_;
}

REAL AmplitudeR::InitialSaturationScaleSqr()
{
    return Q_s0sqr;
}

REAL AmplitudeR::X0()
{
    return x0;
}
