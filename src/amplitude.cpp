/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
#include <tools/tools.hpp>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

using std::cout; using std::endl;

AmplitudeR::AmplitudeR()
{
    bdep=false;
    Csqr=1.0;
    minr=1e-9;
    lambdaqcd=0.241;
	maxalphas=0.7;
	alphas_freeze_c = 0;
    
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


    AddRapidity(0);

    for (int rind=0; rind<RPoints(); rind++)
    {
        for (int bind=0; bind<BPoints(); bind++)
        {
            for (int thetaind=0; thetaind < ThetaPoints(); thetaind++)
            {
                //n[0][rind][bind][thetaind]
                //    = InitialCondition(RVal(rind), BVal(bind) );
                n[0][rind][bind][thetaind] 
					= initial_condition->DipoleAmplitude(RVal(rind), BVal(bind));
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
 * \lambda_{QCD}^2 and the scaling factor alphas_scaling
 * \alpha_s ~ 1/(log (4C^2/(r^2 lambdaqcd))), and C^2=Csqr
 * if scaling is given, it overrides saved alphas_scaling
 */
const double alphas_mu0=2.5;
//double alphas_scl = 1.26095;
REAL AmplitudeR::Alpha_s_ic(REAL rsqr, REAL scaling)
{
	if (std::abs(scaling-1.0)>0.0001)
		cerr << "You are using alpas-scaling, are you sure???? Check code! " << LINEINFO << endl;

	///TODO! VÄLIAIKAINEN T.L. ANALYYSIIN
	
	if (alphas_freeze_c < 0.00001)  // sharp cutoff at maxalphas
	{
	    double scalefactor=0;
		if (std::abs(scaling-1.0)>0.0001)   // Don't use stored value
			scalefactor = scaling;
		else
			scalefactor = 4.0*Csqr;
		
		if (scalefactor/(rsqr*lambdaqcd*lambdaqcd) < 1.0) return maxalphas;
		double alpha = 12.0*M_PI/( (33.0-2.0*Nf)*std::log(scalefactor/ (rsqr*lambdaqcd*lambdaqcd) ) );
		if (alpha>maxalphas)
			return maxalphas;
		return alpha;
	}

	// Smooth cutoff
	return 4.0*M_PI / ( 9.0 * std::log(
		std::pow( std::pow(alphas_mu0, 2.0/alphas_freeze_c) + std::pow(4.0*Csqr/(rsqr*lambdaqcd*lambdaqcd), 1.0/alphas_freeze_c), alphas_freeze_c)	
		) );
		
/*		
       ///TODO: Väliaikainen: ipsat
       //double musqr = 4.0/rsqr + 1.1699999;
       //double ipsatalphas= 12.0*M_PI/( (33.0-2.0*Nf)*std::log(musqr/(lambdaqcd*lambdaqcd) ));
       

    
    return alpha;
*/
}

std::string AmplitudeR::Alpha_s_str()
{
	std::stringstream ss;
	if (alphas_freeze_c < 0.0001)
		ss << "\\alpha_s = 12 \\pi / [ (33.0 - 2.0*Nf) * log(4.0*C^2/(r^2*lambdaqcd^2) ], C^2=" << Csqr << ", lambdaqcd=" << lambdaqcd << " GeV, maxalphas=" << maxalphas <<", Nf=" << Nf;
	else
		ss << "\\alpha_s = 12 \\pi / [ (33.0 - 2.0*Nf) * log[ ( mu0^(1/c) + (4C^2/(lambdaqcd^2*r^2))^(1/c) )^c) ] ], C^2=" << Csqr <<", lambdaqcd=" << lambdaqcd << " GeV, mu0=" << alphas_mu0 <<", Nf=" << Nf;
	return ss.str();
}


int AmplitudeR::RPoints()
{
    return 400;
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
}

REAL AmplitudeR::RMultiplier()
{
    //return 1.08;
    double max = 50;
    //if (max > initial_condition->MaxR())
	//	max = 0.9999 * initial_condition->MaxR();
    return std::pow(max/MinR(), 1.0/(RPoints()-1));
}

REAL AmplitudeR::MaxR()
{
    double max = MinR()*std::pow(RMultiplier(), RPoints()-1);
    //if (max > initial_condition->MaxR())
	//	max = 0.9999 * initial_condition->MaxR();
	return max;
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

void AmplitudeR::SetInitialCondition(InitialCondition* ic_)
{
	initial_condition=ic_;
}

void AmplitudeR::SetAlphasScaling(REAL scaling)
{
    Csqr=scaling;
}

double AmplitudeR::GetAlphasScaling()
{
	return Csqr;
}

double AmplitudeR::GetAlphasFreeze()
{
	return alphas_freeze_c;
}

void AmplitudeR::SetAlphasFreeze(REAL c)
{
	alphas_freeze_c = c;
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
    
	if (minr_ < initial_condition->MinR())
	{
		cout << "# NOTICE: smallest r allowed by IC is " << initial_condition->MinR()
			<< ": most likely you are using datafile IC, and N at smaller r is set to 0" << endl;
		/*cerr << "Can't set minr to " << minr_ << " as the smallest r allowed "
		<< "by the IC is " << initial_condition->MinR() 
		<<", setting minr=" << initial_condition->MinR()  << endl;
		minr_ = initial_condition->MinR()*1.000000001;*/
	}
	
    minr=minr_;
}

InitialCondition* AmplitudeR::GetInitialCondition()
{
	return initial_condition;
}

void AmplitudeR::SetLambdaQcd(double lambda)
{
	lambdaqcd=lambda;
}

double AmplitudeR::GetLambdaQcd()
{
	return lambdaqcd;
}
