/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitude.hpp"
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

using std::cout; using std::endl;
// Heavy quark masses, m_c = 1.27 GeV; m_b=4.2 GeV
/// TODO: Define somewhere else
double heavyqmasses[2] = { 1.27, 4.2 } ;

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

    double lqcd=lambdaqcd;
    double nf = Nf;
    double b0 = 33.0/3.0 - 2.0/3.0*nf;

    /* Varying n_f scheme (heavy quarks are included), compute effective Lambda_QCD
     * (such that alphas(r) is continuous), see 1012.4408 sec. 2.2.
     */
     
    if (alphas_flavours == HEAVYQ)
    {
        cerr << "Check piis ja b0:t!! " << LINEINFO << endl;

        double dipolescale = 4.0*Csqr / rsqr;

        if (dipolescale < SQR(heavyqmasses[0]))
            nf=3;
        else if (dipolescale < SQR(heavyqmasses[1]))
            nf=4;
        else
            nf = 5;

        b0 = 11.0/3.0 - 2.0/3.0*nf;
        // Now we compute "effective" Lambda by requiring that we get the experimental value for alpha_s
        // at the Z0 mass, alphas(Z0)=0.1184, m(Z0)=91.1876
        double a0=0.1184;
        double mz = 91.1876;
        double b5 = 11.0 - 2.0/3.0*5.0; // at Z mass all 5 flavors are active
        double b4 = 11.0 - 2.0/3.0*4.0;
        double b3 = 11.0 - 2.0/3.0*3.0;
        double lambda5 = mz * std::exp(-2.0*M_PI / (a0 * b5) );
        double lambda4 = std::pow( heavyqmasses[1], 1.0 - b5/b4) * std::pow(lambda5, b5/b4);
        double lambda3 = std::pow( heavyqmasses[0], 1.0 - b4/b3) * std::pow(lambda4, b4/b3);

        if (nf==5) lqcd=lambda5;
        else if (nf==4) lqcd=lambda4;
        else if (nf==3) lqcd=lambda3;
        else
            cerr << "WTF, nf=" << nf <<" at " << LINEINFO << endl;
        
    
    }

	
	if (alphas_freeze_c < 0.00001)  // sharp cutoff at maxalphas
	{
	    double scalefactor=0;
		if (std::abs(scaling-1.0)>0.0001)   // Don't use stored value
			scalefactor = scaling;
		else
			scalefactor = 4.0*Csqr;
		
		if (scalefactor/(rsqr*lqcd*lqcd) < 1.0) return maxalphas;
		double alpha = 4.0*M_PI/(  b0 * std::log(scalefactor/ (rsqr*lqcd*lqcd) ) );
		if (alpha>maxalphas)
			return maxalphas;
		return alpha;
	}

	// Smooth cutoff
    ///TODO: varying nf
    if (alphas_flavours==HEAVYQ)
        cerr << "Smooth cutoff alphas does not support (yet) heavy quarks " << LINEINFO << endl;
	return 4.0*M_PI / ( b0 * std::log(
		std::pow( std::pow(alphas_mu0, 2.0/alphas_freeze_c) + std::pow(4.0*Csqr/(rsqr*lambdaqcd*lambdaqcd), 1.0/alphas_freeze_c), alphas_freeze_c)	
		) );
		

}

std::string AmplitudeR::Alpha_s_str()
{

	std::stringstream ss;
	if (alphas_freeze_c < 0.0001)
		ss << "\\alpha_s = 12 \\pi / [ (33.0 - 2.0*Nf) * log(4.0*C^2/(r^2*lambdaqcd^2) ], C^2=" << Csqr << ", lambdaqcd=" << lambdaqcd << " GeV, maxalphas=" << maxalphas;
	else
		ss << "\\alpha_s = 12 \\pi / [ (33.0 - 2.0*Nf) * log[ ( mu0^(1/c) + (4C^2/(lambdaqcd^2*r^2))^(1/c) )^c) ] ], C^2=" << Csqr <<", lambdaqcd=" << lambdaqcd << " GeV, mu0=" << alphas_mu0;

    if (alphas_flavours == LIGHTQ)
        ss << " nf=" << Nf;
    else
        ss << " heavy quarks included, masses: m_c = " << heavyqmasses[0] << ", m_b = " << heavyqmasses[1];

    return ss.str();
}


int AmplitudeR::RPoints()
{
    return RPOINTS;
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
    double max = MAXR;

    return std::pow(max/MinR(), 1.0/(RPoints()-1));
}

REAL AmplitudeR::MaxR()
{
	// If we have already initialized everything, we can take 
	// maxr from the table
	
	if (rvals.size()>0)
		return rvals[rvals.size()-1];
    double max = MinR()*std::pow(RMultiplier(), RPoints()-1);

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
	if (bdep)
		cerr << "WTF, IMPACT PARAMETER DEPENDECE?????" << endl;
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

void AmplitudeR::SetAlphasFlavours(AlphasFlavours f)
{
    alphas_flavours = f;
}

AlphasFlavours AmplitudeR::GetAlphasFlavours()
{
    return alphas_flavours;
}
