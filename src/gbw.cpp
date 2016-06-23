#include "gbw.hpp"
#include <sstream>
#include <string>
#include <cmath>
#include <tools/config.hpp>

using Amplitude::SQR;

double GBW::DipoleAmplitude(double r, double b)
{
	if (b>1e-10)
		cerr << "Impact parameter is not supported!" << LINEINFO << endl;
	///TODO: some algorithm to determina small r, e.g. when one has to linearize
    if (r < 2e-6) 
            return std::pow(SQR(r)*qs0sqr, anomalous_dimension)/4.0;
    return 1.0 - std::exp(-std::pow(SQR(r)*qs0sqr, anomalous_dimension)/4.0);
}

void GBW::SetQsqr(double qsqr)
{
	qs0sqr=qsqr;
}

void GBW::SetAnomalousDimension(double gamma_)
{
	anomalous_dimension=gamma_;
}

void GBW::SetLambdaQcd(double lambda)
{
	lambdaqcd=lambda;
}

std::string GBW::GetString()
{
	std::stringstream ss;
	ss << "GBW model, Q_s0^2 = " << qs0sqr << " GeV^2, \\gamma = " << anomalous_dimension
		<< ", x0=" << x0 << ", \\Lambda_QCD = " << lambdaqcd << " GeV";
	return ss.str();
}

/*
 * Set some reasonable parameters
 * That is, MV1 model for nucleus
 */
GBW::GBW()
{
	qs0sqr = 0.72;
	x0=0.007;
	lambdaqcd=LAMBDAQCD;
	anomalous_dimension=1;
}
