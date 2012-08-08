#include "mv.hpp"
#include <sstream>
#include <string>
#include <cmath>
#include <tools/config.hpp>

double MV::DipoleAmplitude(double r, double b)
{
	if (b>1e-10)
		cerr << "Impact parameter is not supported!" << LINEINFO << endl;
	const double e = 2.7182818;
	///TODO: some algorithm to determina small r, e.g. when one has to linearize
    if (r < 2e-6)   ///NOTE: factor 1/4 "correctly", not as in AAMS paper
            return std::pow(SQR(r)*qs0sqr, anomalous_dimension)/4.0
            * std::log( 1.0/(r*lambdaqcd) + e);
    return 1.0 - std::exp(-std::pow(SQR(r)*qs0sqr, anomalous_dimension)/4.0
			* std::log( 1.0/(r*lambdaqcd) + e) );
}

void MV::SetQsqr(double qsqr)
{
	qs0sqr=qsqr;
}

void MV::SetAnomalousDimension(double gamma_)
{
	anomalous_dimension=gamma_;
}

void MV::SetLambdaQcd(double lambda)
{
	lambdaqcd=lambda;
}

std::string MV::GetString()
{
	std::stringstream ss;
	ss << "MV model, Q_s0^2 = " << qs0sqr << " GeV^2, \\gamma = " << anomalous_dimension
		<< ", x0=" << x0 << ", \\Lambda_QCD = " << lambdaqcd << " GeV";
	return ss.str();
}

/*
 * Set some reasonable parameters
 * That is, MV1 model for nucleus
 */
MV::MV()
{
	qs0sqr = 0.72;
	x0=0.007;
	lambdaqcd=LAMBDAQCD;
	anomalous_dimension=1;
}
