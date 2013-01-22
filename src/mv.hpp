
/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012
 */

#ifndef _BK_MV_IC_
#define _BK_MV_IC_

/*
 * MV and MV^\gamma initial conditions
 */

#include <string>
#include <tools/config.hpp>
#include "ic.hpp"
#include "amplitude.hpp"


class MV : public InitialCondition
{
	public:
		MV();		// Set MV1 parameters
		double DipoleAmplitude(double r, double b=0);
		void SetQsqr(double qsqr);
		void SetAnomalousDimension(double gamma_);
		void SetLambdaQcd(double lambda);
		void SetE(double ec);	// coefficient c of e in Log[1/r\Lambda + cE]
		std::string GetString();
	private:
		double qs0sqr;	// Q_{s0}^2 in GeV^2
		double anomalous_dimension;	// anomalous dimension
		double lambdaqcd;
		double ec;	// coefficient c of e in Log[1/r\Lambda + cE]
};


#endif
