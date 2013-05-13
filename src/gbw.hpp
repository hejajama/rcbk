
/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012
 */

#ifndef _BK_GBW_IC_
#define _BK_GBW_IC_

/*
 * GBW and GBW^\gamma initial conditions
 */

#include <string>
#include "ic.hpp"
#include "amplitude.hpp"


class GBW : public InitialCondition
{
	public:
		GBW();		// Set MV1 parameters
		double DipoleAmplitude(double r, double b=0);
		void SetQsqr(double qsqr);
		void SetAnomalousDimension(double gamma_);
		void SetLambdaQcd(double lambda);
		std::string GetString();
	private:
		double qs0sqr;	// Q_{s0}^2 in GeV^2
		double anomalous_dimension;	// anomalous dimension
		double lambdaqcd;
};


#endif
