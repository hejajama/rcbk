#include "ic_special.hpp"
#include <sstream>
#include <string>
#include <cmath>
#include <tools/config.hpp>

double IC_Special::DipoleAmplitude(double r, double b)
{
	double Qsqr = 0.2; double lambdaqcd = 0.241;
	const double beta = 0.01;	// quadratic action, MV1: beta=0
	double nr = Qsqr * r*r / 4.0 * std::log(1.0/(r*lambdaqcd ) + 2.7182)
		- beta * Qsqr * r*r * std::pow( std::log(1.0/(r*lambdaqcd) + 2.7182), 3);
	
	if (nr < 1e-5) return nr;
	else return 1.0 - std::exp(-nr);
	
}

std::string IC_Special::GetString()
{
	return "MV1 with quadratic action" ;
}
