#include "ic_special.hpp"
#include <sstream>
#include <string>
#include <cmath>


double IC_Special::DipoleAmplitude(double r, double b)
{
    /*
	double qs0sqr = 0.165; double anomalous_dimension = 1.135; 
	double lambdaqcd = 0.241;
	double sigma0 = 32.895 * 2.568; // in GeV
	double A = 208;
	InitializeWSDistribution(208);
	const double e = 2.7182818;
	
	
	double proton_n = std::pow(SQR(r)*qs0sqr, anomalous_dimension)/4.0
            * std::log( 1.0/(r*lambdaqcd) + e);
	
	double s = A/2.0 * T_A(11.57, 208) * sigma0 * proton_n;
	if (s<1e-5) return s;
	else return 1.0-std::exp(-s);
     */
	
	// Quadratic action
	/*double Qsqr = 0.168; double lambdaqcd = 0.241;
	const double beta = 0.01;	// quadratic action, MV1: beta=0
	double nr = Qsqr * r*r / 4.0 * std::log(1.0/(r*lambdaqcd ) + 2.7182)
		- 0*beta * Qsqr * r*r * std::pow( std::log(1.0/(r*lambdaqcd) + 2.7182), 3);
	
	if (nr < 1e-5) return nr;
	else return 1.0 - std::exp(-nr);
	*/
}

std::string IC_Special::GetString()
{
	return "MV1 with quadratic action" ;
}
