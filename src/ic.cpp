
/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */


/*
 * Initial condition class, all other ic's are derived from this
 */
 
#include "ic.hpp"
#include <string>

std::string InitialCondition::GetString()
{
	return "String is not impelmented for the current IC";
}

double InitialCondition::X0()
{
	return x0;
}

void InitialCondition::SetX0(double x0_)
{
	x0=x0_;
}

double InitialCondition::MinR()
{
	return 1e-99;
}

double InitialCondition::MaxR()
{
	return 1e99;
}
