
/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012
 */

#ifndef _BK_IC_
#define _BK_IC_

#include <string>
#include <tools/config.hpp>

/*
 * Initial condition class, all other ic's are derived from this
 */
 
class InitialCondition
{
	public:
		virtual double DipoleAmplitude(double r, double b=0)=0;
		virtual std::string GetString();
		double X0();
		void SetX0(double x0_);
	protected:
		double x0;
};


#endif
