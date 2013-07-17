
/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012
 */

#ifndef _BK_MV_SPECIAL_
#define _BK_MV_SPECIAL_

/*
 * Special IC
 */

#include <string>
#include <tools/config.hpp>
#include "ic.hpp"
#include "amplitude.hpp"


/*
 * Special class which can be used when BK equation is solved using
 * some hardcoded IC
 */

class IC_Special : public InitialCondition
{
	public:
		double DipoleAmplitude(double r, double b=0);
		std::string GetString();
};



#endif
