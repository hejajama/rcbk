
/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012
 */

#ifndef _BK_IC_DATAFILE_
#define _BK_IC_DATAFILE_

#include <string>
#include <tools/config.hpp>
#include <string>
#include <tools/interpolation.hpp>
#include "ic.hpp"

/*
 * Initial condition, where data is read from a given file
 * File syntax is: r x
 */

class IC_datafile : public InitialCondition
{
	public:
		double DipoleAmplitude(double r, double b=0);
		int LoadFile(std::string file);
		IC_datafile();
		~IC_datafile();
		double MinR();
		double MaxR();
	private:
		Interpolator *interpolator;
};

#endif
