
/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012
 */

#ifndef _BK_IC_DATAFILE_
#define _BK_IC_DATAFILE_

#include <string>
#include <string>
#include "interpolation.hpp"
#include "ic.hpp"
#include "config.hpp"

/*
 * Initial condition, where data is read from a given file
 * File syntax is: r x
 * Lines starting with one # are comments
 * Lines starting with three ### include parameters,
 * supported parameters are:
 * ###alphas_scaling:value   (Notice: not ln)
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
        double GetAlphasScaling();
	private:
		Interpolator *interpolator;
        double alphas_scaling;
};

#endif
