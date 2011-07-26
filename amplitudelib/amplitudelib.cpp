/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "amplitudelib.hpp"
#include "datafile.hpp"
#include "../src/tools.hpp"
#include "../src/config.hpp"
#include <string>
#include <cmath>


/*
 * Calculate amplitude interpolating data
 * Interpolate rapidity linearly and r using spline
 */
REAL AmplitudeLib::N(REAL r, REAL y)
{
    if (r < MinR() or r > MaxR() )
    {
        cerr << "r must be between limits [" << MinR() << ", " << MaxR() << "]"
            << endl;
        return 0;
    }

    if (y<0 or y>yvals[yvals.size()-1] )
    {
        cerr << "r must be between limits [" << 0 << ", "
            << yvals[yvals.size()-1] << "]" << endl;
        return 0;
    }

    REAL lnr = std::log(r);
    int yind = FindIndex(y, yvals);
    int rind = FindIndex(lnr, lnrvals);

    int interpolation_points = INTERPOLATION_POINTS;

    int interpolation_start, interpolation_end;
    if (rind - interpolation_points/2 < 0)
    {
		interpolation_start=0;
		interpolation_end=interpolation_points;
	}
	else if (rind + interpolation_points/2 > RPoints()-1 )
	{
		interpolation_end = RPoints()-1;
		interpolation_start = RPoints()-interpolation_points-3;
	}
	else
	{
		interpolation_start = rind - interpolation_points/2;
		interpolation_end = rind + interpolation_points/2;
	}
    int interpo_points = interpolation_end - interpolation_start+1;
    REAL *tmparray = new REAL[interpo_points];
    REAL *tmpxarray = new REAL[interpo_points];
    for (int i=interpolation_start; i<= interpolation_end; i++)
    {
		tmpxarray[i-interpolation_start]=lnrvals[i];

        tmparray[i-interpolation_start] = n[yind][i];

        // Interpolate in y if possible
		if (yind < yvals.size()-1 )
        {
            tmparray[i-interpolation_start]=n[yind][i] 
                + (y - yvals[yind]) * (n[yind+1][i] - n[yind][i])
                / (yvals[yind+1]-yvals[yind]);
		} 
    }
    
    Interpolator interp(tmpxarray, tmparray, interpo_points);
    interp.Initialize();
    REAL result = interp.Evaluate(lnr);
    delete[] tmparray;
    delete[] tmpxarray;

    return result;

}

AmplitudeLib::AmplitudeLib(std::string datafile)
{
    DataFile data(datafile);
    data.GetData(n, yvals);
    minr = data.MinR();
    rmultiplier = data.RMultiplier();
    rpoints = data.RPoints();

    for (int i=0; i<rpoints; i++)
    {
        lnrvals.push_back(std::log(minr*std::pow(rmultiplier, i)));
    }

    cout << "# Data read from file " << datafile << ", minr: " << minr
        << " maxr: " << MaxR() << " rpoints: " << rpoints << " maxy "
        << yvals[yvals.size()-1] << endl;
}

int AmplitudeLib::RPoints()
{
    return rpoints;
}
REAL AmplitudeLib::MinR()
{
    return minr;
}

REAL AmplitudeLib::MaxR()
{
    return minr*std::pow(rmultiplier, rpoints-1);
}

int AmplitudeLib::YVals()
{
    return yvals.size();
}
