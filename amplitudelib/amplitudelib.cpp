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
 *
 * By default der=0, der=1 is 1st derivative w.r.t r, 2 2nd
 * if (sqr), then work with N^2, by default sqr=false
 */
REAL AmplitudeLib::N(REAL r, REAL y, int der, bool sqr)
{
    if (der>2 or der<0)
    {
        cerr << "Derivative degree " << der << " is not 0, 1 or 2! " << LINEINFO
         << endl;
         return 0;
    }
    if (r < MinR() or r > MaxR() )
    {
        cerr << "r must be between limits [" << MinR() << ", " << MaxR() << "]"
            << " asked r=" << r << " " << LINEINFO
            << endl;
        return 0;
    }

    if (y<0 or y>yvals[yvals.size()-1] )
    {
        cerr << "r must be between limits [" << 0 << ", "
            << yvals[yvals.size()-1] << "], asked y=" << y << " "
            << LINEINFO << endl;
        return 0;
    }
    REAL lnr = std::log(r);

    /// Use already initialized interpolator
    if (std::abs(y - interpolator_y) < 0.01)
    {
        REAL result;
        if (der==0)
            result = interpolator->Evaluate(lnr);
        if (der==1)
        {
            result = interpolator->Derivative(lnr);
            result /= r;    // dN / d ln r = r dN/dr
        }
        if (der==2)
        {
            result = interpolator->Derivative2(lnr);
            result /= SQR(r);   // d^2N / d lnr^2 = r^2 d^2 N / dr^2
        }
        return result;

    }


    /// Initialize new interpolator and use it
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
        if (sqr)
            tmparray[i-interpolation_start] = SQR(tmparray[i-interpolation_start]);
    }
    
    Interpolator interp(tmpxarray, tmparray, interpo_points);
    interp.Initialize();
    REAL result;
    if (der==0)
        result = interp.Evaluate(lnr);
    if (der==1)
    {
        result = interp.Derivative(lnr);
        result /= r;    // dN / d ln r = r dN/dr
    }
    if (der==2)
    {
        result = interp.Derivative2(lnr);
        result /= SQR(r);   // d^2N / d lnr^2 = r^2 d^2 N / dr^2
    }
    
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

    tmplnrarray = new REAL[data.RPoints()];
    tmpnarray = new REAL[data.RPoints()];
    for (int i=0; i<rpoints; i++)
    {
        lnrvals.push_back(std::log(minr*std::pow(rmultiplier, i)));
        tmplnrarray[i] = lnrvals[i];
    }

    cout << "# Data read from file " << datafile << ", minr: " << minr
        << " maxr: " << MaxR() << " rpoints: " << rpoints << " maxy "
        << yvals[yvals.size()-1] << endl;

    interpolator_y=-1;  // if >=0, interpolator is initialized, must free
    // memory (delete tmprarray and tmpnarray at the end)
}

/*
 * Initializes interpolation method with all read data points at given y
 */
void AmplitudeLib::InitializeInterpoaltion(REAL y)
{
    if (interpolator_y>=0)
    {
        delete interpolator;
    }
    for (int i=0; i<rpoints; i++)
    {
        REAL tmpr = std::exp(tmplnrarray[i]);
        if (i==0) tmpr*=1.001; if (i==rpoints-1) tmpr*=0.999;
        tmpnarray[i] = N(tmpr, y);
    }
    interpolator = new Interpolator(tmplnrarray, tmpnarray, rpoints);
    interpolator->Initialize();
    interpolator_y = y;
}

AmplitudeLib::~AmplitudeLib()
{
    if (interpolator_y>=0)
    {
        delete interpolator;
    }
    delete[] tmpnarray;
    delete[] tmplnrarray;
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
