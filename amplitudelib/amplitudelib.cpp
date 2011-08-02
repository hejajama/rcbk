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
#include <gsl/gsl_roots.h>


extern "C"
{
    #include "fourier/fourier.h"
}

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
    //REAL lnr = std::log(r);

    /// Use already initialized interpolator
    if (std::abs(y - interpolator_y) < 0.01)
    {
        REAL result=0;
        if (der==0)
            result = interpolator->Evaluate(r);
        if (der==1)
        {
            result = interpolator->Derivative(r);
            //result /= r;    // dN / d ln r = r dN/dr
        }
        if (der==2)
        {
            result = interpolator->Derivative2(r);
            //result /= SQR(r);   // d^2N / d lnr^2 = r^2 d^2 N / dr^2
        }
        return result;

    }


    /// Initialize new interpolator and use it
    int yind = FindIndex(y, yvals);
    int rind = FindIndex(r, rvals);

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
		tmpxarray[i-interpolation_start]=rvals[i];

        tmparray[i-interpolation_start] = n[yind][i];

        // Interpolate in y if possible
		if (yind < yvals.size()-1 )
        {
            tmparray[i-interpolation_start]=n[yind][i] 
                + (y - yvals[yind]) * (n[yind+1][i] - n[yind][i])
                / (yvals[yind+1]-yvals[yind]);
		}
        if (sqr)
            cerr << "sqr is not working!\n";
    }
    
    Interpolator interp(tmpxarray, tmparray, interpo_points);
    interp.Initialize();
    REAL result=0;
    if (der==0)
        result = interp.Evaluate(r);
    if (der==1)
    {
        result = interp.Derivative(r);
        //result /= r;    // dN / d ln r = r dN/dr
    }
    if (der==2)
    {
        result = interp.Derivative2(r);
        //result /= SQR(r);   // d^2N / d lnr^2 = r^2 d^2 N / dr^2
    }
    
    delete[] tmparray;
    delete[] tmpxarray;

    return result;

}

/*
 * FT the amplitude to the k-space
 * N(k) = \int d^2 r/(2\pi r^2) exp(ik.r)N(r)
 *  = \int dr/r BesselJ[0,k*r] * N(r)
 * 
 * Note: for performance reasons it is probably a good idea to
 * call AmplitudeLib::InitializeInterpolation(y) before this
 */
struct N_k_helper
{
    REAL y; REAL kt;
    AmplitudeLib* N;
};
REAL N_k_helperf(REAL r, void* p);
REAL AmplitudeLib::N_k(REAL kt, REAL y)
{
    // Some initialisation stuff -
    set_fpu_state();
    init_workspace_fourier(1500);   // number of bessel zeroes, max 2000
    
    N_k_helper par;
    par.y=y; par.N=this; par.kt=kt;
    REAL result = fourier_j0(kt,N_k_helperf,&par);
    return result;
}

REAL N_k_helperf(REAL r, void* p)
{
    N_k_helper* par = (N_k_helper*) p;
    if (r < par->N->MinR()) return 0;
    else if (r > par->N->MaxR()) return 1.0/r;
    return 1.0/r*par->N->N(r, par->y);
}

AmplitudeLib::AmplitudeLib(std::string datafile)
{
    DataFile data(datafile);
    data.GetData(n, yvals);
    minr = data.MinR();
    rmultiplier = data.RMultiplier();
    rpoints = data.RPoints();

    tmprarray = new REAL[data.RPoints()];
    tmpnarray = new REAL[data.RPoints()];
    for (int i=0; i<rpoints; i++)
    {
        lnrvals.push_back(std::log(minr*std::pow(rmultiplier, i)));
        rvals.push_back(std::exp(lnrvals[i]));
        tmprarray[i] = std::exp(lnrvals[i]);
    }

    cout << "# Data read from file " << datafile << ", minr: " << minr
        << " maxr: " << MaxR() << " rpoints: " << rpoints << " maxy "
        << yvals[yvals.size()-1] << endl;

    interpolator_y=-1;  // if >=0, interpolator is initialized, must free
    // memory (delete tmprarray and tmpnarray at the end)
}

/*
 * d ln N / d ln r^2 = 1/N * d N / d ln r^2 = 1/N d N / dr^2  dr^2 / d ln r^2
 *  = 1/N d N / dr^2 r^2 = 1/N d N / dr  dr / dr^2  r^2 = r/(2N) * dN/dr
 */
REAL AmplitudeLib::LogLogDerivative(REAL r, REAL y)
{
    REAL dndr = N(r,y,1);
    return r/(2.0*N(r,y))*dndr;
}

/*
 * Calculate saturation scale, definition is
 * N(r_s, y) = Ns = (for example) 0.5
 */
struct SatscaleSolverHelper
{
    REAL y;
    REAL Ns;
    AmplitudeLib* N;
};
REAL SatscaleSolverHelperf(double r, void* p)
{
    SatscaleSolverHelper* par = (SatscaleSolverHelper*)p;
    return par->N->N(r, par->y) - par->Ns;
}

REAL AmplitudeLib::SaturationScale(REAL y, REAL Ns)
{
    const int MAX_ITER = 300;
    const REAL ROOTFINDACCURACY = 0.01;

    SatscaleSolverHelper helper;
    helper.y=y; helper.Ns=Ns; helper.N=this;
    
    gsl_function f;
    f.params = &helper;
    f.function = &SatscaleSolverHelperf;

    const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    
    gsl_root_fsolver_set(s, &f, MinR()*1.0001, MaxR()*0.999);
    int iter=0; int status; double min,max;
    do
    {
        iter++;
        gsl_root_fsolver_iterate(s);
        min = gsl_root_fsolver_x_lower(s);
        max = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(min, max, 0, ROOTFINDACCURACY);    

    } while (status == GSL_CONTINUE and iter < MAX_ITER);

    if (iter>=MAX_ITER)
    {
        cerr << "Solving failed at y=" << y << endl;
    }

    REAL res = gsl_root_fsolver_root(s);

    gsl_root_fsolver_free(s);

    return res;
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
        REAL tmpr = tmprarray[i];
        if (i==0) tmpr*=1.001; if (i==rpoints-1) tmpr*=0.999;
        tmpnarray[i] = N(tmpr, y);
    }
    interpolator = new Interpolator(tmprarray, tmpnarray, rpoints);
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
    delete[] tmprarray;
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

REAL AmplitudeLib::MaxY()
{
    return yvals[yvals.size()-1];
}
