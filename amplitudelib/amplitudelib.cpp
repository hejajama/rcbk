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
#include <gsl/gsl_integration.h>
#include <algorithm>


extern "C"
{
    #include "fourier/fourier.h"
}

/*
 * Calculate amplitude interpolating data
 * Interpolate rapidity linearly and r using spline
 *
 * By default der=0, der=1 is 1st derivative w.r.t r, 2 2nd
 */
REAL AmplitudeLib::N(REAL r, REAL y, int der)
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

/*
 * Unintegrated gluon density
 * \psi(k, y) = C_F/( (2\pi)^3 \alphas(k) ) * \int d^2 r exp(-ik.r) \nabla^2 N_G,
 * \nabla^2 N_G = (1/r \der_r + \der^2_r)(2N-N^2)
 *  = 2/r \der N + 2 \der^2 N - 2N/r \der N - 2(\der r N)^2 - 2N \der^2 N
 *
 * For performance reasons it is probably a good idea to call
 * InitializeInterpoaltion before this
 */

struct UGDHelper
{
    REAL y;
    AmplitudeLib* N;
};
REAL UGDHelperf(REAL x, void* p);

REAL AmplitudeLib::UGD(REAL k, REAL y)
{
    if (k < UGD_IR_CUTOFF) return 0;
    
    set_fpu_state();
    init_workspace_fourier(1500);   // number of bessel zeroes, max 2000

    UGDHelper par;
    par.y=y; par.N=this; 
    REAL result = fourier_j0(k,UGDHelperf,&par);

    REAL Cf = (SQR(Nc)-1.0)/(2.0*Nc);
    result *= Cf/Alpha_s(SQR(k));  //Todo: scaling?
    result /= SQR(2.0*M_PI);

    return result;

}

REAL UGDHelperf(REAL r, void* p)
{
    UGDHelper* par = (UGDHelper*) p;
    REAL result=0;
    REAL dern=0, der2n=0, n=0;
    if (r > par->N->MaxR()) // N==1 for all r>MaxR()
    {
        n=1; dern=0; der2n=0;
    }
    else if (r < 1e-5)       //  TODO: par->N->MinR())
        return 0;
    else
    {       
        dern = par->N->N(r, par->y, 1);
        der2n = par->N->N(r, par->y, 2);
        n = par->N->N(r, par->y);
    }

    result = 2.0/r*dern + 2*der2n - 2.0*n/r*dern - 2.0*SQR(dern) - 2.0*n*der2n;
    return r*result;
}

/*
 * k_T factorization, calculate d\sigma^{A+B -> g} / (dy d^2 p_T)
 * Ref. 1011.5161 eq. (8)
 * Normalization is arbitrary
 */
struct Inthelper_dsigmadyd2pt
{
    AmplitudeLib* N;
    REAL pt, x1, x2;
    REAL kt;
    REAL sqrts,y;
};
REAL Inthelperf_dsigmadyd2pt(REAL kt, void* p);
REAL Inthelperf_dsigmadyd2pt_thetaint(REAL kt, void* p);
REAL AmplitudeLib::dSigmadyd2pt(REAL pt, REAL x1, REAL x2)
{
    const int KTINTPOINTS = 100;
    const REAL KTINTACCURACY = 0.05;
    
    Inthelper_dsigmadyd2pt helper;
    helper.N=this; helper.pt=pt; helper.x1=x1; helper.x2=x2;

    gsl_function fun;
    fun.params = &helper;
    fun.function = Inthelperf_dsigmadyd2pt;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(KTINTPOINTS);

    REAL minkt = 0;
    REAL maxkt = pt;

    int status; REAL result, abserr;
    status=gsl_integration_qag(&fun, minkt, maxkt,
            0, KTINTACCURACY, KTINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        cerr << "k_T integration failed at " << LINEINFO <<", pt=" << pt
            << ", x1=" << x1 <<", x2=" << x2 <<", result=" << result
            <<" relerr " << std::abs(abserr/result) << endl;
    }
    return result/(4.0*SQR(pt));    // 4.0 is in the integration measure d^2k/4
}

REAL Inthelperf_dsigmadyd2pt(REAL kt, void* p)
{
    const int THETAINTPOINTS = 80;
    const REAL THETAINTACCURACY = 0.05;
    
    Inthelper_dsigmadyd2pt* par = (Inthelper_dsigmadyd2pt*) p;
    par->kt=kt;
    gsl_function fun; fun.params=par;
    fun.function=Inthelperf_dsigmadyd2pt_thetaint;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);

    int status; REAL result, abserr;
    status=gsl_integration_qag(&fun, 0, M_PI,
            0, THETAINTACCURACY, THETAINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        cerr << "thetaintegration failed at " << LINEINFO <<", pt=" << par->pt
            << ", kt=" << kt <<", x1=" << par->x1 <<", x2=" << par->x2
            <<", result=" << result
            <<" relerr " << std::abs(abserr/result) << endl;
    }

    return 2.0*result;  // 2.0 from int. limits [0,2\pi] -> [0,\pi]
    
}

REAL Inthelperf_dsigmadyd2pt_thetaint(REAL theta, void* p)
{
    Inthelper_dsigmadyd2pt* par = (Inthelper_dsigmadyd2pt*) p;
    REAL costheta = std::cos(theta);
    REAL ktpluspt = SQR(par->kt) + SQR(par->pt) + 2.0*par->kt*par->pt*costheta;
    ktpluspt = std::sqrt(ktpluspt);
    REAL ktminuspt = SQR(par->kt) + SQR(par->pt) - 2.0*par->kt*par->pt*costheta;
    ktminuspt = std::sqrt(ktminuspt);
    
    REAL Q = std::max(ktpluspt/2.0, ktminuspt/2.0);
    REAL ugd1 = par->N->UGD(ktpluspt/2.0, par->x1);
    REAL ugd2 = par->N->UGD(ktminuspt/2.0, par->x2);

    if (isnan(ugd1) or isnan(ugd2) or Q<1e-10)
    {
        cerr << "???\n";
    }

    return Alpha_s(SQR(Q))*ugd1*ugd2;
}

/*
 * Rapidity distribution
 * d\sigma/dy = \int d^2 pt dSigmadydp2t
 */
REAL Inthelperf_dsigmady(REAL pt, void* p);
REAL AmplitudeLib::dSigmady(REAL y, REAL sqrts)
{
    REAL PTINTPOINTS = 400;
    REAL PTINTACCURACY = 0.05;
    
    REAL minpt = 0;
    REAL maxpt = 12;    // As in ref. 1011.5161

    Inthelper_dsigmadyd2pt helper;
    helper.N=this; helper.sqrts=sqrts; helper.y=y;

    gsl_function fun;
    fun.params = &helper;
    fun.function = Inthelperf_dsigmadyd2pt;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(PTINTPOINTS);

    int status; REAL result, abserr;
    status=gsl_integration_qag(&fun, minpt, maxpt,
            0, PTINTACCURACY, PTINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        cerr << "p_T integration failed at " << LINEINFO 
            << ", y=" << y <<", result=" << result
            <<" relerr " << std::abs(abserr/result) << endl;
    }

    return 0;
}

REAL Inthelperf_dsigmady(REAL pt, void* p)
{
    Inthelper_dsigmadyd2pt* par = (Inthelper_dsigmadyd2pt*)p;
    REAL x1 = pt/par->sqrts*std::exp(par->y);
    REAL x2 = pt/par->sqrts*std::exp(-par->y);

    return par->N->dSigmadyd2pt(pt, x1, x2);    

}

/*
 * Load data from a given file
 * Format is specified in file bk/README
 */
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
