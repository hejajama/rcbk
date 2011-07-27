/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "solver.hpp"
#include "tools.hpp"
#include "interpolation.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h> // odeiv2 Requires GSL 1.15
#include <cmath>

Solver::Solver(AmplitudeR* _N)
{
    N=_N;
    deltay=0.15;
}

/*
 * Solve the BK
 */

struct EvolutionHelperR
{
    AmplitudeR* N;
    Solver* S;
};

int EvolveR(REAL y, const REAL amplitude[], REAL dydt[], void *params);

void Solver::Solve(REAL maxy)
{
    /*
     * In order to be able to use Runge Kutta we create one large array
     * which is basically just vector that we evolve
     * Array size is RPoints()*BPoints()*ThetaPoints(), and indexes are
     * ind = rind*BPoitns()*ThetaPoints()+bind*ThetaPoints()+Theta
     *
     * Values are ln of amplitude
     */

    int vecsize = N->BPoints()*N->RPoints()*N->ThetaPoints();
    REAL *vec = new REAL[vecsize];
    for (int rind=0; rind<N->RPoints(); rind++)
    {
        for (int bind=0; bind<N->BPoints(); bind++)
        {
            for (int thetaind=0; thetaind < N->ThetaPoints(); thetaind++)
            {
                 vec[rind*N->BPoints()*N->ThetaPoints() + bind*N->ThetaPoints()+thetaind]
                    = N->Ntable(0, rind, bind, thetaind);
            }
        }
    }

    REAL y=0; REAL step = deltay;  // We have always solved up to y

    // Intialize GSL
    EvolutionHelperR help; help.N=N; help.S=this;
    gsl_odeiv_system sys = {EvolveR, NULL, vecsize, &help};
        
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rk2;

    gsl_odeiv_step * s    = gsl_odeiv_step_alloc (T, vecsize);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (0.0, 0.05);    //abserr relerr
    gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (vecsize);
    REAL h = step;  // Initial ODE solver step size
    
    do
    {
        REAL nexty = y+step;
        while (y<nexty)
        {
            int status = gsl_odeiv_evolve_apply(e, c, s, &sys,
                &y, nexty, &h, vec);
            if (status != GSL_SUCCESS) {
                cerr << "Error in gsl_odeiv_evolve_apply at " << LINEINFO
                << ": " << gsl_strerror(status) << " (" << status << ")"
                << " y=" << y << ", h=" << h << endl;
            }
            cout << "Evolved up to " << y << "/" << nexty << ", h=" << h << endl;
        }

        
        int yind = N->AddRapidity(nexty);
        for (int rind=0; rind<N->RPoints(); rind++)
        {
            for (int bind=0; bind<N->BPoints(); bind++)
            {
                for (int thetaind=0; thetaind < N->ThetaPoints(); thetaind++)
                {
                    int tmpind = rind*N->BPoints()*N->ThetaPoints()
                            + bind*N->ThetaPoints()+thetaind;
                    N->AddDataPoint(yind, rind, bind, thetaind, vec[tmpind]);
                }
            }
        }
    } while (y < maxy);

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);

    

    delete[] vec;


}

int EvolveR(REAL y, const REAL amplitude[], REAL dydt[], void *params)
{
    cout << "Evolving with y=" <<y << endl;
    EvolutionHelperR* par = (EvolutionHelperR*)params;
    //int vecsize = par->N->BPoints()*par->N->RPoints()*par->N->ThetaPoints();
    #pragma omp parallel for
    for (int rind=0; rind < par->N->RPoints(); rind++)
    {
        for (int thetaind=0; thetaind < par->N->ThetaPoints(); thetaind++)
        {
            for (int bind=0; bind < par->N->BPoints(); bind++)
            {
                int tmpind = rind*par->N->BPoints()*par->N->ThetaPoints()
                    + bind*par->N->ThetaPoints()+thetaind;
                REAL tmplnr = par->N->LogRVal(rind);
                REAL tmplnb = par->N->LogBVal(bind);
                REAL tmptheta = par->N->ThetaVal(thetaind);
                dydt[tmpind] = par->S->RapidityDerivative(y, tmplnr, tmplnb, tmptheta,
                    amplitude);
                /*if (tmpind % 1 == 0) cout << "tmpind " << tmpind << " maxrind "
                 << par->N->RPoints()-1 << " r=" << par->N->RVal(tmpind) <<
                 " amplitude " << par->S->InterpolateN(tmplnr, 0, 0, amplitude) <<
                 " dydt " << dydt[tmpind] << endl;*/
            }
        }
    }
    return GSL_SUCCESS;
    
}

/*
 * Calculate \partial_Y N
 * Data contais the previously calculated data points at rapidity y, format
 * ind = rind*BPoitns()*ThetaPoints()+bind*ThetaPoints()+Theta
 */
struct Inthelper_rthetaint
{
    AmplitudeR* N;
    Solver* Solv;
    REAL y;
    REAL lnb01;   // impact parameter of the parent dipole
    REAL lnr01;   // parent dipole
    REAL lnr02;   // integrated over
    REAL n01;
    REAL thetab; // angle between r01 and b01
    bool bdep;     // take into account impact parameter dependency
    const REAL* data;
};
REAL Inthelperf_rint(REAL lnr, void* p);
REAL Inthelperf_thetaint(REAL theta, void* p);

REAL Solver::RapidityDerivative(REAL y,
            REAL lnr01, REAL lnb01, REAL thetab, const REAL* data)
{
    const int RINTPOINTS = 1000;
    const REAL RINTACCURACY = 0.01;
    
    // Integrate first over r, then over r, then \theta
    Inthelper_rthetaint helper;
    helper.N=N; helper.Solv=this;
    helper.y=y; helper.lnb01=lnb01; helper.lnr01=lnr01;
    helper.thetab = thetab; helper.data=data;
    helper.bdep = N->ImpactParameter();
    helper.n01 = InterpolateN(lnr01, lnb01, thetab, data);

    gsl_function fun;
    fun.params = &helper;
    fun.function = Inthelperf_rint;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(RINTPOINTS);

    REAL minr = 0.5*(N->LogRVal(0) + N->LogRVal(1));
    REAL maxr = 0.5*(N->LogRVal(N->RPoints()-1) + N->LogRVal(N->RPoints()-2));

    int status; REAL result, abserr;
    status=gsl_integration_qag(&fun, minr,
            maxr, 0, RINTACCURACY, RINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        cerr << "RInt failed at " << LINEINFO << ": lnr=" << lnr01 <<", thetab="
        << thetab <<", lnb01=" << lnb01 <<", result " << result << ", relerr "
        << std::abs(abserr/result) << endl;
    }

    
    if (result<0) cerr << "Amplitude seems to decrease from " << helper.n01
        << " at r01 = " << std::exp(lnr01)
        << ", y=" << y << " result " << result << " " << LINEINFO << endl;

    return result;
}

REAL Inthelperf_rint(REAL lnr, void* p)
{
    Inthelper_rthetaint* par = (Inthelper_rthetaint*)p;

    const int THETAINTPOINTS = 1000;
    const REAL THETAINTACCURACY = 0.05;

    par->lnr02=lnr;
    gsl_function fun;
    fun.function = Inthelperf_thetaint;
    fun.params = par;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);

     REAL mintheta = 0.0001;
     REAL maxtheta = 2.0*M_PI-0.0001;

     if (!par->N->ImpactParameter()) maxtheta=M_PI-0.0001;
     

    int status; REAL result, abserr;
    status = gsl_integration_qag(&fun, mintheta,
            maxtheta, 0, THETAINTACCURACY, THETAINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    if (status and std::abs(result)>1e-6)
    {
        cerr << "Thetaint failed at " << LINEINFO <<": r=" << std::exp(lnr) 
            << " n01 " << par->n01 <<
            " result " << result << " relerr " << std::abs(abserr/result) << endl;
    }

    result *= std::exp(2.0*lnr);        // Change of variables r->ln r

    if (!par->N->ImpactParameter()) result*=2.0;
        // As integration limits are [0,Pi]

    return result;

}

REAL Inthelperf_thetaint(REAL theta, void* p)
{
    Inthelper_rthetaint* par = (Inthelper_rthetaint*)p;

    // No impact parameter dependency, "easy" case
    if (!par->bdep)
    {
        REAL r01 = std::exp(par->lnr01);
        REAL r02 = std::exp(par->lnr02);
        REAL r12sqr = SQR(r01)+SQR(r02)-2.0*r01*r02*std::cos(theta);
        REAL n02 = par->Solv->InterpolateN(par->lnr02, 0, 0, par->data);
        REAL n12 = par->Solv->InterpolateN(0.5*std::log(r12sqr), 0, 0, par->data);
        REAL n01 = par->n01;

        REAL result = n02 + n12 - n01 - n02*n12;

        result *= par->Solv->Kernel(r01, r02, std::sqrt(r12sqr));

        return result;
    }

    ///TODO: use amplitude, not ln of it
    /*
    REAL r12sqr = std::exp(2.0*par->lnr01)+std::exp(2.0*par->lnr02)
        - 2.0*std::exp(par->lnr01+par->lnr02)*std::cos(theta);
    REAL b02sqr = std::exp(2.0*par->lnb01) + 0.25*std::exp(2.0*par->lnr01)
        + 0.25*std::exp(2.0*par->lnr02)
        - std::exp(par->lnb01+par->lnr01)*std::cos(par->thetab)
        + std::exp(par->lnb01 + par->lnr02)*std::cos(par->thetab - theta)
        - 0.5*std::exp(par->lnr01+par->lnr02)*std::cos(theta);
    REAL b12sqr = std::exp(2.0*par->lnb01)+0.25*std::exp(2.0*par->lnr02)
        + std::exp(par->lnb01+par->lnr02)*std::cos(par->thetab - theta);

    REAL thetab12 = std::asin( std::sqrt(std::exp(2.0*par->lnr02)/r12sqr
                                    *std::sin(theta)) ) - (M_PI - par->thetab);

    REAL  n02 = par->Solv->InterpolateN(par->lnr02, std::log(std::sqrt(b02sqr)),
            par->thetab - theta, par->data);
    REAL n12 = par->Solv->InterpolateN(std::log(std::sqrt(r12sqr)),
            std::log(std::sqrt(b12sqr)), thetab12, par->data );
    
    REAL result = 0;/*std::exp(n02) + std::exp(n12) - std::exp(par->ln_n01)
            - std::exp(n02+n12);*/ ///FIXME

    //result *= par->Solv->Kernel(std::exp(par->lnr01), std::exp(par->lnb01), par->thetab,
    //        std::exp(par->lnr02), theta, par->y);
    return 0;
}

/*
 * BK kernel
 * r01: parent dipole size, b01: parent dipole impact parameter,
 * thetab: angle between vecs r01 and b01, r02 integration variable, pos. of
 * point 2, theta2 is angle between r01 and r02, also integrated over
 * Rapidity y may be used with modified kernel
 * By default y=b01=tetab=theta=0
 *
 * Running coupling is also handled here
 */
REAL Solver::Kernel(REAL r01, REAL r02, REAL r12, REAL y,
        REAL b01, REAL thetab, REAL theta2)
{
    if (r12<1e-5 or r02 < 1e-5)
        return 0;
    REAL result;

    // Ref for different prescriptions: 0704.0612
    // Convention: r01=r, r02 = r1, r12=r2
    switch(rc)
    {
        case CONSTANT:
            result = ALPHAS*Nc/(2.0*SQR(M_PI))
                    * SQR(r01) / ( SQR(r12) * SQR(r02) );
            break;
        case PARENT:
            result = Alpha_s_r(SQR(r01), alphas_scaling)*Nc/(2.0*SQR(M_PI))
                    * SQR(r01) / ( SQR(r12) * SQR(r02) );
            break;
        case BALITSKY:
            result = Nc/(2.0*SQR(M_PI))*Alpha_s_r(SQR(r01), alphas_scaling)
            * (
            SQR(r01) / ( SQR(r12) * SQR(r02) )
            + 1.0/SQR(r02)*(Alpha_s_r(SQR(r02), alphas_scaling)
                        /Alpha_s_r(SQR(r12), alphas_scaling) - 1.0)
            + 1.0/SQR(r12)*(Alpha_s_r(SQR(r12), alphas_scaling)
                        /Alpha_s_r(SQR(r02), alphas_scaling) - 1.0)
            );

         /*   if (std::abs(result -  Alpha_s_r(SQR(r01), alphas_scaling)*Nc/(2.0*SQR(M_PI))
                    * SQR(r01) / ( SQR(r12) * SQR(r02) ) ) > 0.1)
                cout << " r01 " << r01 << " r02 " << r02 << " r12 " << r12
                 << " parent " << Alpha_s_r(SQR(r01), alphas_scaling)*Nc/(2.0*SQR(M_PI))
                    * SQR(r01) / ( SQR(r12) * SQR(r02) )
                    << " balitsky " << result << endl; */


            break;
        case KW:
            cerr << "KW RC is not implemented yet!" << endl;
            result=0;
            break;
    }
    return result;
        
}

/*
 * Interpolate data array
 */
REAL Solver::InterpolateN(REAL lnr, REAL lnb, REAL thetab, const REAL* data)
{
    // array format: ind = rind*BPoitns()*ThetaPoints()+bind*ThetaPoints()+Theta
    if (lnr > N->MaxLnR()) lnr = N->MaxLnR()*0.999;
    else if (lnr < N->MinLnR()) lnr = N->MinLnR()*0.999;  ///TODO: Same check for lnb
    int rind = FindIndex(lnr, N->LogRVals());
    if (rind<0) rind=0;
    int bind, thetaind;
    if (N->ImpactParameter())
    {
        bind = FindIndex(lnb, N->LogBVals());
        while (thetab<0) thetab += 2.0*M_PI;
        while (thetab > 2.0*M_PI) thetab-=2.0*M_PI;
        thetaind = FindIndex(thetab, N->ThetaVals());
        if (bind<0) bind=0;
        ///TODO: interpolate, exp
        return data[rind*N->BPoints()*N->ThetaPoints()+bind*N->ThetaPoints()+thetaind];
    }
    else
    {
        thetaind=0; bind=0;
        int size = N->RPoints();
        int interpolation_start, interpolation_end;
        if (rind - INTERPOLATION_POINTS/2 < 0)
        {
            interpolation_start=0; interpolation_end=INTERPOLATION_POINTS;
        }
        else if (rind + INTERPOLATION_POINTS/2 >= size)
        {
            interpolation_start = size-INTERPOLATION_POINTS/2-1;
            interpolation_end = size-1;
        }
        else
        {
            interpolation_start = rind - INTERPOLATION_POINTS/2;
            interpolation_end = rind + INTERPOLATION_POINTS/2;
        }
        int points = interpolation_end - interpolation_start + 1;
        REAL *tmpnarray = new REAL[points];
        REAL *tmplnrarray = new REAL[points];
        for (int i=interpolation_start; i<=interpolation_end; i++)
        {
            tmpnarray[i-interpolation_start] = data[i];
            tmplnrarray[i-interpolation_start] = N->LogRVal(i);
        }
        Interpolator interp(tmplnrarray, tmpnarray, points);
        interp.Initialize();
        REAL result = interp.Evaluate(lnr);

        if (result < 0 or result>1.1)
        {
            cerr << "Interpolation result is " << result << ", lnr = "
            << lnr << ". " << LINEINFO << endl;
        }

        delete[] tmpnarray;
        delete[] tmplnrarray;

        return result;
        
    }
}


void Solver::SetRunningCoupling(RunningCoupling rc_)
{
    rc=rc_;
}

RunningCoupling Solver::GetRunningCoupling()
{
    return rc;
}

void Solver::SetAlphasScaling(REAL scaling)
{
    alphas_scaling = scaling;
}

void Solver::SetDeltaY(REAL dy)
{
    deltay = dy;
}
