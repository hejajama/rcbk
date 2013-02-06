/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "solver.hpp"
#include <tools/tools.hpp>
#include <tools/interpolation.hpp>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h> // odeiv2 Requires GSL 1.15
#include <gsl/gsl_sf_bessel.h>
#include <cmath>


/// Accuracy parameters
const int THETAINTPOINTS = 200;
const double THETAINTACCURACY = 0.005;
const int RINTPOINTS = 200;
const REAL RINTACCURACY = 0.005;
const double DESOLVEACCURACY = 0.005; // orig 0.01, relative accuracy of de solver

Solver::Solver(AmplitudeR* _N)
{
    N=_N;
    deltay=0.2;
    bfkl=false;
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
        
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45; //2; //f45;

    gsl_odeiv_step * s    = gsl_odeiv_step_alloc (T, vecsize);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (0.0, DESOLVEACCURACY);    //abserr relerr
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
    } while (y <= maxy);

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);

    

    delete[] vec;


}

int EvolveR(REAL y, const REAL amplitude[], REAL dydt[], void *params)
{
    //cout << "Evolving with y=" <<y << endl;
    EvolutionHelperR* par = (EvolutionHelperR*)params;
    //int vecsize = par->N->BPoints()*par->N->RPoints()*par->N->ThetaPoints();

    // Intialize interpolator (works only if there is no b-dep.)
    REAL *tmprarray = new REAL[par->N->RPoints()];
    REAL *tmpyarray = new REAL[par->N->RPoints()];
    for (int i=0; i<par->N->RPoints(); i++)
    {
        tmprarray[i] = par->N->RVal(i);
        tmpyarray[i] = amplitude[i];
    }
    Interpolator interp(tmprarray, tmpyarray, par->N->RPoints());
    interp.Initialize();
    interp.SetFreeze(true);
    interp.SetUnderflow(0);
    interp.SetOverflow(1.0);
    
    #pragma omp parallel for schedule(dynamic, 15) // firstprivate(interp) 
    for (int rind=0; rind < par->N->RPoints(); rind++)
    {
        for (int thetaind=0; thetaind < par->N->ThetaPoints(); thetaind++)
        {
            for (int bind=0; bind < par->N->BPoints(); bind++)
            {
                int tmpind = rind*par->N->BPoints()*par->N->ThetaPoints()
                    + bind*par->N->ThetaPoints()+thetaind;

                ////THIS IS NOT TRUE IF THERE IS b-DEPENDENCE
                // optimize: as we know that the amplitude saturates to N=1,
                // we don't have to evolve it at large r
                if (rind>10 and !par->S->GetBfkl())
                {
                    if (amplitude[rind-2]>0.99999 and amplitude[rind-1]>0.99999
                        and amplitude[rind]>0.99999)
                    {
                        dydt[tmpind]=0;
                        /*#pragma omp critical
                        {
                            cout << "Skipping r=" << par->N->RVal(rind) << endl;
                        }*/
                        continue;
                    }
                }
                REAL tmplnr = par->N->LogRVal(rind);
                REAL tmplnb = par->N->LogBVal(bind);
                REAL tmptheta = par->N->ThetaVal(thetaind);
                /*cout << "tmpind " << tmpind << " maxrind "
                 << par->N->RPoints()-1 << " r=" << par->N->RVal(tmpind) <<
                 " amplitude " << par->S->InterpolateN(tmplnr, 0, 0, amplitude) ;*/
                dydt[tmpind] = par->S->RapidityDerivative(y, tmplnr, tmplnb, tmptheta,
                    amplitude, &interp);
                //cout << " reldydt " << dydt[tmpind]/par->S->InterpolateN(tmplnr, 0, 0, amplitude) << endl;
                
            }
        }
    }

    delete[] tmprarray;
    delete[] tmpyarray;
    
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
    REAL n02;
    REAL thetab; // angle between r01 and b01
    bool bdep;     // take into account impact parameter dependency
    REAL alphas_r01;
    REAL alphas_r02;
    const REAL* data;
    Interpolator* interp;
};
REAL Inthelperf_rint(REAL lnr, void* p);
REAL Inthelperf_thetaint(REAL theta, void* p);

REAL Solver::RapidityDerivative(REAL y,
            REAL lnr01, REAL lnb01, REAL thetab, const REAL* data,
            Interpolator *interp)
{

    if (lnr01 <= N->LogRVal(0)) lnr01*=0.999;
    else if (lnr01 >= N->LogRVal(N->RPoints()-1)) lnr01*=0.999;
    
    // Integrate first over r, then over \theta
    Inthelper_rthetaint helper;
    helper.N=N; helper.Solv=this;
    helper.y=y; helper.lnb01=lnb01; helper.lnr01=lnr01;
    helper.thetab = thetab; helper.data=data;
    helper.interp = interp;
    helper.bdep = N->ImpactParameter();
    //helper.n01 = InterpolateN(lnr01, lnb01, thetab, data);
    helper.n01 = interp->Evaluate(std::exp(lnr01));

    if (rc!=CONSTANT)
        helper.alphas_r01 = N->Alpha_s_ic(std::exp(2.0*lnr01));
    else
        helper.alphas_r01=0;

    gsl_function fun;
    fun.params = &helper;
    fun.function = Inthelperf_rint;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(RINTPOINTS);

    REAL minlnr = std::log( 1.000001*N->RVal(0) );
    REAL maxlnr = std::log( 0.99999*N->RVal(N->RPoints()-1) );

    int status; REAL result, abserr;
    status=gsl_integration_qag(&fun, minlnr,
            maxlnr, 0, RINTACCURACY, RINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        #pragma omp critical
        cerr << "RInt failed at " << LINEINFO << ": lnr=" << lnr01 <<", r= " << std::exp(lnr01) <<", thetab="
        << thetab <<", lnb01=" << lnb01 <<", result " << result << ", relerr "
        << std::abs(abserr/result) << endl;
    }

    return result;
}

REAL Inthelperf_rint(REAL lnr, void* p)
{
    Inthelper_rthetaint* par = (Inthelper_rthetaint*)p;


    par->lnr02=lnr;
    //par->n02 = par->Solv->InterpolateN(lnr, 0, 0, par->data);
    par->n02 = par->interp->Evaluate(std::exp(lnr));
    gsl_function fun;
    fun.function = Inthelperf_thetaint;
    fun.params = par;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);

    REAL mintheta = 0.0000001;
    REAL maxtheta = 2.0*M_PI-0.0001;

    if (!par->N->ImpactParameter()) maxtheta=M_PI-0.000001;
    if (par->Solv->GetRunningCoupling()!=CONSTANT)
         par->alphas_r02 = par->N->Alpha_s_ic(std::exp(2.0*lnr));
    else
        par->alphas_r02=0;
        
        //par->alphas_r02 = Alpha_s_r(std::exp(2.0*lnr),
         //   par->Solv->GetAlphasScaling());


    int status; REAL result, abserr;
    status = gsl_integration_qag(&fun, mintheta,
            maxtheta, 0, THETAINTACCURACY, THETAINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    if (status and std::abs(result)>1e-3) 
    {
        if (! (std::exp(par->lnr01) > 0.1 and std::exp(lnr) < 1e-5) )
            #pragma omp critical
            cerr << "Thetaint failed at " << LINEINFO <<": r=" << std::exp(lnr) 
                << " n01 " << par->n01 << " r01 " << std::exp(par->lnr01) <<
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

        REAL n12=0;
        if (r12sqr < SQR(par->N->MinR())) n12=0;
        else if (r12sqr > SQR(par->N->MaxR())) n12=1;
        else n12 = par->interp->Evaluate(std::sqrt(r12sqr));
        //else n12=par->Solv->InterpolateN(0.5*std::log(r12sqr), 0, 0, par->data);
        REAL alphas_r12=0;
        if (par->Solv->GetRunningCoupling()!=CONSTANT
            and par->Solv->GetRunningCoupling() != PARENT)
                alphas_r12 = par->N->Alpha_s_ic(r12sqr);
        
        REAL n02 = par->n02;
        //REAL n12 = par->Solv->InterpolateN(0.5*std::log(r12sqr), 0, 0, par->data);
        REAL n01 = par->n01;

        REAL result = n02 + n12 - n01;
        if (!par->Solv->GetBfkl()) result -= n02*n12;

        result *= par->Solv->Kernel(r01, r02, std::sqrt(r12sqr), par->alphas_r01,
            par->alphas_r02, alphas_r12, par->y, theta);

        return result;
    }

	cerr << "Impact parameter dependence is not supported!" << endl;
	exit(1);

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
REAL Solver::Kernel(REAL r01, REAL r02, REAL r12, REAL alphas_r01,
        REAL alphas_r02, REAL alphas_r12, REAL y,
        REAL theta2, REAL b01, REAL thetab)
{
    if (r12<N->MinR()/10 or r02 < N->MinR()/10)
        return 0;
    REAL result=0;

    // Ref for different prescriptions: 0704.0612
    // Convention: r01=r, r02 = r1, r12=r2
    REAL Rsqr, z, zsqrt;
    REAL costheta2; REAL r02dotr12;

    switch(rc)
    {
        case CONSTANT:
            result = ALPHABAR_s/(2.0*M_PI)
                    * SQR(r01) / ( SQR(r12) * SQR(r02) );
            break;
        case PARENT:
            result = alphas_r01*Nc/(2.0*SQR(M_PI))
                    * SQR(r01) / ( SQR(r12) * SQR(r02) );
            break;
        case BALITSKY:
            result = Nc/(2.0*SQR(M_PI))*alphas_r01
            * (
            SQR(r01) / ( SQR(r12) * SQR(r02) )
            + 1.0/SQR(r02)*(alphas_r02/alphas_r12 - 1.0)
            + 1.0/SQR(r12)*(alphas_r12/alphas_r02 - 1.0)
            );
			return result;
            break;
        case JIMWLK_SQRTALPHA:  // Similar than BALITSKY, slightly simpler (TL&HM JIMWLK paper)
            result = Nc/(2.0*SQR(M_PI))
            * (
            alphas_r02/SQR(r02) + alphas_r12/SQR(r12)
				+ std::sqrt(alphas_r02*alphas_r12)* ( SQR(r01/(r02*r12)) - 1.0/SQR(r02) - 1.0/SQR(r12)) 	// = -2 (x-z).(y-z)/(x-z)^2(y-z)^2
            );
            break;
        case KW:
            if (std::abs(SQR(r02)-SQR(r12)) < SQR(N->MinR())/1000) return 0;

            costheta2 = std::cos(theta2);
            r02dotr12 = SQR(r02) - r02*r01*costheta2;

            /* calculate dot product using geometry, slower than vector
             * calculateion used above
             if (theta2>M_PI) theta2=2.0*M_PI-theta2;

            theta012 = std::acos(
                -0.9999*( SQR(r02) - SQR(r01) - SQR(r12) ) / (2.0*r01*r12) );

            cosr1r2 = -std::cos(theta2 + theta012);
            r02dotr12 = r02*r12*cosr1r2;   // r02 \cdot r12
            */
            

            Rsqr = r02*r12*std::pow((r12/r02),
                (SQR(r02)+SQR(r12))/(SQR(r02)-SQR(r12))
                    - 2.0*SQR(r02)*SQR(r12)/(r02dotr12*(SQR(r02)-SQR(r12)) )
                );
            /*if (std::abs(r02-r12) < 1e-3 and std::abs(r02-r12)>1e-5)
                cout << "r1 " << r02 << " r2 " << r12 << " theta " << theta2 << " Rsqr " <<
                Rsqr << " exp " << r02*r12*std::exp(-1.0+1.0/costheta2) << endl;*/
            if (Rsqr<1e-50) return 0;

            result = Nc/(2.0*SQR(M_PI)) * (
                alphas_r02 / SQR(r02)
                - 2.0*alphas_r02*alphas_r12 / N->Alpha_s_ic(Rsqr)
                 * r02dotr12 / (SQR(r02)*SQR(r12) )
                + alphas_r12 / SQR(r12)
                );
            break;
        case MS:
            // Motyka & Staśto, 0901.4949: kinematical constraint, bessel kernel
            ///FIXME: 0.01 factor????
            cerr << "Motyka&Stasto kin. constraint kernel is most likely wrong..."
				<<" Check before you use this! " << LINEINFO << endl;
				
            z = 0.01*std::exp(-y); zsqrt = std::exp(-0.5*y);

            // r02 dot r12
            costheta2 = std::cos(theta2);
            r02dotr12 = SQR(r02) - r02*r01*costheta2;

            REAL besselr02r01 = gsl_sf_bessel_K1(r02/r01*zsqrt);
            REAL besselr12r01 = gsl_sf_bessel_K1(r12/r01*zsqrt);

            result = z/SQR(r01) * (
                  SQR( besselr02r01 )
                + SQR( besselr12r01 )
                -  2.0*besselr02r01*besselr12r01
                  * r02dotr12 / (r02*r12)
                );

            // Parent dipole RC
            result *= Nc/(2.0*SQR(M_PI))*alphas_r01;
            break;
            
    }
    return result;
        
}

/*
 * Interpolate data array
 */
REAL Solver::InterpolateN(REAL lnr, REAL lnb, REAL thetab, const REAL* data)
{
    cerr << "Solver::InterpolateN called, why?" << endl;
    // array format: ind = rind*BPoitns()*ThetaPoints()+bind*ThetaPoints()+Theta
    if (lnr > N->MaxLnR()) lnr = N->MaxLnR()*0.999;
    else if (lnr < N->MinLnR()) lnr = N->MinLnR()*0.999; 
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
        REAL *tmprarray = new REAL[points];
        for (int i=interpolation_start; i<=interpolation_end; i++)
        {
            tmpnarray[i-interpolation_start] = data[i];
            tmprarray[i-interpolation_start] = N->RVal(i);
        }
        Interpolator interp(tmprarray, tmpnarray, points);
        interp.Initialize();
        REAL result = interp.Evaluate(std::exp(lnr));

        if (result < -1e-3 or result>1.001)
        {
            cerr << "Interpolation result is " << result << ", lnr = "
            << lnr << ". " << LINEINFO << endl;
        }

        delete[] tmpnarray;
        delete[] tmprarray;

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

void Solver::SetDeltaY(REAL dy)
{
    deltay = dy;
}

void Solver::SetBfkl(bool bfkl_)
{
    bfkl=bfkl_;
}

bool Solver::GetBfkl()
{
    return bfkl;
}
