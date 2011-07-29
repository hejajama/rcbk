/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include "interpolation.hpp"


/*
 * Intialize interpolation
 * Returns -1 in case of error, 0 otherwise
 */
int Interpolator::Initialize()
{
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            acc = gsl_interp_accel_alloc();
            spline = gsl_spline_alloc(gsl_interp_cspline, points);
            gsl_spline_init(spline, xdata, ydata, points);
            break;
        case INTERPOLATE_BSPLINE:
            gsl_vector *x = gsl_vector_alloc(points);
            gsl_vector *y = gsl_vector_alloc(points);
            gsl_vector *w = gsl_vector_alloc(points);

            for (int i=0; i< points; i++)
            {
                gsl_vector_set(x, i, xdata[i]);
                gsl_vector_set(y, i, ydata[i]);
                gsl_vector_set(w, i, 1.0);
            }
     
            /* allocate a cubic bspline workspace (k = 4) */
            bw = gsl_bspline_alloc(k, nbreak);
            derbw = gsl_bspline_deriv_alloc(k);
            B = gsl_vector_alloc(ncoeffs);
       
            X = gsl_matrix_alloc(points, ncoeffs);
            c = gsl_vector_alloc(ncoeffs);
       
            cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
            mw = gsl_multifit_linear_alloc(points, ncoeffs);
     
     
            // use uniform breakpoints
            gsl_bspline_knots_uniform(xdata[0], xdata[points-1], bw);
     
            /* construct the fit matrix X */
            for (int i = 0; i < points; ++i)
            {
               double xi = gsl_vector_get(x, i);
             
               /* compute B_j(xi) for all j */
               gsl_bspline_eval(xi, B, bw);
             
               /* fill in row i of X */
               for (int j = 0; j < ncoeffs; ++j)
               {
                  double Bj = gsl_vector_get(B, j);
                  gsl_matrix_set(X, i, j, Bj);
               }
            }
     
            /* do the fit */
            REAL chisq;
            gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);

            gsl_vector_free(x);
            gsl_vector_free(y);
            gsl_vector_free(w);

            break;
    }
    ready=true;
    return 0;   //ok, there is no error handling at the moment...
}


REAL Interpolator::Evaluate(REAL x)
{
    if (!ready)
    {
        cerr << "Interpolator is not ready! Did you forget to call Interpolator::Initialize()?" << endl;
        return 0;
    }
    REAL res, yerr; int status;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_e(spline, x, acc, &res);
            if (status)
            {
                cerr << "Interpolatioin failed at " << LINEINFO << ", error " << gsl_strerror(status)
                 << " (" << status << "), x=" << x << ", minx=" << xdata[0]
                 << ", maxx=" << xdata[points-1] << ", result=" << res << endl;
            }
            return res;
            break;
        case INTERPOLATE_BSPLINE:
            gsl_bspline_eval(x, B, bw);
            gsl_multifit_linear_est(B, c, cov, &res, &yerr);

            /*if (std::abs(yerr/res)>0.05 )
            {
                cerr << "Interpolation failed at " << LINEINFO << ": bspline result "
                << res << " pm " << yerr << " relerr " << std::abs(yerr/res) << endl;
            }*/
            return res;
            break;
    }

    cerr << "Interpolation method is invalid! " << LINEINFO << endl;
    return 0;   //Shoudn't end up here
}

REAL Interpolator::Derivative(REAL x)
{
    REAL res; int status;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_deriv_e(spline, x, acc, &res);
            break;
        case INTERPOLATE_BSPLINE:
            gsl_matrix* mat = gsl_matrix_alloc(nbreak+k-2, 2);
            gsl_bspline_deriv_eval(x, 1, mat, bw, derbw);
            for (int i=0; i<ncoeffs; i++)
            {
                res += gsl_vector_get(c, i)*gsl_matrix_get(mat, i, 1);
            }
            gsl_matrix_free(mat);
            return res;
            
    }

    return res;
}

REAL Interpolator::Derivative2(REAL x)
{
    REAL res; int status;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_deriv2_e(spline, x, acc, &res);
            break;
        case INTERPOLATE_BSPLINE:
            cerr << "2nd derivative is not implemented for BSPLINE interpolation!"
            << " " << LINEINFO << endl;
            break;
    }
    return res;

}

Interpolator::Interpolator(REAL *x, REAL *y, int p)
{
    points=p;
    xdata=x;
    ydata=y;
    method = INTERPOLATE_SPLINE;
    ready=false;
}

void Interpolator::SetMethod(INTERPOLATION_METHOD m)
{
    method = m;
}

void Interpolator::Clear()
{
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
            break;
        case INTERPOLATE_BSPLINE:
            gsl_bspline_free(bw);
            gsl_bspline_deriv_free(derbw);
            gsl_vector_free(B);
            gsl_matrix_free(X);
            gsl_vector_free(c);
            
            gsl_matrix_free(cov);
            gsl_multifit_linear_free(mw);
            break;
    }
}

Interpolator::~Interpolator()
{
    Clear();

}

REAL* Interpolator::GetXData()
{
    return xdata;
}
REAL* Interpolator::GetYData()
{
    return ydata;
}
int Interpolator::GetNumOfPoints()
{
    return points;
}
INTERPOLATION_METHOD Interpolator::GetMethod()
{
    return method;
}

// Copy data from given class and initialize this, as this is
// the copy constructor
Interpolator::Interpolator(Interpolator& inter)
{
    points=inter.GetNumOfPoints();
    xdata = inter.GetXData();
    ydata = inter.GetYData();
 
    method = inter.GetMethod();
    Initialize();
}
