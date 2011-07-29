/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _INTERPOLATION_H
#define _INTERPOLATION_H

/*
 * Interpolation helper
 * Interpolates given data using spline (goes trough every data point)
 * or bspline (=noisy data)
 * Uses GSL
 */

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include "config.hpp"
#include "tools.hpp"

enum INTERPOLATION_METHOD {
    INTERPOLATE_SPLINE,
    INTERPOLATE_BSPLINE
};

class Interpolator
{
    public:
        Interpolator(REAL* x, REAL* y, int p);
        Interpolator(Interpolator& inter);
        ~Interpolator();
        void Clear();
        REAL Evaluate(REAL x);
        REAL Derivative(REAL x);    // 1st derivative
        REAL Derivative2(REAL x);   // 2nd derivative
        void SetMethod(INTERPOLATION_METHOD m);
        int Initialize();

        REAL* GetXData();
        REAL* GetYData();
        int GetNumOfPoints();
        INTERPOLATION_METHOD GetMethod();

    private:
        INTERPOLATION_METHOD method;
        REAL* xdata, *ydata;
        int points;
        bool ready;
        
        // spline
        gsl_interp_accel *acc;
        gsl_spline *spline;

        // bspline
        gsl_bspline_workspace *bw;
        gsl_bspline_deriv_workspace *derbw;
        gsl_vector *B;
        gsl_vector *c;
        gsl_matrix *X;
        gsl_matrix *cov;
        gsl_multifit_linear_workspace *mw;

        static const int k=4;
        static const int ncoeffs = 12;
        static const int nbreak = ncoeffs-k+2;


};




#endif
