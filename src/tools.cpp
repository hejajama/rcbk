/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "tools.hpp"
#include "config.hpp"
#include <string>
#include <sstream>
#include <cmath>
#include <vector>

/*
 * Str to REAL/int
 */
REAL StrToReal(std::string str)
{
    std::stringstream buff(str);
    REAL tmp;
    buff >> tmp;
    return tmp;
}

int StrToInt(std::string str)
{
    std::stringstream buff(str);
    int tmp;
    buff >> tmp;
    return tmp;
}

// GSL Error handler
int errors;
void ErrHandler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno)
{
    
    // Errors related to convergence of integrals are handled when
    // gsl_integration functions are called, don't do anything with them here
    // 14 = failed to reach tolerance
    // 18 = roundoff error prevents tolerance from being achieved
    // 11 = maximum number of subdivisions reached
    if (gsl_errno == 14 or gsl_errno == 18 or gsl_errno == 11)
        return;

    // 15 in expint.c:363: gsl_sf_gamma_inc underflows
    if (gsl_errno == 15 and string(file)=="expint.c" and line==363) return;

    errors++;
    std::cerr << file << ":"<< line <<": Error " << errors << ": " <<reason
            << " (code " << gsl_errno << ")." << std::endl;
}

/* 
 * Q^2 dependent strong coupling constant
 * Takes into account only u,d ands s quarks
 *
 * Reqularized to avoind infrared divergence following e.g.
 * Berger&Stasto 1010.0671 [hep-ph]
 * Reqularization is set in file config.hpp
 */
REAL Alpha_s(REAL Qsqr, REAL scaling)
{
    if (scaling*Qsqr < LAMBDAQCD2)
        return MAXALPHA;
    REAL alpha = 12.0*M_PI/( (33.0-2.0*Nf)*log(scaling*Qsqr/LAMBDAQCD2) );
    if (alpha > MAXALPHA)
        return MAXALPHA;
    return alpha;
}

REAL Alphabar_s(REAL Qsqr, REAL scaling)
{
    return Alpha_s(Qsqr, scaling)*Nc/M_PI;
}

std::string Alpha_s_str()
{
    std::stringstream stream;
    stream << "\\alpha_s is freezed at " << MAXALPHA;
    return stream.str();
}

/* Returns index i for which
 * vec[i]<=val
 * If such index can't be found, returns -1
 */

int FindIndex(REAL val, std::vector<REAL> &vec)
{
    if (val < vec[0]) return -1;
    
    int ind=-1;
    
    uint start=0; uint end=vec.size()-1;
    while(end-start>10)
    {
        int tmp = static_cast<int>((start+end)/2.0);
        
        if (vec[tmp]>=val)
            end=tmp;
        else
            start=tmp;
    }
    
    
    for (uint i=start; i<=end; i++)
    {
        if (vec[i]<=val and vec[i+1]>val)
        {
            ind=i;
            break;
        }
    }
    if (ind == -1) return vec.size()-1;
    return ind;
}
