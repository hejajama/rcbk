/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2013
 */

#ifndef _BKCONFIG_HPP
#define _BKCONFIG_HPP

#include <string>

typedef double REAL;
const int Nc=3;
const int Nf=3;
const double LAMBDAQCD=0.241;



#ifndef LINEINFO
    #define LINEINFO __FILE__ << ":" << __LINE__
#endif




double StrToReal(std::string str);

int StrToInt(std::string str);



#endif

