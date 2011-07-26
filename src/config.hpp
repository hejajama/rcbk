/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _CONFIG_HPP
#define _CONFIG_HPP 


#include <iostream>
using std::cout;
using std::cerr;
using std::cin;
using std::endl;
#include <string>
using std::string;
#include <cmath>

typedef double REAL;
typedef unsigned int uint;

// Physical constants
const REAL LAMBDAQCD2 = 0.21416*0.21416;   // GeV^2
const int Nf=3;
const int Nc=3;
const REAL ALPHA_e = 1.0/137.035999679; 
const REAL e = sqrt(4.0*M_PI*ALPHA_e);

// Reqularization of the running coupling
// Berger&Stasto 1010.0671: 0.3
// Fit to HERA: 0.7 0902.1112
const REAL MAXALPHA = 0.7;  
const REAL ALPHAS = 0.2;       // \alpha_s if RC=constant

// Other constants

const REAL eps=0.000001;

// Inline functions

#define LINEINFO __FILE__ << ":" << __LINE__

inline const REAL SQR(const REAL x) { return x*x; }


#endif
