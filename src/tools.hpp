/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _TOOLS_H
#define _TOOLS_H

#include <string>
#include <vector>
#include "config.hpp"

REAL StrToReal(std::string str);
int StrToInt(std::string str);
// GSL error handler
void ErrHandler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno);

REAL Alpha_s(REAL Qsqr, REAL scaling=1.0);
REAL Alpha_s_r(REAL rsqr, REAL scaling=1.0);
REAL Alphabar_s(REAL Qsqr, REAL scaling=1.0); // \alpha_s N_C / Pi

int FindIndex(REAL val, std::vector<REAL> &array);

std::string Alpha_s_str();

#endif

