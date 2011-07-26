/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "../src/tools.hpp"
#include "amplitudelib.hpp"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>

const std::string version = "v. 0.1  2011-xx-xx";

int main(int argc, char* argv[])
{
    cout << "# Reading data from file " << argv[1] << endl;
    AmplitudeLib N(argv[1]);

    REAL y = StrToReal(argv[2]);
    cout << "# y = " << y << endl;
    cout << "# r [1/GeV]     Amplitude" << endl;
    for (REAL r=N.MinR(); r<N.MaxR(); r*=1.1)
    {
        cout << r << " " << N.N(r, y) << endl;

    }
    

    return 0;
}
