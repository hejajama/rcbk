/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "../src/tools.hpp"
#include "ugd.hpp"
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
    gsl_set_error_handler(&ErrHandler);
    cout << "# Reading data from file " << argv[1] << endl;
    AmplitudeLib N(argv[1]);
    REAL y = StrToReal(argv[2]);
    N.InitializeInterpoaltion(y);
    /*UGD ugd(&N);

    for (REAL k=0.1; k<10; k*=1.1)
    {
        cout << k << " " << ugd.Evaluate(k, y) << endl;
    }
    return 0;*/
    
    
    cout << "# y = " << y << endl;
    cout << "# r [1/GeV]     Amplitude   \\partial_r   \\partial2 r" << endl;
    for (REAL r=N.MinR(); r<N.MaxR(); r*=1.1)
    {
        cout << r << " " << N.N(r, y) <<  " "
         << N.N(r,y,1) << " " << N.N(r,y,2) << endl;

    }
    

    return 0;
}
