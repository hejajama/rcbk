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
using std::string;
const string version = "v. 0.1  2011-xx-xx";

enum Mode
{
    X,
    K,
    SATSCALE,
    GD
};

int main(int argc, char* argv[])
{
    std::stringstream infostr;
    infostr << "#";
    for (int i=0; i<argc; i++)
        infostr << argv[i] << " ";
    cout << infostr.str() << endl;
    
    gsl_set_error_handler(&ErrHandler);

    Mode mode=X;
    REAL Ns=0.5;
    REAL y=0;
    string datafile="output.dat";

    if (string(argv[1])=="-help")
    {
        cout << "-y y: set rapidity" << endl;
        cout << "-data datafile" << endl;
        cout << "-x: print amplitude in x space" << endl;
        cout << "-k: print amplitude in k space" << endl;
        cout << "-satscale Ns, print satscale r_s defined as N(r_s)=Ns" << endl;
        return 0;
    }
    
    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-y")
            y = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-data")
            datafile = argv[i+1];
        else if (string(argv[i])=="-x")
            mode=X;
        else if (string(argv[i])=="-k")
            mode=K;
        else if (string(argv[i])=="-satscale")
        {
            mode=SATSCALE;
            Ns = StrToReal(argv[i+1]);
        }
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }
    }

    cout << "# Reading data from file " << datafile << endl;
    AmplitudeLib N(datafile);
    N.InitializeInterpoaltion(y);
    cout << "# y = " << y << endl;

    if (mode==K)
    {
        cout << "# k [GeV]     Amplitude  " << endl;
        for (REAL k=1.0/N.MaxR(); k<1.0/N.MinR(); k*=1.1)
        {
            cout << k << " " << N.N_k(k, y) << endl;

        }
    } else if (mode==X)
    {
        cout << "# r [1/GeV]     Amplitude   \\partial_r   \\partial2"
         << " r d ln N / d ln r^2" << endl;
        for (REAL r=N.MinR()*1.01; r<N.MaxR(); r*=1.1)
        {
            cout << r << " " << N.N(r, y) <<  " "
             << N.N(r,y,1) << " " << N.N(r,y,2) <<
             " " << N.LogLogDerivative(r,y) << endl;
        }
    }
    else if (mode==SATSCALE)
    {
        cout <<"# Saturation scale N(r_s) = " << Ns << endl;
        cout <<"# y    r_s [1/GeV]" << endl;
        for (REAL y=0; y < N.MaxY(); y+=0.1)
        {
            cout << y << " " << N.SaturationScale(y, Ns) << endl;
        }
    }
    else if (mode==GD)
    {
        UGD ugd(&N);

        for (REAL k=0.1; k<10; k*=1.1)
        {
            cout << k << " " << ugd.Evaluate(k, y) << endl;
        }
    }
    else
    {
        cerr << "Unkown mode " << argv[3] << endl;
        return -1;
    }
    return 0;
    

}
