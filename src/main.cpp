/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "tools.hpp"
#include "amplitude.hpp"
#include "solver.hpp"
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
    REAL maxy=1;
    std::string output="output.dat";
    RunningCoupling rc=CONSTANT;
    REAL alphas_scaling=1.0;
    
    if (argc>1)  if (string(argv[1])=="-help")
    {
        cout << "Usage: " << endl;
        cout << "-maxy y: set max rapidity" << endl;
        cout << "-output file: save output to given file" << endl;
        cout << "-rc [CONSTANT,PARENT,BALITSKY,KW]: set RC prescription" << endl;
        cout << "-alphas_scaling factor: scale \\lambdaQCD^2 by given factor" << endl;
    }

    /*******************
     * Handle parameters
     ******************/

    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-maxy")
            maxy = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-output")
            output = argv[i+1];
        else if (string(argv[i])=="-rc")
        {
            if (string(argv[i+1])=="CONSTANT")
                rc = CONSTANT;
            else if (string(argv[i+1])=="PARENT")
                rc = PARENT;
            else if (string(argv[i+1])=="BALITSKY")
                rc = BALITSKY;
            else if (string(argv[i+1])=="KW")
                rc = KW;
            else
            {
                cerr << "Unknown running coupling " << argv[i+1] << endl;
                return -1;
            }
        }
        else if (string(argv[i])=="-alphas_scaling")
            alphas_scaling = StrToReal(argv[i+1]);
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }
    }
    
    AmplitudeR N;
    Solver s(&N);
    s.SetRunningCoupling(rc);
    s.SetAlphasScaling(alphas_scaling);

    std::stringstream infostr;
    infostr << "#";
    for (int i=0; i<argc; i++)
        infostr << argv[i] << " ";
    infostr << endl << "# BK equation solver " << version << endl;
    infostr << "# Running coupling: ";
    if (rc==CONSTANT) infostr << "constant";
    if (rc==PARENT) infostr << "parent dipole";
    if (rc==BALITSKY) infostr << "Balitsky";
    if (rc==KW) infostr << "KW";
    infostr << endl;
    if (rc!=CONSTANT)
        infostr << "Scaling factor in alpha_s: " << alphas_scaling << endl;

    infostr << "# Solving BK equation up to y=" << maxy << endl;
    infostr << "# r limits: " << N.MinR() << " - " << N.MaxR() << " points "
    << N.RPoints() << " multiplier " << N.RMultiplier() << endl;
    cout << infostr.str() ;

    s.Solve(maxy);

    cout << "Saving data in file " << output << endl;

    /*
     * Save data into a file
     * Syntax is described in file bk/README
     */

    std::ofstream out;
    out.open(output.c_str());
    out << infostr.str();

    out << "###" << std::scientific << std::setprecision(15) << N.MinR() << endl;
    out << "###" << std::scientific << std::setprecision(15) <<
        N.RMultiplier()  << endl;
    out << "###" << N.RPoints() << endl;

    for (int yind=0; yind<N.YPoints(); yind++)
    {
        out << "###" << std::scientific << std::setprecision(15)
            << N.YVal(yind) << endl;
        for (int rind=0; rind<N.RPoints(); rind++)
        {
            out << std::scientific << std::setprecision(15)
                << N.Ntable(yind, rind) << endl;
        }
    }
    out.close();


    return 0;
}
