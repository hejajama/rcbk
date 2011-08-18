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
#include <csignal>


const std::string version = "v. 0.1  2011-xx-xx";

// We need global variables so that the singla handler works
std::string output="output.dat";
std::stringstream infostr;
AmplitudeR* N;

void SaveData();
void SigIntHandler(int param);

int main(int argc, char* argv[])
{
    gsl_set_error_handler(&ErrHandler);
    std::signal(SIGINT, SigIntHandler);
    
    REAL maxy=1;
    REAL dy = 0.2;  // ystep
    RunningCoupling rc=CONSTANT;
    REAL alphas_scaling=1.0;
    InitialConditionR ic = GBW;
    REAL minr = 1e-9;
    bool bfkl=false;
    
    if (argc>1)  if (string(argv[1])=="-help")
    {
        cout << "Usage: " << endl;
        cout << "-maxy y: set max rapidity" << endl;
        cout << "-minr minr: set smallest dipole size for the grid" << endl;
        cout << "-output file: save output to given file" << endl;
        cout << "-rc [CONSTANT,PARENT,BALITSKY,KW,MS]: set RC prescription" << endl;
        cout << "-ic [GBW, MV, MV1, AN06]: set initial condition" << endl;
        cout << "-alphas_scaling factor: scale \\lambdaQCD^2 by given factor" << endl;
        cout << "-ystep step: set rapidity step size" << endl;
        cout << "-bfkl: solve bfkl equation, no bk" << endl;
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
            else if (string(argv[i+1])=="MS")
                rc = MS;
            else
            {
                cerr << "Unknown running coupling " << argv[i+1] << endl;
                return -1;
            }
        }
        else if (string(argv[i])=="-ic")
        {
            if (string(argv[i+1])=="GBW")
                ic = GBW;
            else if (string(argv[i+1])=="MV")
                ic = MV;
            else if (string(argv[i+1])=="MV1")
                ic = MV1;
            else if (string(argv[i+1])=="AN06")
                ic = AN06;
            else
            {
                cerr << "Unknown initial condition " << argv[i+1] << endl;
                return -1;
            }
        }
        else if (string(argv[i])=="-alphas_scaling")
            alphas_scaling = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-ystep")
            dy = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-minr")
            minr = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-bfkl")
            bfkl=true;
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }
    }
    
    N = new AmplitudeR();
    N->SetInitialCondition(ic);
    N->SetMinR(minr);
    N->Initialize();
    Solver s(N);
    s.SetBfkl(bfkl);
    s.SetRunningCoupling(rc);
    s.SetAlphasScaling(alphas_scaling);
    s.SetDeltaY(dy);

    infostr << "#";
    for (int i=0; i<argc; i++)
        infostr << argv[i] << " ";
    infostr << endl << "# BK equation solver " << version << endl;
    infostr << "# Running coupling: ";
    if (rc==CONSTANT) infostr << "constant, alphabar=" << ALPHABAR_s;
    if (rc==PARENT) infostr << "parent dipole";
    if (rc==BALITSKY) infostr << "Balitsky";
    if (rc==KW) infostr << "KW";
    if (rc==MS) infostr << "Motyka&Stasto, kin. constraing";
    infostr << endl;
    if (rc!=CONSTANT)
        infostr << "# Scaling factor in alpha_s: " << alphas_scaling << endl;
    
    infostr << "# Initial condition is ";
    if (ic==GBW)
        infostr << "GBW 1-exp(-r^2Q_s^2/4)";
    else if (ic == MV)
        infostr << "MV 1-exp(-(r^2 Q_s^2/4)^\\gamma log(1/r\\lambda_QCD + e) )";
    else if (ic == MV1)
        infostr << "MV 1-exp(-r^2 Q_s^2/4 log(1/r\\lambda_QCD + e) )";
    else if (ic==AN06)
        infostr << "AN06 1-exp(-(r^2Q_s^2)^\\gamma/4)";
    infostr << endl;
    infostr <<"# Initial saturation scale Q_s^2=" << N->InitialSaturationScaleSqr()
        << " GeV^2" << endl;
        
    if (bfkl) infostr <<"# Solving BFKL equation ";
    else infostr << "# Solving BK equation "; 
    infostr << "up to y=" << maxy << endl;
    infostr << "# r limits: " << N->MinR() << " - " << N->MaxR() << " points "
    << N->RPoints() << " multiplier " << N->RMultiplier() << endl;
    infostr << "# maxy " << maxy << " ystep " << dy << endl;
    cout << infostr.str() ;

    s.Solve(maxy);

    SaveData();

    cout << "# Done!" << endl;

    delete N;
    return 0;
}

void SaveData()
{
    cout << "Saving data in file " << output << endl;

    /*
     * Save data into a file
     * Syntax is described in file bk/README
     */

    std::ofstream out;
    out.open(output.c_str());
    out << infostr.str();

    out << "###" << std::scientific << std::setprecision(15) << N->MinR() << endl;
    out << "###" << std::scientific << std::setprecision(15) <<
        N->RMultiplier()  << endl;
    out << "###" << N->RPoints() << endl;

    for (int yind=0; yind<N->YPoints(); yind++)
    {
        out << "###" << std::scientific << std::setprecision(15)
            << N->YVal(yind) << endl;
        for (int rind=0; rind<N->RPoints(); rind++)
        {
            out << std::scientific << std::setprecision(15)
                << N->Ntable(yind, rind) << endl;
        }
    }
    out.close();
}

// User pressed ctrl+c or killed the program, save data
void SigIntHandler(int param)
{
    cerr << endl << "# Received SigInt signal, trying to save data..." << endl;
    infostr << "# Received SigInt signal, trying to save data..." << endl;

    SaveData();

    delete N;

    exit(1);
}
