/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2012
 */

#include <tools/tools.hpp>
#include <tools/config.hpp>
#include "amplitude.hpp"
#include "solver.hpp"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <csignal>

#include "ic.hpp"
#include "mv.hpp"
#include "ic_datafile.hpp"


const std::string version = "v. 0.9b  2012-xx-xx";

// We need global variables so that the signal handler works
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
    REAL minr = 1e-9;
    bool bfkl=false;
    
    if (argc>1)  if (string(argv[1])=="-help")
    {
        cout << "Usage: " << endl;
        cout << "-maxy y: set max rapidity" << endl;
        cout << "-minr minr: set smallest dipole size for the grid" << endl;
        cout << "-output file: save output to given file" << endl;
        cout << "-rc [CONSTANT,PARENT,BALITSKY,KW,MS]: set RC prescription" << endl;
        cout << "-ic [MV, FILE] params, MV params: qsqr anomalous_dim x0" << endl;
        cout <<"                        FILE params: filename x0" << endl;
        cout << "-alphas_scaling factor: scale \\lambdaQCD^2 by given factor" << endl;
        cout << "-ystep step: set rapidity step size" << endl;
        cout << "-bfkl: solve bfkl equation, no bk" << endl;
    }

    /*******************
     * Handle parameters
     ******************/
     
    N = new AmplitudeR();
    InitialCondition *ic=NULL;

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
            {
				cerr << "GBW is is not supported" << endl;
            }
            else if (string(argv[i+1])=="MV")
            {
				double qsqr, x0, gamma;
				qsqr = StrToReal(argv[i+2]);
				gamma = StrToReal(argv[i+3]);
				x0 = StrToReal(argv[i+4]);
				MV *tmpic = new MV();
				tmpic->SetQsqr(qsqr);
				tmpic->SetAnomalousDimension(gamma);
				tmpic->SetX0(x0);
				tmpic->SetLambdaQcd(0.241);
				N->SetLambdaQcd(0.241);
				N->SetInitialCondition(tmpic);  
				ic=tmpic;
			}
			else if (string(argv[i+1])=="FILE")
			{
				std::string fname = argv[i+2];
				IC_datafile *tmpic=new IC_datafile();
				tmpic->LoadFile(fname);
				tmpic->SetX0(StrToReal(argv[i+3]));
				N->SetInitialCondition(tmpic);
				ic = tmpic;
				
			}
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
        infostr << "# Scaling factor C^2 in alpha_s: " << N->GetAlphasScaling() << endl;
    
    infostr << "# Initial condition is " << N->GetInitialCondition()->GetString()
		<< endl;
     
    if (s.GetBfkl()) infostr <<"# Solving BFKL equation ";
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
    if (ic!=NULL)
		delete ic;
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
    out << "###" << N->GetInitialCondition()->X0() << endl;

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
