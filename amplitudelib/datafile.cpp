/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "../src/tools.hpp"
#include "../src/config.hpp"
#include "datafile.hpp"
#include <fstream>
#include <sstream>
using std::ifstream;
using std::getline;
using std::stringstream;

DataFile::DataFile(string fname)
{
    filename=fname;
    ifstream file(fname.c_str());
    if (!file.is_open())
    {
        cerr << "ERROR! Coudn't read file " << fname << endl;
        return;
    }
    int confid=0;
    while(!file.eof() and confid < 3)
    {
        string line;
        getline(file, line);
        if (line.substr(0, 3)=="###")
        {                    
            switch (confid)
            {
                case 0:
                    minr = StrToReal(line.substr(3,line.length()-3));
                    break;
                case 1:
                    r_multiplier = StrToReal(line.substr(3,line.length()-3));
                    break;
                case 2:
                    rpoints = StrToInt(line.substr(3,line.length()-3));
                    break;
                default:
                    cerr << "File " << fname << " is formatted incorrectly!" << endl;
                    break;
            }
            confid++; 
        }
    }

    //TODO: It's impossible that this condition doesn't hold
    if (confid < 3)
        cerr << "File " << fname << " doesn't have enough metadata!" << endl;

    // Ok, configurations are read, then read all yvals
    REAL y=-1;
    std::vector<REAL> tmpvec;
    while (!file.eof())
    {
        string line;
        getline(file, line);

        // New rapidity?
        if (line.substr(0,3)=="###")
        {
            if (tmpvec.size()>0)
                data.push_back(tmpvec);

            if (tmpvec.size()>0 and tmpvec.size() != rpoints)
            {
                cerr << "File " << fname << ": read " << tmpvec.size() << " points, but "
                << "there should have been " << rpoints << " points, y=" 
                << y << ". " << LINEINFO << endl;
            }

            y = StrToReal(line.substr(3,line.length()-3));
            yvals.push_back(y);
            tmpvec.clear();
            continue;   // Next line is probably new amplitude value
        }
        else if (line.substr(0,1)=="#")
            continue;   // Comment

        // Ok, so this a new amplitude value
        tmpvec.push_back(StrToReal(line));
    }

    // Add last entry
    ///TODO: tmpvec[rpoints+1] \approx 0 ???
    data.push_back(tmpvec);

    if (data[0].size() != rpoints)
    {
        cerr << "File " << fname << ": read " << data.size() << " rpoints, but "
        << "there should have been " << rpoints << " points???" << LINEINFO << endl;
    }
}

void DataFile::GetData(std::vector< std::vector<REAL> > &n,
                        std::vector<REAL> &rapidities)
{
	n.clear();
    rapidities.clear();
    // Return vector where indexes are vec[y][r] containing amplitude

    for (uint yind=0; yind < data.size(); yind++)
    {
        std::vector<REAL> tmpvec;
        for (uint rind=0; rind<rpoints; rind++)
        {
            tmpvec.push_back(data[yind][rind]);
        }
        n.push_back(tmpvec);

        rapidities.push_back(yvals[yind]);
        
    }
    
}

REAL DataFile::MinR()
{
    return minr;
}

REAL DataFile::RMultiplier()
{
    return r_multiplier;
}

int DataFile::RPoints()
{
        return rpoints;
}

REAL DataFile::MaxY()
{
	return yvals[yvals.size()-1];
}

