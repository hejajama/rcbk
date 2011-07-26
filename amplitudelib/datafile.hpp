/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _DATAFILE_H
#define _DATAFILE_H

#include "../src/tools.hpp"
#include "../src/config.hpp"
#include <sstream>
#include <vector>
/* Read data from datafiles, file format is defined in file README
 */

class DataFile
{
    public:
        DataFile(string fname);
        REAL MinR();
        REAL RMultiplier();
        int RPoints();
		REAL MaxY();

        void GetData(std::vector< std::vector<REAL> > &n,
            std::vector<REAL> &rapidities);

    private:
        string filename;
        std::vector<std::vector <REAL> > data;
        std::vector<REAL> yvals;
        REAL minr;
        REAL r_multiplier;
        int rpoints;

};

#endif
