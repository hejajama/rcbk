#include "config.hpp"
#include <iostream>
#include <sstream>


double StrToReal(std::string str)
{
    std::stringstream buff(str);
    double tmp;
    buff >> tmp;
    return tmp;
}

int StrToInt(std::string str)
{
    std::stringstream buff(str);
    int tmp;
    buff >> tmp;
    return tmp;
}

