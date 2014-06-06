#ifndef TOOLS_H
#define TOOLS_H


#include "Log.h"

#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>


std::vector<std::string> split_string(const std::string& str, char delim);

template<typename T> T convert_string(const std::string& str)
{
    std::istringstream iss(str);

    T value;
    iss >> value;

    if (iss.fail()) {
        log(Error) << "Could not convert <<" << str
            << ">> to the proper type.\n";
        std::exit(1);
    }

    return value;
}


template<typename T> T sqr(const T& x)
{
    return x * x;
}


#endif
