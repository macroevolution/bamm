#include "Tools.h"

#include <vector>
#include <string>
#include <sstream>


std::vector<std::string> split_string(const std::string& str, char delim)
{
    std::vector<std::string> elems;
    std::stringstream ss(str);
    std::string item;

    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }

    return elems;
}
