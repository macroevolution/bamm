#include "Tools.h"


void ConfigFileReader::read(std::istream& in) const
{
    while (in) {
        std::string line = readLine(in);

        line = stripWhiteSpace(line);
        line = stripComments(line);

        if (line.size() > 0) {
            _settings.push_back(parseConfigLine(line));
        }
    }
}


std::string ConfigFileReader::readLine(std::istream& in) const
{
    std::string line;
    std::getline(in, line);
    return line;
}


std::string ConfigFileReader::stripWhiteSpace(std::string line) const
{
    line.erase(std::remove_if(line.begin(), line.end(),
        (int(*)(int))isspace), line.end());
    return line;
}


std::string ConfigFileReader::stripComments(std::string line) const
{
    std::istringstream iss(line);
    getline(iss, line, '#');
    return line;
}


RawSetting ConfigFileReader::parseConfigLine(const std::string& line) const
{
    const std::vector<std::string>& tokens = split_string(line, '=');

    if (tokens.size() != 2) {
        throw InvalidLineError();
    }

    return RawSetting(tokens[0], tokens[1]);
}
