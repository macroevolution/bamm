#ifndef COMMAND_LINE_PROCESSOR_H
#define COMMAND_LINE_PROCESSOR_H


#include <vector>
#include <string>

typedef std::pair<std::string, std::string> UserParameter;


class CommandLineProcessor
{
public:

    CommandLineProcessor(int argc, char* argv[]);

    const std::string& controlFileName() const;
    const std::vector<UserParameter>& parameters() const;

private:

    std::string usageText() const;
    std::string versionText() const;
    std::string missingArgumentValueText() const;
    bool startsWithTwoHyphens(const std::string& str) const;
    std::string invalidArgumentNameText() const;

    std::string _controlFileName;
    std::vector<UserParameter> _parameters;
};


#endif
