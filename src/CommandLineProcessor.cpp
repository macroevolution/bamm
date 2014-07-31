#include "CommandLineProcessor.h"
#include "Log.h"

#include <vector>
#include <string>


CommandLineProcessor::CommandLineProcessor(int argc, char* argv[])
{
    // Start at argv[1] because argv[0] is the program name
    for (int i = 1; i < argc; i += 2) {
        std::string argName(argv[i]);

        // Print help message
        if (argName == "-h" || argName == "--help") {
            exitWithMessage(usageText());
        }

        if (argName == "--version") {
            exitWithMessage(versionText());
        }

        // Every argument name must be followed by its value
        if (i + 1 == argc) {
            exitWithError(missingArgumentValueText());
        }

        std::string argValue(argv[i + 1]);

        // Control file argument is special, don't add it to parameter list
        if (argName == "-c" || argName == "--control") {
            _controlFileName = argValue;
        } else {
            // Every parameter name must start with two hyphens
            if (!startsWithTwoHyphens(argName)) {
                exitWithError(invalidArgumentNameText());
            }

            argName = argName.substr(2);    // Cut out the "--"
            _parameters.push_back(UserParameter(argName, argValue));
        }
    }

    // Print help message if control file was not specified
    if (_controlFileName == "") {
        exitWithMessage(usageText());
    }
}


std::string CommandLineProcessor::usageText() const
{
    return "Usage: bamm -c <control-file> "
        "[--<parameter-name> <parameter-value> ...]";
}


std::string CommandLineProcessor::versionText() const
{
    // "BAMM " is in a std::string so that its operator+() is used
    return std::string("BAMM ") + BAMM_VERSION +
        " (" + BAMM_VERSION_DATE + ")";
}


std::string CommandLineProcessor::missingArgumentValueText() const
{
    return "Missing value after parameter name";
}


bool CommandLineProcessor::startsWithTwoHyphens(const std::string& str) const
{
    return (str.length() > 2) && (str[0] == '-') && (str[1] == '-');
}


std::string CommandLineProcessor::invalidArgumentNameText() const
{
    return "Parameter name must start with \"--\"";
}


const std::string& CommandLineProcessor::controlFileName() const
{
    return _controlFileName;
}


const std::vector<UserParameter>& CommandLineProcessor::parameters() const
{
    return _parameters;
}
