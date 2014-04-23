#include "gtest/gtest.h"
#include "CommandLineProcessor.h"

#include <string>
#include <vector>
#include <cstring>


void tokenizeArgumentString(const std::string& args, int& argc, char** argv[]);
std::vector<std::string> tokenize(const std::string& str, const char* delim);
char* copyToNewCharArray(const std::string& str);
char** copyToNewArrayOfCharArrays(const std::vector<std::string>& tokens);
void deleteArguments(int argc, char* argv[]);


TEST(CommandLineProcessor, Processing)
{
    int argc = 0;
    char** argv = NULL;

    // Convert the "command line" into argc and argv
    std::string argumentString
        ("./bamm -c divcontrol.txt --seed 1979 --overwrite 1");
    tokenizeArgumentString(argumentString, argc, &argv);

    // Simple test of the conversion
    EXPECT_EQ(argc, 7);

    CommandLineProcessor cmdLineProcessor(argc, argv);
    deleteArguments(argc, argv);

    // Test that the control file name was processed
    EXPECT_EQ("divcontrol.txt", cmdLineProcessor.controlFileName());

    // Test that parameters were processed
    const std::vector<UserParameter>& parameters =
        cmdLineProcessor.parameters();
    EXPECT_EQ("seed", parameters[0].first);
    EXPECT_EQ("1979", parameters[0].second);
    EXPECT_EQ("overwrite", parameters[1].first);
    EXPECT_EQ("1", parameters[1].second);
}


void tokenizeArgumentString(const std::string& args, int& argc, char** argv[])
{
    const std::vector<std::string>& tokens = tokenize(args, " \t");
    argc = tokens.size();
    *argv = copyToNewArrayOfCharArrays(tokens);
}


std::vector<std::string> tokenize(const std::string& str, const char* delim)
{
    std::vector<std::string> tokens;

    char* cStr = copyToNewCharArray(str);

    char* token = std::strtok(cStr, delim);
    while (token != NULL) {
        tokens.push_back(token);
        token = std::strtok(NULL, delim);
    }

    delete[] cStr;

    return tokens;
}


char* copyToNewCharArray(const std::string& str)
{
    char* cStr = new char[str.length() + 1];
    std::strcpy(cStr, str.c_str());
    return cStr;
}


char** copyToNewArrayOfCharArrays(const std::vector<std::string>& tokens)
{
    char** arrayOfCharArrays = new char*[tokens.size()];
    for (int i = 0; i < tokens.size(); i++) {
        arrayOfCharArrays[i] = copyToNewCharArray(tokens[i]);
    }
    return arrayOfCharArrays;
}


void deleteArguments(int argc, char* argv[])
{
    for (int i = 0; i < argc; i++) {
        delete[] argv[i];
    }

    delete[] argv;
}
