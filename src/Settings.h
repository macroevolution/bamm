#ifndef SETTINGS_H
#define SETTINGS_H

#include "SettingsParameter.h"
#include "Log.h"

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <cstdlib>


typedef std::pair<std::string, SettingsParameter> Parameter;
typedef std::map<std::string, SettingsParameter> ParameterMap;

typedef std::pair<std::string, std::string> UserParameter;


class Settings
{
public:

    Settings(const std::string& controlFilename,
        const std::vector<UserParameter>& commandLineParameters);

    std::string get(const std::string& name) const;
    template<typename T> T get(const std::string& name) const;

    void set(const std::string& name, const std::string& value);
  
    void printCurrentSettings(std::ostream& out = std::cout) const;

private:

    void readControlFile(const std::string& controlFilename);

    void initializeGlobalSettings();
    void initializeSpeciationExtinctionSettings();
    void initializeTraitSettings();
    void initializeSettingsWithUserValues();

    void checkAllSettingsAreUserDefined() const;
    void checkAllOutputFilesAreWriteable() const;

    void assertNotUserDefined(const SettingsParameter& parameter) const;
    void addParameter(const std::string& name, const std::string& value,
        UserDefinedStatus userDefined = Required,
        DeprecationStatus deprecated = NotDeprecated);

    void attachPrefixToOutputFiles();
    std::string attachPrefix
        (const std::string& prefix, const std::string& str) const;
    std::string extractDir(const std::string& path) const;
    std::string extractFileName(const std::string& path) const;

    bool anyOutputFileExists() const;
    bool fileExists(const std::string& filename) const;

    void exitWithErrorNoControlFile() const;
    void exitWithErrorInvalidLine(const std::string& line) const;
    void exitWithErrorUndefinedParameter(const std::string& name) const;
    void exitWithErrorInvalidModelType() const;
    void exitWithErrorParametersNotFound
        (const std::vector<std::string>& paramsNotFound) const;
    void exitWithErrorParameterIsDeprecated(const std::string& param) const;
    void exitWithErrorDuplicateParameter(const std::string& param) const;
    void exitWithErrorOutputFileExists() const;

    static const size_t NumberOfParamsToPrefix = 10;
 
    // Parameters that settings knows about
    ParameterMap _parameters;
    
    // Parameters read from the control file
    std::vector<UserParameter> _userParameters;

    // Parameters read from the command line
    std::vector<UserParameter> _commandLineParameters;
    
    // function to handle the validation of settings for
    //   expanded oct 2015 options
    void validateSettings(void);
    
    
};


inline std::string Settings::get(const std::string& name) const
{
    return get<std::string>(name);
}


template<typename T>
inline T Settings::get(const std::string& name) const
{
    ParameterMap::const_iterator it = _parameters.find(name);
    if (it != _parameters.end()) {
        return (it->second).value<T>();
    } else {
        log(Error) << "Parameter <<" << name << ">> does not exist.\n";
        std::exit(1);
    }
}





#endif
