#ifndef SETTINGS_PARAMETER_H
#define SETTINGS_PARAMETER_H


#include <string>
#include <sstream>


class SettingsParameter
{

public:

    SettingsParameter(const std::string& name, const std::string& defaultValue,
        bool mustBeUserDefined = true);

    SettingsParameter(const SettingsParameter& other);
    SettingsParameter& operator=(const SettingsParameter& other);

    const std::string& name() const;
    // No setter: the parameter name should never be changed.

    template <typename T> T value() const;
    template <typename T> void setValue(const T& value);

    void setStringValue(const std::string& value);

    // By default, parameters must be user-defined, but some parameters
    // are optional and should instead use their default value
    bool mustBeUserDefined() const;

    // The parameter becomes user-defined when the user sets its value
    // (even if it is the same as the default value).
    bool isUserDefined() const;

private:

    void exitWithErrorWrongType() const;

    std::string _name;
    std::string _value;

    bool _mustBeUserDefined;
    bool _isUserDefined;

};


template <typename T>
T SettingsParameter::value() const
{
    // Convert string value to the specified type
    std::istringstream valueStream(_value);
    T typedValue;
    valueStream >> typedValue;

    if (valueStream.fail())
        exitWithErrorWrongType();

    return typedValue;
}


template <> inline    // Must be inline for linker to work
std::string SettingsParameter::value<>() const
{
    return _value;
}


template <typename T>
void SettingsParameter::setValue(const T& value)
{
    // Convert value to a string
    std::ostringstream valueStream;
    valueStream << value;

    if (valueStream.fail())
        exitWithErrorWrongType();

    _value = valueStream.str();
    _isUserDefined = true;
}


template <> inline    // Must be inline for linker to work
void SettingsParameter::setValue<>(const std::string& value)
{
    _value = value;
    _isUserDefined = true;
}


#endif
