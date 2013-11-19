#include "SettingsParameter.h"
#include <iostream>
#include <string>
#include <cstdlib>


SettingsParameter::SettingsParameter(const std::string& name,
    const std::string& defaultValue, bool mustBeUserDefined) :
        _name(name), _value(defaultValue),
        _mustBeUserDefined(mustBeUserDefined), _isUserDefined(false)
{
}


SettingsParameter::SettingsParameter(const SettingsParameter& other) :
    _name(other._name), _value(other._value),
    _mustBeUserDefined(other._mustBeUserDefined),
    _isUserDefined(other._isUserDefined)
{
}


SettingsParameter& SettingsParameter::operator=(const SettingsParameter &other)
{
    _name = other._name;
    _value = other._value;
    _mustBeUserDefined = other._mustBeUserDefined;
    _isUserDefined = other._isUserDefined;
    return *this;
}


const std::string& SettingsParameter::name() const
{
    return _name;
}


void SettingsParameter::setStringValue(const std::string& value)
{
    _value = value;
    _isUserDefined = true;
}


bool SettingsParameter::mustBeUserDefined() const
{
    return _mustBeUserDefined;
}


bool SettingsParameter::isUserDefined() const
{
    return _isUserDefined;
}


void SettingsParameter::exitWithErrorWrongType() const
{
    std::cout << "ERROR: Parameter " << _name << " has the wrong type.\n";
    std::cout << "Fix by assigning the parameter a value of the right type.\n";
    std::exit(1);
}
