#include "Log.h"
#include <iostream>


Log Log::_logger;


Log& Log::instance()
{
    return _logger;
}


std::ostream& Log::outputStream(LogType logType, std::ostream& out)
{
    if (logType == Message) {
        // Do nothing
    } else if (logType == Warning) {
        startWarning(out);
    } else if (logType == Error) {
        startError(out);
    }

    return out;
}


void Log::startWarning(std::ostream& out)
{
    out << TEXT_COLOR_WARNING << "\nWARNING: " << TEXT_COLOR_DEFAULT;
}


void Log::startError(std::ostream& out)
{
    out << TEXT_COLOR_ERROR << "\nERROR: " << TEXT_COLOR_DEFAULT;
}


std::ostream& log(LogType logType)
{
    std::ostream* out;

    if (logType == Message) {
        out = &std::cout;
    } else if (logType == Warning) {
        out = &std::cerr;
    } else if (logType == Error) {
        out = &std::cerr;
    } else {    // Should never happen
        out = &std::cout;
    }

    return log(logType, *out);
}


std::ostream& log(std::ostream& out)
{
    return log(Message, out);
}


std::ostream& log(LogType logType, std::ostream& out)
{
    return Log::instance().outputStream(logType, out);
}
