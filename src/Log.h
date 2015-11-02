#ifndef LOG_H
#define LOG_H

#include <iostream>
#include <string>


#define TEXT_COLOR_DEFAULT "\x1b[39m"    // Default color (gray)
#define TEXT_COLOR_MESSAGE "\x1b[39m"    // Default color (gray)
#define TEXT_COLOR_WARNING "\x1b[33m"    // Yellow
#define TEXT_COLOR_ERROR   "\x1b[31m"    // Red


enum LogType
{
    Message,
    Warning,
    Error
};


class Log
{

public:

    static Log& instance();

    std::ostream& outputStream(LogType logType, std::ostream& out);

private:

    Log() {};

    void startMessage(std::ostream& out);
    void startWarning(std::ostream& out);
    void startError(std::ostream& out);

    static Log _logger;

};


std::ostream& log(LogType logType = Message);
std::ostream& log(std::ostream& out);
std::ostream& log(LogType logType, std::ostream& out);

void exitWithMessage(const std::string& message);
void exitWithError(const std::string& message);


#endif
