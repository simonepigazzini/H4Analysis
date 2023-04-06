#include "interface/PreProcessorBase.h"

void PreProcessorBase::Log(std::string message, LoggerLevel lv)
{
    //---Format message with preProcessor type and instance name
    std::map<LoggerLevel, std::string> colors = { {INFO, "\033[1;36m"}, 
                                                  {WARN, "\033[1;33m"},
                                                  {ERR,  "\033[1;31m"} };
    std::cout << colors[lv] << "["+preProcessorType_+"::"+instanceName_+"]: "
              << "\033[0m" << message << std::endl;

    return;
}
