#include "interface/PluginBase.h"

//**********Getters***********************************************************************
vector<SharedData> PluginBase::GetSharedData(string tag, string type, bool permanent)
{
    vector<SharedData> outData;
    for(auto& data : data_)
    {
        if(data.obj && data.permanent==permanent && (tag=="" || data.tag==tag) &&
           (type=="" || type==data.obj->IsA()->GetName()))
            outData.push_back(data);
    }

    return outData;
}

//**********Utils (private)***************************************************************
void PluginBase::RegisterSharedData(TObject* obj, string tag, bool isPermanent)
{
    if(obj)
    {
        std::cout << "Registering Shared Data " << this->instanceName_+"_"+tag << std::endl;
        data_.push_back(SharedData{isPermanent, this->instanceName_+"_"+tag, obj});
    }
    return;
}

void PluginBase::Log(std::string message, LoggerLevel lv)
{
    //---Format message with plugin type and instance name
    std::map<LoggerLevel, std::string> colors = { {INFO, "\033[1;36m"}, 
                                                  {WARN, "\033[1;33m"},
                                                  {ERR,  "\033[1;31m"} };
    std::cout << colors[lv] << "["+pluginType_+"::"+instanceName_+"]: "
              << "\033[0m" << message << std::endl;

    return;
}
