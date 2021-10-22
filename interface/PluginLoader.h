#ifndef __PLUGIN_LOADER__
#define __PLUGIN_LOADER__

#include <iostream>
#include <sstream>
#include <dlfcn.h>
#include <string>

using namespace std;

template<class P> class PluginLoader
{
public:
    //---ctors---
    PluginLoader(string plugin_name);

    //---dtor
    ~PluginLoader();

    //---utils---
    void Create() {
        instance = pluginCreator(); 
        instance->SetPluginType(pluginName_); 
    };
    P*   CreateInstance(string name) {
        instance->SetInstanceName(name); 
        return instance;
    };
    void Destroy() {
        pluginDestroyer(instance);
    };
    
private:
    string pluginName_;
    void* pluginHandle;
    P*    (*pluginCreator)(void);
    void  (*pluginDestroyer)(P*);
    P*    instance;
};

//----------Constructor-------------------------------------------------------------------
template<class P> inline PluginLoader<P>::PluginLoader(string plugin_name)
{
    pluginName_ = plugin_name;
    char* error;
    pluginHandle = dlopen(("lib"+plugin_name+".so").c_str(), RTLD_LAZY);
    if(!pluginHandle)
    {
        cout << ">>> PluginLoader: " << dlerror() << endl;
        exit(1);
    }
    pluginCreator = (P* (*)(void))dlsym(pluginHandle, "create");
    if((error = dlerror()) != NULL)
    {
        cout << ">>> PluginLoader: " << error << endl;
        exit(1);
    }
    pluginDestroyer = (void (*)(P*))dlsym(pluginHandle, "destroy");
    if((error = dlerror()) != NULL)
    {
        cout << ">>> PluginLoader: " << error << endl;
        exit(1);
    }
}    

//----------Destructor--------------------------------------------------------------------
template<class P> inline PluginLoader<P>::~PluginLoader()
{
    dlclose(pluginHandle);
}

#endif
