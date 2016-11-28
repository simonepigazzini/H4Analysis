CXX = g++
CXXFLAGS = -std=c++1y -fPIC
SOFLAGS = -shared -O3
INCLUDE = -I"./" 
LIB = -L"./lib/" -L"./DynamicTTree/lib" -L"./CfgManager/lib" -Wl,-rpath=lib/:DynamicTTree/lib/:CfgManager/lib/ -lDTT 

ROOT_LIB := `root-config --libs --glibs`
ROOT_FLAGS := `root-config --cflags --ldflags` -lMathCore -lMathMore

DEPS = Makefile DynamicTTree/interface/DynamicTTreeBase.h DynamicTTree/interface/DynamicTTreeInterface.h \
	CfgManager/interface/CfgManager.h CfgManager/interface/CfgManagerT.h \
	interface/WFClass.h
DEPS_OBJS = lib/utils.o lib/WFClass.o lib/WFClassNINO.o lib/FFTClass.o lib/WFViewer.o lib/MCPAnalyzer.o \
	lib/H4Tree.o lib/RecoTree.o lib/DigiTree.o lib/WFTree.o lib/PositionTree.o lib/PluginBase.o lib/H4Dict.so
PLUG_DEPS = lib/PluginBase.o
PLUG_OBJS = lib/libDigitizerReco.so lib/libMakeCovarianceMatrix.so \
	lib/libHodoReco.so lib/libHodoBTFReco.so lib/libWireChamberReco.so \
	lib/libInfoTreeMaker.so lib/libADCReco.so lib/libWFAnalyzer.so lib/libFFTAnalyzer.so
DICT_OBJS = lib/WFViewer.o lib/MCPAnalyzer.o

MAIN = bin/H4Reco 

all: dynamicTree cfgManager $(DEPS_OBJS) $(PLUG_OBJS) $(MAIN) lib/LinkDef.cxx 

lib/%.o: src/%.cc interface/%.h $(DEPS)
	@echo " CXX $<"
	@$ $(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDE) $(ROOT_LIB) $(ROOT_FLAGS)

lib/lib%.so: plugins/%.cc plugins/%.h $(PLUG_DEPS) $(DEPS_OBJS)
	@echo " CXX $<"
	@$ $(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $< $(INCLUDE) $(ROOT_LIB) $(ROOT_FLAGS) $(LIB)

lib/LinkDef.cxx: interface/WFViewer.h interface/MCPAnalyzer.h interface/LinkDef.h 
	@$ rootcling -f $@ -c $^
	@echo " rootcling $<"

lib/H4Dict.so: lib/LinkDef.cxx $(DICT_OBJS)
	@echo " CXX $<"
	@$ $(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $^ $(INCLUDE) $(ROOT_LIB) $(ROOT_FLAGS) $(LIB)

bin/%: main/%.cpp $(DEPS_OBJS) interface/PluginLoader.h CfgManager/lib/libCFGMan.so
	@echo " CXX $<"
	@$ $(CXX) $(CXXFLAGS) -ldl -o $@ $^ $(INCLUDE) $(ROOT_LIB) $(ROOT_FLAGS) $(LIB)

CfgManager/lib/libCFGMan.so: cfgManager

dynamicTree:
	cd DynamicTTree && $(MAKE)

cfgManager:
	cd CfgManager && $(MAKE)

clean:
	rm -fr tmp/*
	rm -fr lib/*
	rm -fr bin/*
	cd DynamicTTree && $(MAKE) clean
	cd CfgManager && $(MAKE) clean
