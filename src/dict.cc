// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME srcdIdict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "interface/TrackTree.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *TrackPar_Dictionary();
   static void TrackPar_TClassManip(TClass*);
   static void *new_TrackPar(void *p = 0);
   static void *newArray_TrackPar(Long_t size, void *p);
   static void delete_TrackPar(void *p);
   static void deleteArray_TrackPar(void *p);
   static void destruct_TrackPar(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TrackPar*)
   {
      ::TrackPar *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::TrackPar));
      static ::ROOT::TGenericClassInfo 
         instance("TrackPar", "", 20,
                  typeid(::TrackPar), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &TrackPar_Dictionary, isa_proxy, 4,
                  sizeof(::TrackPar) );
      instance.SetNew(&new_TrackPar);
      instance.SetNewArray(&newArray_TrackPar);
      instance.SetDelete(&delete_TrackPar);
      instance.SetDeleteArray(&deleteArray_TrackPar);
      instance.SetDestructor(&destruct_TrackPar);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TrackPar*)
   {
      return GenerateInitInstanceLocal((::TrackPar*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TrackPar*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *TrackPar_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::TrackPar*)0x0)->GetClass();
      TrackPar_TClassManip(theClass);
   return theClass;
   }

   static void TrackPar_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_TrackPar(void *p) {
      return  p ? new(p) ::TrackPar : new ::TrackPar;
   }
   static void *newArray_TrackPar(Long_t nElements, void *p) {
      return p ? new(p) ::TrackPar[nElements] : new ::TrackPar[nElements];
   }
   // Wrapper around operator delete
   static void delete_TrackPar(void *p) {
      delete ((::TrackPar*)p);
   }
   static void deleteArray_TrackPar(void *p) {
      delete [] ((::TrackPar*)p);
   }
   static void destruct_TrackPar(void *p) {
      typedef ::TrackPar current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TrackPar

namespace ROOT {
   static TClass *vectorlETrackPargR_Dictionary();
   static void vectorlETrackPargR_TClassManip(TClass*);
   static void *new_vectorlETrackPargR(void *p = 0);
   static void *newArray_vectorlETrackPargR(Long_t size, void *p);
   static void delete_vectorlETrackPargR(void *p);
   static void deleteArray_vectorlETrackPargR(void *p);
   static void destruct_vectorlETrackPargR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TrackPar>*)
   {
      vector<TrackPar> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TrackPar>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TrackPar>", -2, "vector", 214,
                  typeid(vector<TrackPar>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETrackPargR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TrackPar>) );
      instance.SetNew(&new_vectorlETrackPargR);
      instance.SetNewArray(&newArray_vectorlETrackPargR);
      instance.SetDelete(&delete_vectorlETrackPargR);
      instance.SetDeleteArray(&deleteArray_vectorlETrackPargR);
      instance.SetDestructor(&destruct_vectorlETrackPargR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TrackPar> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TrackPar>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETrackPargR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TrackPar>*)0x0)->GetClass();
      vectorlETrackPargR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETrackPargR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETrackPargR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TrackPar> : new vector<TrackPar>;
   }
   static void *newArray_vectorlETrackPargR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TrackPar>[nElements] : new vector<TrackPar>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETrackPargR(void *p) {
      delete ((vector<TrackPar>*)p);
   }
   static void deleteArray_vectorlETrackPargR(void *p) {
      delete [] ((vector<TrackPar>*)p);
   }
   static void destruct_vectorlETrackPargR(void *p) {
      typedef vector<TrackPar> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TrackPar>

namespace {
  void TriggerDictionaryInitialization_dict_Impl() {
    static const char* headers[] = {
0    };
    static const char* includePaths[] = {
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.14.04-0d8dc/x86_64-slc6-gcc62-opt/include",
"/afs/cern.ch/work/m/meridian/MTD/H4Analysis/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class TrackPar;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#ifndef __TRACK_TREE__
#define __TRACK_TREE__

#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"

typedef unsigned long long int uint64;

using namespace std;

class TrackPar
{
 public:
  TrackPar() {};
  ~TrackPar() {};

  double x() { return  value[0]; }
  double y() { return  value[1]; }
  double alpha() { return  value[2]; }
  double beta() { return  value[3]; }

  double err_x() { return  sqrt(covariance[0]); }
  double err_y() { return  sqrt(covariance[2]); }
  double err_alpha() { return  sqrt(covariance[5]); }
  double err_beta() { return  sqrt(covariance[9]); }

  double corr_x_alpha()  { return covariance[3]/(err_x() * err_alpha()); }
  double corr_y_beta()  { return covariance[7]/(err_y() * err_beta()); }
  double corr_x_y() { return covariance[1]/(err_x() * err_y()); }
  double corr_alpha_beta() { return covariance[8]/(err_alpha() * err_beta()); }

  std::vector<double> value;
  std::vector<double> covariance;
};

class TrackTree
{
public:


  //---ctors---
  TrackTree(){};
  TrackTree(uint64* idx, TTree* tree=NULL);
  //---dtor---
  ~TrackTree(){};
  
  //---utils---
  void Init();
  void Clear() 
  {
    trackHits.clear();
    trackChi2.clear();
    fitStatus.clear();
    fitResult.clear();
  };

  void Fill() {tree_->Fill();};
  
  TTree*  tree_; 
  
  uint64* index;
  int n_tracks;
  std::vector<int> trackHits;
  std::vector<float> trackChi2;
  std::vector<int> fitStatus;
  std::vector<TrackPar> fitResult;
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TrackPar", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_dict() {
  TriggerDictionaryInitialization_dict_Impl();
}
