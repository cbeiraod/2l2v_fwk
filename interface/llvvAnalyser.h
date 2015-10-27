// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2015-10-27</date>
// <summary>Header file for the Analyser class</summary>
//
// <description>
//  Header file with the declarations of the llvvAnalyser class and related classes.
//  The base analyser from which other analyses can be derived from is declared within.
// </description>

#ifndef LLVV_ANALYSER_H
#define LLVV_ANALYSER_H

#include <map>
#include <vector>
#include <string>
#include <exception>
#include <iostream>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/FWLite/interface/Handle.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/llvvObjects.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/ValueWithSystematics.h"

class AnalyserException: public std::exception
{
public:
  AnalyserException(std::string mess): message(mess)
  {
  }

private:
  std::string message;

  virtual const char* what() const throw()
  {
    return message.c_str();
  }

protected:
};

class EventInfo
{
public:
  EventInfo();
  EventInfo(const EventInfo&) = delete; // Delete copy constructor
  EventInfo(EventInfo&&) = delete;      // Delete move constructor
  void Reset();
  inline void Lock()
  {
    isLocked = true;
    for(auto & kv : eventDoubles)
      kv.second.Lock();
    for(auto & kv : eventInts)
      kv.second.Lock();
    for(auto & kv : eventBools)
      kv.second.Lock();
  };
  #ifdef WITH_UNLOCK
  inline void Unlock()
  {
    isLocked = false;
    for(auto & kv : eventDoubles)
      kv.second.Unlock();
    for(auto & kv : eventInts)
      kv.second.Unlock();
    for(auto & kv : eventBools)
      kv.second.Unlock();
  };
  #endif

  ValueWithSystematics<double>& AddDouble(std::string name, double defaultVal);
  ValueWithSystematics<double>& GetDouble(std::string name);
  ValueWithSystematics<int>&    AddInt   (std::string name, int defaultVal);
  ValueWithSystematics<int>&    GetInt   (std::string name);
  ValueWithSystematics<bool>&   AddBool  (std::string name, bool defaultVal);
  ValueWithSystematics<bool>&   GetBool  (std::string name);
  
  void OutputEventListHeader(ofstream& file, const std::vector<std::string>& priority = std::vector<std::string>(0)) const;
  void OutputEventList(ofstream& file, const std::vector<std::string>& priority = std::vector<std::string>(0)) const;
  void SetSummaryTreeBranches(TTree* const tree);

private:
protected:
  bool isLocked;
  std::map<std::string,ValueWithSystematics<double>> eventDoubles;
  std::map<std::string,ValueWithSystematics<int>>    eventInts;
  std::map<std::string,ValueWithSystematics<bool>>   eventBools;
  
  template<class T>
  void OutputValueListHeader(ofstream& file, const ValueWithSystematics<T>& val, const std::string& name) const;

  template<class T>
  void OutputValueList(ofstream& file, const ValueWithSystematics<T>& val) const;

  template<class T>
  void AddBranch(TTree* const tree, ValueWithSystematics<T>& val, std::string name);

};

class Analyser
{
public:
  Analyser(std::string cfgFile);
  virtual ~Analyser();

  virtual void Setup();
  virtual void LoopOverEvents();

  inline void SetEventLimit(size_t val)      { limitEvents = val; };
  inline void SetDebugEvent(bool val)        { debugEvent = val; };
  inline void SetSkipEvents(int val)         { skipEvents = val; };
  inline void SetEventlistSelected(bool val) { eventlistSelected = val; };

  inline void RedirectCout() { saveRedirect = true; };
  inline void DiscardCout() { saveRedirect = false; };
  
  inline void KeepAllEvents() { keepAllEvents = true; };
  inline void KeepTriggeredEvents() { keepAllEvents = false; };
  
  inline void MergeBoostedTaus() { mergeBoostedTaus = true; };
  inline void DontMergeBoostedTaus() { mergeBoostedTaus = false; };

private:

protected:
  size_t limitEvents;
  bool debugEvent;
  int skipEvents;
  bool isSetup;
  std::ostream analyserCout;
  bool saveRedirect;
  bool keepAllEvents;
  bool mergeBoostedTaus;
  bool eventlistSelected;

  edm::ParameterSet cfgOptions;
  std::string outFile;
  std::string summaryOutFile;
  TFile* summaryOutTFile;
  TTree* summaryTree;
  std::string eventlistOutFile;
  ofstream eventListFile;
  SmartSelectionMonitor histMonitor;
  std::stringstream analyserBuffer;

  bool isMC;
  double crossSection;
  std::vector<std::string> fileList;
  std::string baseDir;
  std::string outDir;
  std::string jecDir;
  std::string pdfDir;
  bool runSystematics;
  bool saveSummaryTree;
  bool applyScaleFactors;
  bool debug;
  bool doDDBkg;
  bool outputEventList;
  bool isV0JetsMC;
  std::vector<double> pileupDistribution;
  
  llvvGenEvent genEv;
  fwlite::Handle<LHEEventProduct> LHEHandle;
  std::vector<bool> triggerBits;
  llvvGenParticleCollection gen;
  llvvLeptonCollection leptons;
  llvvTauCollection taus;
  llvvTauCollection boostedTaus;
  llvvJetCollection jets_;
  llvvJetExtCollection jets;
  llvvMet metVec;
  std::vector<int> triggerPrescales;
  
  ValueWithSystematics<double> xsecWeight;
  ValueWithSystematics<double> crossSection_;
  double PUNorm[3];
  edm::LumiReWeighting* LumiWeights;
  utils::cmssw::PuShifter_t PuShifters;
  MuScleFitCorrector* muCor;
  FactorizedJetCorrector* jesCor;
  JetCorrectionUncertainty* totalJESUnc;
  
  std::vector<TString> chTags;
  
  EventInfo eventContent;

  virtual void LoadCfgOptions();
  virtual void InitHistograms();
  virtual void EventContentSetup();
  virtual void ProcessEvent(size_t iev);
  virtual void FillHistograms();
  
  virtual void UserLoadCfgOptions() = 0;
  virtual void UserSetup() = 0;
  virtual void UserProcessEvent(size_t iev) = 0;
  virtual void UserInitHistograms() = 0;
  virtual void UserEventContentSetup() = 0;
  virtual void UserFillHistograms() = 0;

  template<class T>
  void loadSystematics(std::vector<std::string>& list, ValueWithSystematics<T> variable);

  ValueWithSystematics<TLorentzVector> getMETvariations();
  ValueWithSystematics<TLorentzVector> MET;
  
  std::vector<PDFInfo*> pdfVariations;

};

#include "UserCode/llvv_fwk/interface/llvvAnalyser.hpp"

#endif

