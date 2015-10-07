#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for SVfit


#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
//#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"
#include "UserCode/llvv_fwk/interface/llvvObjects.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"


#include "TROOT.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TNtuple.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRotation.h"

#include <iostream>
#include <ios>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <boost/shared_ptr.hpp>
#include <Math/VectorUtil.h>
#include <bitset>
#include <cctype>
#include <cmath>
#include <exception>
#include <algorithm>

// Include MT2 library:
// http://particle.physics.ucdavis.edu/hefti/projects/doku.php?id=wimpmass    ** Code from here
// http://www.hep.phy.cam.ac.uk/~lester/mt2/    ** Other libraries
#include "UserCode/llvv_fwk/interface/mt2_bisect.h"


//struct MatchPathSeparator
//{
//  bool operator()( char ch ) const
//  {
//    return ch == '/';
//  }
//};

//std::string basename( std::string const& pathname )
//{
//  return std::string( 
//        std::find_if( pathname.rbegin(), pathname.rend(),
//                      MatchPathSeparator() ).base(),
//        pathname.end() );
//}


#ifndef DEBUG_EVENT
#define DEBUG_EVENT true
#endif
//#define WITH_UNLOCK

#define NAN_WARN(X) if(std::isnan(X)) std::cout << "  Warning: " << #X << " is nan" << std::endl;
#define EVENTLISTWIDTH 15

class AnalyserException: public exception
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

template<class T>
class ValueWithSystematicsInternal
{
public:
  ValueWithSystematicsInternal(T val = T(0));
  ValueWithSystematicsInternal(const ValueWithSystematicsInternal<T>& val); // Copy constructor
  virtual void Reset();
  inline void Lock() { isLocked = true; };
  #ifdef WITH_UNLOCK
  inline void Unlock() { isLocked = false; };
  #endif
  
  bool AddMetadata(const std::string& key, const std::string& value);
  std::string GetMetadata(std::string& key);
  std::string GetMetadata(const std::string& key) const;

  inline T& Value() { return value; };
  inline const T& Value() const { return value; };
  inline T& DefaultValue() { return defaultValue; };
  T& Systematic(const std::string& name);
  const T& Systematic(const std::string& name) const;
  T& GetSystematicOrValue(const std::string& name);
  const T& GetSystematicOrValue(const std::string& name) const;
  
  // These next operators are where the magic happens
  // ---------  Function  operator  ---------
  T& operator()(const std::string& name);
  // ---------  Casting   operator  ---------
  explicit operator T () const; // Can only be explicitly called, ie. it disables implicit calling
  // --------- Assignment operators ---------
  ValueWithSystematicsInternal<T>& operator= (const T& val);
  ValueWithSystematicsInternal<T>& operator= (const ValueWithSystematicsInternal<T>& val);
  //       Compound Assignment Operators
  ValueWithSystematicsInternal<T>& operator+=(const T& val);
  ValueWithSystematicsInternal<T>& operator+=(const ValueWithSystematicsInternal<T>& val);
  ValueWithSystematicsInternal<T>& operator-=(const T& val);
  ValueWithSystematicsInternal<T>& operator-=(const ValueWithSystematicsInternal<T>& val);
  ValueWithSystematicsInternal<T>& operator*=(const T& val);
  ValueWithSystematicsInternal<T>& operator*=(const ValueWithSystematicsInternal<T>& val);
  ValueWithSystematicsInternal<T>& operator/=(const T& val);
  ValueWithSystematicsInternal<T>& operator/=(const ValueWithSystematicsInternal<T>& val);
  // --------- Arithmetic operators ---------
  const ValueWithSystematicsInternal<T> operator+(const T& val) const;
  const ValueWithSystematicsInternal<T> operator+(const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematicsInternal<T> operator-(const T& val) const;
  const ValueWithSystematicsInternal<T> operator-(const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematicsInternal<T> operator*(const T& val) const;
  const ValueWithSystematicsInternal<T> operator*(const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematicsInternal<T> operator/(const T& val) const;
  const ValueWithSystematicsInternal<T> operator/(const ValueWithSystematicsInternal<T>& val) const;
  // --------- Comparison operators ---------
  const ValueWithSystematicsInternal<bool> operator==(const T& val) const;
  const ValueWithSystematicsInternal<bool> operator==(const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematicsInternal<bool> operator!=(const T& val) const;
  const ValueWithSystematicsInternal<bool> operator!=(const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematicsInternal<bool> operator> (const T& val) const;
  const ValueWithSystematicsInternal<bool> operator> (const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematicsInternal<bool> operator< (const T& val) const;
  const ValueWithSystematicsInternal<bool> operator< (const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematicsInternal<bool> operator>=(const T& val) const;
  const ValueWithSystematicsInternal<bool> operator>=(const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematicsInternal<bool> operator<=(const T& val) const;
  const ValueWithSystematicsInternal<bool> operator<=(const ValueWithSystematicsInternal<T>& val) const;
  // ---------  Logical  operators  ---------
  const ValueWithSystematicsInternal<T> operator! () const;
  const ValueWithSystematicsInternal<T> operator&&(const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematicsInternal<T> operator||(const ValueWithSystematicsInternal<T>& val) const;
  // ---------   Unary  operators   ---------
  const ValueWithSystematicsInternal<T> operator-() const;
  ValueWithSystematicsInternal<T>& operator++();
  ValueWithSystematicsInternal<T> operator++(int);
  ValueWithSystematicsInternal<T>& operator--();
  ValueWithSystematicsInternal<T> operator--(int);
  
  inline typename std::map<std::string, T>::iterator& begin() { return systematics.begin(); };
  inline typename std::map<std::string, T>::iterator& end() { return systematics.end(); };
  std::vector<std::string> Systematics() const;

private:
protected:
  bool isLocked;
  T defaultValue;
  T value;
  mutable std::map<std::string, T> systematics;
  std::map<std::string,std::string> metadata;

};

template<class T>
ValueWithSystematicsInternal<T>::ValueWithSystematicsInternal(T val): isLocked(false), defaultValue(val), value(val)
{
}

template<class T>
ValueWithSystematicsInternal<T>::ValueWithSystematicsInternal(const ValueWithSystematicsInternal<T>& val): isLocked(false), defaultValue(val.defaultValue), value(val.value), systematics(val.systematics), metadata(val.metadata)
{
}

template<class T>
void ValueWithSystematicsInternal<T>::Reset()
{
  value = defaultValue;
  if(isLocked)
  {
    for(auto& kv : systematics)
      kv.second = defaultValue;
  }
  else
  {
    systematics.clear();
//    metadata.clear();
  }
}

template<class T>
bool ValueWithSystematicsInternal<T>::AddMetadata(const std::string& key, const std::string& value)
{
  if(metadata.count(key) != 0)
    std::cout << "Metadata already exists with that key, it will be overwritten. Old value: \"" << metadata[key] << "\"" << std::endl;
  else
    if(isLocked)
      throw AnalyserException("Unable to add metadata \""+key+":"+value+"\" after locking.");

  metadata[key] = value;
  return true;
}

template<class T>
std::string ValueWithSystematicsInternal<T>::GetMetadata(std::string& key)
{
  if(metadata.count(key) == 0)
    return "";
  return metadata.at(key);
}

template<class T>
std::string ValueWithSystematicsInternal<T>::GetMetadata(const std::string& key) const
{
  if(metadata.count(key) == 0)
    return "";
  return metadata.at(key);
}

template<class T>
T& ValueWithSystematicsInternal<T>::Systematic(const std::string& name)
{
  if(systematics.count(name) == 0)
  {
    if(isLocked)
      throw AnalyserException("Unable to add systematic \""+name+"\" after locking.");
    systematics[name] = value;
  }

  return systematics[name];
}

template<class T>
const T& ValueWithSystematicsInternal<T>::Systematic(const std::string& name) const
{
  if(systematics.count(name) == 0)
  {
    if(isLocked)
      throw AnalyserException("Unable to add systematic \""+name+"\" after locking.");
    systematics[name] = value;
  }

  return systematics[name];
}

template<class T>
T& ValueWithSystematicsInternal<T>::operator()(const std::string& name)
{
  return Systematic(name);
}

template<class T>
ValueWithSystematicsInternal<T>::operator T () const
{
  return value;
}

// Specialized method for the bool type, where the return value is the logical or of all systematics
template<>
ValueWithSystematicsInternal<bool>::operator bool () const
{
  bool retVal = value;
  for(auto& kv: systematics)
    retVal = retVal || kv.second;
  return retVal;
}

template<class T>
ValueWithSystematicsInternal<T>& ValueWithSystematicsInternal<T>::operator=(const T& val)
{
  value = val;
  
  for(auto& kv: systematics)
  {
    kv.second = val;
  }

  return *this;
}

template<class T>
ValueWithSystematicsInternal<T>& ValueWithSystematicsInternal<T>::operator=(const ValueWithSystematicsInternal<T>& val)
{
  if(this == &val) // Check for self assignment
    return *this;

  value = val.value;

  if(isLocked)
  {
    for(auto& kv: systematics)
      if(val.systematics.count(kv.first) == 0)
        kv.second = value;
  }
  else
    systematics.clear();
  
  for(auto& kv: val.systematics)
    Systematic(kv.first) = kv.second;

  return *this;
}

template<class T>
ValueWithSystematicsInternal<T>& ValueWithSystematicsInternal<T>::operator+=(const T& val)
{
  static_assert(!(std::is_same<T, bool>::value), "You can not use addition operators with booleans");

  for(auto& kv: systematics)
  {
    kv.second += val;
  }

  value += val;

  return *this;
}

template<class T>
ValueWithSystematicsInternal<T>& ValueWithSystematicsInternal<T>::operator+=(const ValueWithSystematicsInternal<T>& val)
{
  static_assert(!(std::is_same<T, bool>::value), "You can not use addition operators with booleans");

  for(auto& kv: systematics)
    if(val.systematics.count(kv.first) == 0)
      kv.second += val.value;

  for(auto& kv: val.systematics)
    Systematic(kv.first) += kv.second;

  value += val.value;

  return *this;
}

template<class T>
ValueWithSystematicsInternal<T>& ValueWithSystematicsInternal<T>::operator-=(const T& val)
{
  static_assert(!(std::is_same<T, bool>::value), "You can not use subtraction operators with booleans");

  for(auto& kv: systematics)
  {
    kv.second -= val;
  }

  value -= val;

  return *this;
}

template<class T>
ValueWithSystematicsInternal<T>& ValueWithSystematicsInternal<T>::operator-=(const ValueWithSystematicsInternal<T>& val)
{
  static_assert(!(std::is_same<T, bool>::value), "You can not use subtraction operators with booleans");

  for(auto& kv: systematics)
    if(val.systematics.count(kv.first) == 0)
      kv.second -= val.value;
  
  for(auto& kv: val.systematics)
    Systematic(kv.first) -= kv.second;

  value -= val.value;

  return *this;
}

template<class T>
ValueWithSystematicsInternal<T>& ValueWithSystematicsInternal<T>::operator*=(const T& val)
{
  static_assert(!(std::is_same<T, bool>::value), "You can not use multiplication operators with booleans");

  for(auto& kv: systematics)
  {
    kv.second *= val;
  }

  value *= val;

  return *this;
}

template<class T>
ValueWithSystematicsInternal<T>& ValueWithSystematicsInternal<T>::operator*=(const ValueWithSystematicsInternal<T>& val)
{
  static_assert(!(std::is_same<T, bool>::value), "You can not use multiplication operators with booleans");

  for(auto& kv: systematics)
    if(val.systematics.count(kv.first) == 0)
      kv.second *= val.value;
  
  for(auto& kv: val.systematics)
    Systematic(kv.first) *= kv.second;

  value *= val.value;

  return *this;
}

template<class T>
ValueWithSystematicsInternal<T>& ValueWithSystematicsInternal<T>::operator/=(const T& val)
{
  static_assert(!(std::is_same<T, bool>::value), "You can not use division operators with booleans");

  for(auto& kv: systematics)
  {
    kv.second /= val;
  }

  value /= val;

  return *this;
}

template<class T>
ValueWithSystematicsInternal<T>& ValueWithSystematicsInternal<T>::operator/=(const ValueWithSystematicsInternal<T>& val)
{
  static_assert(!(std::is_same<T, bool>::value), "You can not use division operators with booleans");

  for(auto& kv: systematics)
    if(val.systematics.count(kv.first) == 0)
      kv.second /= val.value;
  
  for(auto& kv: val.systematics)
    Systematic(kv.first) /= kv.second;

  value /= val.value;

  return *this;
}

template<class T>
const ValueWithSystematicsInternal<T> ValueWithSystematicsInternal<T>::operator+(const T& val) const
{
  return ValueWithSystematicsInternal<T>(*this) += val;
}

template<class T>
const ValueWithSystematicsInternal<T> ValueWithSystematicsInternal<T>::operator+(const ValueWithSystematicsInternal<T>& val) const
{
  return ValueWithSystematicsInternal<T>(*this) += val;
}

template<class T>
const ValueWithSystematicsInternal<T> ValueWithSystematicsInternal<T>::operator-(const T& val) const
{
  return ValueWithSystematicsInternal<T>(*this) -= val;
}

template<class T>
const ValueWithSystematicsInternal<T> ValueWithSystematicsInternal<T>::operator-(const ValueWithSystematicsInternal<T>& val) const
{
  return ValueWithSystematicsInternal<T>(*this) -= val;
}

template<class T>
const ValueWithSystematicsInternal<T> ValueWithSystematicsInternal<T>::operator*(const T& val) const
{
  return ValueWithSystematicsInternal<T>(*this) *= val;
}

template<class T>
const ValueWithSystematicsInternal<T> ValueWithSystematicsInternal<T>::operator*(const ValueWithSystematicsInternal<T>& val) const
{
  return ValueWithSystematicsInternal<T>(*this) *= val;
}

template<class T>
const ValueWithSystematicsInternal<T> ValueWithSystematicsInternal<T>::operator/(const T& val) const
{
  return ValueWithSystematicsInternal<T>(*this) /= val;
}

template<class T>
const ValueWithSystematicsInternal<T> ValueWithSystematicsInternal<T>::operator/(const ValueWithSystematicsInternal<T>& val) const
{
  return ValueWithSystematicsInternal<T>(*this) /= val;
}

template<class T>
const ValueWithSystematicsInternal<T> ValueWithSystematicsInternal<T>::operator-() const
{
  ValueWithSystematicsInternal<T> retVal(*this);

  retVal.value = -retVal.value;
  for(auto& kv: retVal.systematics)
    kv.second = -kv.second;

  return retVal;
}

template<class T>
ValueWithSystematicsInternal<T>& ValueWithSystematicsInternal<T>::operator++()
{
  static_assert(!(std::is_same<T, bool>::value), "You can not use increment operators with booleans");
  
  ++value;
  for(auto& kv: systematics)
    ++kv.second;
   
  return *this;
}

template<class T>
ValueWithSystematicsInternal<T> ValueWithSystematicsInternal<T>::operator++(int)
{
  static_assert(!(std::is_same<T, bool>::value), "You can not use increment operators with booleans");

  ValueWithSystematicsInternal<T> retVal(*this);  

  ++value;
  for(auto& kv: systematics)
    ++kv.second;
   
  return retVal;
}

template<class T>
ValueWithSystematicsInternal<T>& ValueWithSystematicsInternal<T>::operator--()
{
  static_assert(!(std::is_same<T, bool>::value), "You can not use increment operators with booleans");
  
  --value;
  for(auto& kv: systematics)
    --kv.second;
   
  return *this;
}

template<class T>
ValueWithSystematicsInternal<T> ValueWithSystematicsInternal<T>::operator--(int)
{
  static_assert(!(std::is_same<T, bool>::value), "You can not use increment operators with booleans");
  
  ValueWithSystematicsInternal<T> retVal(*this);  

  --value;
  for(auto& kv: systematics)
    --kv.second;
   
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<bool> ValueWithSystematicsInternal<T>::operator==(const T& val) const
{
  ValueWithSystematicsInternal<bool> retVal(value == val);
  
  for(auto& kv: systematics)
    retVal(kv.first) = (kv.second == val);
  
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<bool> ValueWithSystematicsInternal<T>::operator==(const ValueWithSystematicsInternal<T>& val) const
{
  ValueWithSystematicsInternal<bool> retVal(value == val.value);
  
  for(auto& kv: systematics)
  {
    if(val.systematics.count(kv.first) == 0)
      retVal.systematics[kv.first] = (kv.second == val.value);
  }
  
  for(auto& kv: val.systematics)
  {
    if(systematics.count(kv.first) == 0)
      retVal.systematics[kv.first] = (value == kv.second);
    else
      retVal.systematics[kv.first] = (systematics.at(kv.first) == kv.second);
  }
  
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<bool> ValueWithSystematicsInternal<T>::operator!=(const T& val) const
{
  ValueWithSystematicsInternal<bool> retVal(value != val);
  
  for(auto& kv: systematics)
    retVal.systematics[kv.first] = (kv.second != val);
  
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<bool> ValueWithSystematicsInternal<T>::operator!=(const ValueWithSystematicsInternal<T>& val) const
{
  ValueWithSystematicsInternal<bool> retVal(value != val.value);
  
  for(auto& kv: systematics)
  {
    if(val.systematics.count(kv.first) == 0)
      retVal.systematics[kv.first] = (kv.second != val.value);
  }
  
  for(auto& kv: val.systematics)
  {
    if(systematics.count(kv.first) == 0)
      retVal.systematics[kv.first] = (value != kv.second);
    else
      retVal.systematics[kv.first] = (systematics.at(kv.first) != kv.second);
  }
  
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<bool> ValueWithSystematicsInternal<T>::operator> (const T& val) const
{
  static_assert(!(std::is_same<T, bool>::value), "You can not test > with booleans");
  ValueWithSystematicsInternal<bool> retVal(value > val);
  
  for(auto& kv: systematics)
    retVal(kv.first) = (kv.second > val);
  
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<bool> ValueWithSystematicsInternal<T>::operator> (const ValueWithSystematicsInternal<T>& val) const
{
  static_assert(!(std::is_same<T, bool>::value), "You can not test > with booleans");
  ValueWithSystematicsInternal<bool> retVal(value > val.value);
  
  for(auto& kv: systematics)
  {
    if(val.systematics.count(kv.first) == 0)
      retVal.systematics[kv.first] = (kv.second > val.value);
  }
  
  for(auto& kv: val.systematics)
  {
    if(systematics.count(kv.first) == 0)
      retVal.systematics[kv.first] = (value > kv.second);
    else
      retVal.systematics[kv.first] = (systematics.at(kv.first) > kv.second);
  }
  
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<bool> ValueWithSystematicsInternal<T>::operator< (const T& val) const
{
  static_assert(!(std::is_same<T, bool>::value), "You can not test < with booleans");
  ValueWithSystematicsInternal<bool> retVal(value < val);
  
  for(auto& kv: systematics)
    retVal.systematics[kv.first] = (kv.second < val);
  
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<bool> ValueWithSystematicsInternal<T>::operator< (const ValueWithSystematicsInternal<T>& val) const
{
  static_assert(!(std::is_same<T, bool>::value), "You can not test < with booleans");
  ValueWithSystematicsInternal<bool> retVal(value < val.value);
  
  for(auto& kv: systematics)
  {
    if(val.systematics.count(kv.first) == 0)
      retVal.systematics[kv.first] = (kv.second < val.value);
  }
  
  for(auto& kv: val.systematics)
  {
    if(systematics.count(kv.first) == 0)
      retVal.systematics[kv.first] = (value < kv.second);
    else
      retVal.systematics[kv.first] = (systematics.at(kv.first) < kv.second);
  }
  
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<bool> ValueWithSystematicsInternal<T>::operator>=(const T& val) const
{
  static_assert(!(std::is_same<T, bool>::value), "You can not test >= with booleans");
  ValueWithSystematicsInternal<bool> retVal(value >= val);
  
  for(auto& kv: systematics)
    retVal.systematics[kv.first] = (kv.second >= val);
  
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<bool> ValueWithSystematicsInternal<T>::operator>=(const ValueWithSystematicsInternal<T>& val) const
{
  static_assert(!(std::is_same<T, bool>::value), "You can not test >= with booleans");
  ValueWithSystematicsInternal<bool> retVal(value >= val.value);
  
  for(auto& kv: systematics)
  {
    if(val.systematics.count(kv.first) == 0)
      retVal.systematics[kv.first] = (kv.second >= val.value);
  }
  
  for(auto& kv: val.systematics)
  {
    if(systematics.count(kv.first) == 0)
      retVal.systematics[kv.first] = (value >= kv.second);
    else
      retVal.systematics[kv.first] = (systematics.at(kv.first) >= kv.second);
  }
  
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<bool> ValueWithSystematicsInternal<T>::operator<=(const T& val) const
{
  static_assert(!(std::is_same<T, bool>::value), "You can not test <= with booleans");
  ValueWithSystematicsInternal<bool> retVal(value <= val);
  
  for(auto& kv: systematics)
    retVal.systematics[kv.first] = (kv.second <= val);
  
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<bool> ValueWithSystematicsInternal<T>::operator<=(const ValueWithSystematicsInternal<T>& val) const
{
  static_assert(!(std::is_same<T, bool>::value), "You can not test <= with booleans");
  ValueWithSystematicsInternal<bool> retVal(value <= val.value);
  
  for(auto& kv: systematics)
  {
    if(val.systematics.count(kv.first) == 0)
      retVal(kv.first) = (kv.second <= val.value);
  }
  
  for(auto& kv: val.systematics)
  {
    if(systematics.count(kv.first) == 0)
      retVal(kv.first) = (value <= kv.second);
    else
      retVal(kv.first) = (systematics.at(kv.first) <= kv.second);
  }
  
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<T> ValueWithSystematicsInternal<T>::operator! () const
{
  static_assert(std::is_same<T, bool>::value, "You can only negate with booleans");
  
  ValueWithSystematicsInternal<T> retVal(!value);
  
  for(auto& kv: systematics)
    retVal.systematics[kv.first] = !kv.second;
  
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<T> ValueWithSystematicsInternal<T>::operator&&(const ValueWithSystematicsInternal<T>& val) const
{
  static_assert(std::is_same<T, bool>::value, "You can only perform a logical AND with booleans");
  
  ValueWithSystematicsInternal<T> retVal(value && val.value);
  
  for(auto& kv: systematics)
    if(val.systematics.count(kv.first) == 0)
      retVal.systematics[kv.first] = kv.second && val.value;

  for(auto& kv: val.systematics)
  {
    if(systematics.count(kv.first) == 0)
      retVal.systematics[kv.first] = value && kv.second;
    else
      retVal.systematics[kv.first] = systematics.at(kv.first) && kv.second;
  }
  
  return retVal;
}

template<class T>
const ValueWithSystematicsInternal<T> ValueWithSystematicsInternal<T>::operator||(const ValueWithSystematicsInternal<T>& val) const
{
  static_assert(std::is_same<T, bool>::value, "You can only perform a logical OR with booleans");
  
  ValueWithSystematicsInternal<T> retVal(value || val.value);
  
  for(auto& kv: systematics)
    if(val.systematics.count(kv.first) == 0)
      retVal.systematics[kv.first] = kv.second || val.value;

  for(auto& kv: val.systematics)
  {
    if(systematics.count(kv.first) == 0)
      retVal.systematics[kv.first] = value || kv.second;
    else
      retVal.systematics[kv.first] = systematics.at(kv.first) || kv.second;
  }
  
  return retVal;
}

template<class T>
std::vector<std::string> ValueWithSystematicsInternal<T>::Systematics() const
{
  std::vector<std::string> retVal;
  
  for(auto& kv: systematics)
    retVal.push_back(kv.first);
  
  return retVal;
}

template<class T>
T& ValueWithSystematicsInternal<T>::GetSystematicOrValue(const std::string& name)
{
  if(systematics.count(name) != 0)
    return systematics[name];
  return value;
}

template<class T>
const T& ValueWithSystematicsInternal<T>::GetSystematicOrValue(const std::string& name) const
{
  if(systematics.count(name) != 0)
    return systematics[name];
  return value;
}

template<class T, typename = void>
class ValueWithSystematics: public ValueWithSystematicsInternal<T>
{
public:
//  using ValueWithSystematicsInternal<T>::ValueWithSystematicsInternal; //Why doesn't this one work?
  ValueWithSystematics(T val = T(0)): ValueWithSystematicsInternal<T>(val) {};
  ValueWithSystematics(const ValueWithSystematics<T>& val): ValueWithSystematicsInternal<T>(val) {}; // Copy constructor
  ValueWithSystematics(const ValueWithSystematicsInternal<T>& val): ValueWithSystematicsInternal<T>(val) {}; // Copy constructor
  
/*  // --------- Assignment operators ---------
  ValueWithSystematics<T>& operator= (const T& val) {return static_cast<ValueWithSystematics<T>(ValueWithSystematicsInternal::operator=(val))>;};
  ValueWithSystematics<T>& operator= (const ValueWithSystematicsInternal<T>& val)  {return static_cast<ValueWithSystematics<T>(ValueWithSystematicsInternal::operator=(val))>;};
  //       Compound Assignment Operators
  ValueWithSystematics<T>& operator+=(const T& val)  {return static_cast<ValueWithSystematics<T>(ValueWithSystematicsInternal::operator+=(val))>;};
  ValueWithSystematics<T>& operator+=(const ValueWithSystematicsInternal<T>& val) {return static_cast<ValueWithSystematics<T>(ValueWithSystematicsInternal::operator+=(val))>;};
  ValueWithSystematics<T>& operator-=(const T& val) {return static_cast<ValueWithSystematics<T>(ValueWithSystematicsInternal::operator-=(val))>;};
  ValueWithSystematics<T>& operator-=(const ValueWithSystematicsInternal<T>& val) {return static_cast<ValueWithSystematics<T>(ValueWithSystematicsInternal::operator-=(val))>;};
  ValueWithSystematics<T>& operator*=(const T& val) {return static_cast<ValueWithSystematics<T>(ValueWithSystematicsInternal::operator*=(val))>;};
  ValueWithSystematics<T>& operator*=(const ValueWithSystematicsInternal<T>& val) {return static_cast<ValueWithSystematics<T>(ValueWithSystematicsInternal::operator*=(val))>;};
  ValueWithSystematics<T>& operator/=(const T& val) {return static_cast<ValueWithSystematics<T>(ValueWithSystematicsInternal::operator/=(val))>;};
  ValueWithSystematics<T>& operator/=(const ValueWithSystematicsInternal<T>& val) {return static_cast<ValueWithSystematics<T>(ValueWithSystematicsInternal::operator/=(val))>;};
  // --------- Arithmetic operators ---------
  const ValueWithSystematics<T> operator+(const T& val) const;
  const ValueWithSystematics<T> operator+(const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematics<T> operator-(const T& val) const;
  const ValueWithSystematics<T> operator-(const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematics<T> operator*(const T& val) const;
  const ValueWithSystematics<T> operator*(const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematics<T> operator/(const T& val) const;
  const ValueWithSystematics<T> operator/(const ValueWithSystematicsInternal<T>& val) const;
  // --------- Comparison operators ---------
  const ValueWithSystematics<bool> operator==(const T& val) const;
  const ValueWithSystematics<bool> operator==(const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematics<bool> operator!=(const T& val) const;
  const ValueWithSystematics<bool> operator!=(const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematics<bool> operator> (const T& val) const;
  const ValueWithSystematics<bool> operator> (const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematics<bool> operator< (const T& val) const;
  const ValueWithSystematics<bool> operator< (const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematics<bool> operator>=(const T& val) const;
  const ValueWithSystematics<bool> operator>=(const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematics<bool> operator<=(const T& val) const;
  const ValueWithSystematics<bool> operator<=(const ValueWithSystematicsInternal<T>& val) const;
  // ---------  Logical  operators  ---------
  const ValueWithSystematics<T> operator! () const;
  const ValueWithSystematics<T> operator&&(const ValueWithSystematicsInternal<T>& val) const;
  const ValueWithSystematics<T> operator||(const ValueWithSystematicsInternal<T>& val) const;
  // ---------   Unary  operators   ---------
  const ValueWithSystematics<T> operator-() const;
  ValueWithSystematics<T>& operator++();
  ValueWithSystematics<T> operator++(int);
  ValueWithSystematics<T>& operator--();
  ValueWithSystematics<T> operator--(int);// */

private:
protected:
};

template<>
class ValueWithSystematics<double>;

template<class T>
class ValueWithSystematics<std::vector<T>>: public ValueWithSystematicsInternal<std::vector<T>>
{
public:
  ValueWithSystematics(): ValueWithSystematicsInternal<std::vector<T>>(std::vector<T>(0)) {};
  ValueWithSystematics(std::vector<T> val): ValueWithSystematicsInternal<std::vector<T>>(val) {};
  ValueWithSystematics(const ValueWithSystematics<std::vector<T>>& val): ValueWithSystematicsInternal<std::vector<T>>(val) {}; // Copy constructor
  ValueWithSystematics(const ValueWithSystematicsInternal<std::vector<T>>& val): ValueWithSystematicsInternal<std::vector<T>>(val) {}; // Copy constructor

  ValueWithSystematics<int> size() const;

private:
protected:
  using ValueWithSystematicsInternal<std::vector<T>>::systematics;
  using ValueWithSystematicsInternal<std::vector<T>>::value;
};

template<class T>
class ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>: public ValueWithSystematicsInternal<T>
{
public:
//  using ValueWithSystematicsInternal<T>::ValueWithSystematicsInternal; //Why doesn't this one work?
  ValueWithSystematics(T val = T(0, 0, 0, 0)): ValueWithSystematicsInternal<T>(val) {};
  ValueWithSystematics(const ValueWithSystematics<T>& val): ValueWithSystematicsInternal<T>(val) {}; // Copy constructor
  ValueWithSystematics(const ValueWithSystematicsInternal<T>& val): ValueWithSystematicsInternal<T>(val) {}; // Copy constructor

  ValueWithSystematics<T>& operator*=(const double& val);
  ValueWithSystematics<T>& operator*=(const ValueWithSystematics<double>& val);
  const ValueWithSystematics<T> operator*(const double& val) const {return ValueWithSystematics<T>(*this) *= val;};
  const ValueWithSystematics<T> operator*(const ValueWithSystematics<double>& val) const {return ValueWithSystematics<T>(*this) *= val;};
  
  ValueWithSystematics<double> Pt() const;
  ValueWithSystematics<double> Phi() const;
  ValueWithSystematics<double> Angle(const ValueWithSystematics<T>& other) const;
  ValueWithSystematics<double> DeltaPhi(const ValueWithSystematics<T>& other) const;
  ValueWithSystematics<double> DeltaR(const ValueWithSystematics<T>& other) const;
  template<class U>
  ValueWithSystematics<double> MinDeltaPhi(const ValueWithSystematics<std::vector<U>>& other) const;
  ValueWithSystematics<double> CosTheta() const;
  ValueWithSystematics<TVector3> BoostVector() const;
  ValueWithSystematics<T>& Boost(const ValueWithSystematics<TVector3>& boostVec);
  ValueWithSystematics<TRotation> RotateTozz() const;
  ValueWithSystematics<T>& Transform(const ValueWithSystematics<TRotation>& transformation);
  ValueWithSystematics<double> M() const;
  
  using ValueWithSystematicsInternal<T>::Systematic;
  using ValueWithSystematicsInternal<T>::GetSystematicOrValue;

private:
protected:
  using ValueWithSystematicsInternal<T>::systematics;
  using ValueWithSystematicsInternal<T>::value;
};

template<class T>
class ValueWithSystematics<T, typename std::enable_if<std::is_base_of<LorentzVectorF, T>::value>::type>: public ValueWithSystematicsInternal<T>
{
public:
//  using ValueWithSystematicsInternal<T>::ValueWithSystematicsInternal; //Why doesn't this one work?
  ValueWithSystematics(T val = T()): ValueWithSystematicsInternal<T>(val) {};
  ValueWithSystematics(const ValueWithSystematics<T>& val): ValueWithSystematicsInternal<T>(val) {}; // Copy constructor
  ValueWithSystematics(const ValueWithSystematicsInternal<T>& val): ValueWithSystematicsInternal<T>(val) {}; // Copy constructor
  
  ValueWithSystematics<double> Pt() const;
  ValueWithSystematics<double> Phi() const;
  template<class U>
  ValueWithSystematics<double> DeltaR(const ValueWithSystematics<U>& other) const;
  ValueWithSystematics<TLorentzVector> ToTLorentzVector() const;

private:
protected:
  using ValueWithSystematicsInternal<T>::systematics;
  using ValueWithSystematicsInternal<T>::value;
};

template<>
class ValueWithSystematics<double>: public ValueWithSystematicsInternal<double>
{
public:
//  using ValueWithSystematicsInternal<T>::ValueWithSystematicsInternal; //Why doesn't this one work?
  ValueWithSystematics(double val = 0): ValueWithSystematicsInternal<double>(val) {};
  ValueWithSystematics(const ValueWithSystematics<double>& val): ValueWithSystematicsInternal<double>(val) {}; // Copy constructor
  ValueWithSystematics(const ValueWithSystematicsInternal<double>& val): ValueWithSystematicsInternal<double>(val) {}; // Copy constructor

  ValueWithSystematics<double> Cos() const;
  ValueWithSystematics<double> Sqrt() const;
  ValueWithSystematics<double> abs() const;
  
  template<class U>
  friend ValueWithSystematics<U>& ValueWithSystematics<U, typename std::enable_if<std::is_base_of<TLorentzVector, U>::value>::type>::operator*=(const ValueWithSystematics<double>& val);

private:
protected:
  using ValueWithSystematicsInternal<double>::systematics;
  using ValueWithSystematicsInternal<double>::value;
};

ValueWithSystematics<double> ValueWithSystematics<double>::Cos() const
{
  ValueWithSystematics<double> retVal = cos(value);
  
  for(auto& kv: systematics)
    retVal(kv.first) = cos(kv.second);
  
  return retVal;
}

ValueWithSystematics<double> ValueWithSystematics<double>::Sqrt() const
{
  ValueWithSystematics<double> retVal = sqrt(value);
  
  for(auto& kv: systematics)
    retVal(kv.first) = sqrt(kv.second);
  
  return retVal;
}

ValueWithSystematics<double> ValueWithSystematics<double>::abs() const
{
  ValueWithSystematics<double> retVal = std::abs(value);
  
  for(auto& kv: systematics)
    retVal(kv.first) = std::abs(kv.second);
  
  return retVal;
}

template<class T>
ValueWithSystematics<int> ValueWithSystematics<std::vector<T>>::size() const
{
  ValueWithSystematics<int> retVal(value.size());

  for(auto& kv: systematics)
    retVal(kv.first) = kv.second.size();

  return retVal;
}

template<class T>
ValueWithSystematics<T>& ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>::operator*=(const double& val)
{
  for(auto& kv: systematics)
  {
    kv.second *= val;
  }

  value *= val;

  return *this;
}

template<class T>
ValueWithSystematics<T>& ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>::operator*=(const ValueWithSystematics<double>& val)
{
  for(auto& kv: systematics)
    if(val.systematics.count(kv.first) == 0)
      kv.second *= val.value;
  
  for(auto& kv: val.systematics)
    Systematic(kv.first) *= kv.second;

  value *= val.value;

  return *this;
}

template<class T>
ValueWithSystematics<double> ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>::Pt() const
{
  ValueWithSystematics<double> retVal = value.Pt();
  
  for(auto& kv: systematics)
    retVal(kv.first) = kv.second.Pt();
  
  return retVal;
}

template<class T>
ValueWithSystematics<double> ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>::Phi() const
{
  ValueWithSystematics<double> retVal = value.Phi();
  
  for(auto& kv: systematics)
    retVal(kv.first) = kv.second.Phi();
  
  return retVal;
}

template<class T>
ValueWithSystematics<double> ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>::Angle(const ValueWithSystematics<T>& other) const
{
  ValueWithSystematics<double> retVal = value.Angle(other.value.Vect());
  
  for(auto& kv: systematics)
    if(other.systematics.count(kv.first) == 0)
      retVal(kv.first) = kv.second.Angle(other.value.Vect());

  for(auto& kv: other.systematics)
  {
    if(systematics.count(kv.first) == 0)
      retVal(kv.first) = value.Angle(kv.second.Vect());
    else
      retVal(kv.first) = systematics.at(kv.first).Angle(kv.second.Vect());
  }
  
  return retVal;
}

template<class T>
ValueWithSystematics<double> ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>::DeltaPhi(const ValueWithSystematics<T>& other) const
{
  ValueWithSystematics<double> retVal = value.DeltaPhi(other.value);
  
  for(auto& kv: systematics)
    if(other.systematics.count(kv.first) == 0)
      retVal(kv.first) = kv.second.DeltaPhi(other.value);

  for(auto& kv: other.systematics)
  {
    if(systematics.count(kv.first) == 0)
      retVal(kv.first) = value.DeltaPhi(kv.second);
    else
      retVal(kv.first) = systematics.at(kv.first).DeltaPhi(kv.second);
  }
  
  return retVal;
}

template<class T>
ValueWithSystematics<double> ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>::DeltaR(const ValueWithSystematics<T>& other) const
{
  ValueWithSystematics<double> retVal = value.DeltaR(other.value);
  
  for(auto& kv: systematics)
    if(other.systematics.count(kv.first) == 0)
      retVal(kv.first) = kv.second.DeltaR(other.value);

  for(auto& kv: other.systematics)
  {
    if(systematics.count(kv.first) == 0)
      retVal(kv.first) = value.DeltaR(kv.second);
    else
      retVal(kv.first) = systematics.at(kv.first).DeltaR(kv.second);
  }
  
  return retVal;
}

template<class T>
template<class U>
ValueWithSystematics<double> ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>::MinDeltaPhi(const ValueWithSystematics<std::vector<U>>& other) const
{
  ValueWithSystematics<double> retVal;
  
  std::vector<std::string> tmpLoop;
  tmpLoop.push_back("Value");
  for(auto& kv: systematics)
    tmpLoop.push_back(kv.first);
  for(auto& syst: other.Systematics())
    if(std::find(tmpLoop.begin(), tmpLoop.end(), syst) == tmpLoop.end())
      tmpLoop.push_back(syst);

  for(auto& val: tmpLoop)
  {
    if(val != "Value")
      retVal(val);

    auto& retVal_ = retVal.GetSystematicOrValue(val);
    auto& vec = GetSystematicOrValue(val);
    auto& list = other.GetSystematicOrValue(val);
    
    retVal_ = 10;
    for(auto& entry: list)
    {
      TLorentzVector vec2(entry.Px(), entry.Py(), entry.Pz(), entry.E());
      
      double temp = vec.DeltaPhi(vec2);
      if(temp < retVal_)
        retVal_ = temp;
    }
  }
  
  return retVal;
}

template<class T>
ValueWithSystematics<double> ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>::CosTheta() const
{
  ValueWithSystematics<double> retVal = value.CosTheta();
  
  for(auto& kv: systematics)
    retVal(kv.first) = kv.second.CosTheta();
  
  return retVal;
}

template<class T>
ValueWithSystematics<TVector3> ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>::BoostVector() const
{
  ValueWithSystematics<TVector3> retVal = value.BoostVector();
  
  for(auto& kv: systematics)
    retVal(kv.first) = kv.second.BoostVector();
  
  return retVal;
}

template<class T>
ValueWithSystematics<T>& ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>::Boost(const ValueWithSystematics<TVector3>& boostVec)
{
  for(auto& kv: systematics)
  {
//    std::cerr << "      Doing syst:" << kv.first << "\n";
    auto tmpVec = boostVec.Systematics();
//    std::cerr << "      BoostVec systs:\n";
//    for(auto& temp: tmpVec)
//    {
//      std::cerr << "        " << temp << "\n";
//    }
//    std::cerr << "      Doing loop\n";
    if(std::find(tmpVec.begin(), tmpVec.end(), kv.first) == tmpVec.end())
      kv.second.Boost(boostVec.Value());
//    std::cerr << "      Loop done\n";
  }
//  std::cerr << "Between loops\n";

  for(auto& syst: boostVec.Systematics())
  {
//    std::cerr << "Creating holder for systematic " << syst << " if it doesn't exist yet\n";
    if(systematics.count(syst) == 0)
      systematics[syst] = value;
//    std::cerr << "Holder for " << syst << " ready to be used\nBoosting now\n";
    systematics[syst].Boost(boostVec.Systematic(syst));
//    std::cerr << "Done boosting\n";
  }

  value.Boost(boostVec.Value());

//  std::cerr << "End of function" << std::endl;

  return *this;
}

template<class T>
ValueWithSystematics<TRotation> ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>::RotateTozz() const
{
  TRotation tmpRot;
  ValueWithSystematics<TRotation> retVal(tmpRot);
  std::vector<std::string> tmpLoop;
  tmpLoop.push_back("Value");

  for(auto& kv: systematics)
  {
    retVal(kv.first);
    tmpLoop.push_back(kv.first);
  }

  for(auto& val: tmpLoop)
  {
    auto& SQA = GetSystematicOrValue(val);

    TVector3 newZAxis = SQA.Vect().Unit();
    TVector3 targetZaxis(0, 0, 1);
    TVector3 rotAxis = targetZaxis.Cross(newZAxis);
    double rotAngle = targetZaxis.Angle(newZAxis);

    retVal.GetSystematicOrValue(val).Rotate(-rotAngle, rotAxis);
  }

  return retVal;
}

template<class T>
ValueWithSystematics<T>& ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>::Transform(const ValueWithSystematics<TRotation>& transformation)
{
  for(auto& kv: systematics)
  {
    auto tmp = transformation.Systematics();
    if(std::find(tmp.begin(), tmp.end(), kv.first) == tmp.end())
      kv.second.Transform(transformation.Value());
  }

  for(auto& syst: transformation.Systematics())
  {
    if(systematics.count(syst) == 0)
      systematics[syst] = value;
    systematics[syst].Transform(transformation.Systematic(syst));
  }

  value.Transform(transformation.Value());
  
  return *this;
}

template<class T>
ValueWithSystematics<double> ValueWithSystematics<T, typename std::enable_if<std::is_base_of<TLorentzVector, T>::value>::type>::M() const
{
  ValueWithSystematics<double> retVal = value.M();
  
  for(auto& kv: systematics)
    retVal(kv.first) = kv.second.M();
  
  return retVal;
}

template<class T>
ValueWithSystematics<double> ValueWithSystematics<T, typename std::enable_if<std::is_base_of<LorentzVectorF, T>::value>::type>::Pt() const
{
  ValueWithSystematics<double> retVal = value.Pt();
  
  for(auto& kv: systematics)
    retVal(kv.first) = kv.second.Pt();
  
  return retVal;
}

template<class T>
ValueWithSystematics<double> ValueWithSystematics<T, typename std::enable_if<std::is_base_of<LorentzVectorF, T>::value>::type>::Phi() const
{
  ValueWithSystematics<double> retVal = value.Phi();
  
  for(auto& kv: systematics)
    retVal(kv.first) = kv.second.Phi();
  
  return retVal;
}

template<class T>
template<class U>
ValueWithSystematics<double> ValueWithSystematics<T, typename std::enable_if<std::is_base_of<LorentzVectorF, T>::value>::type>::DeltaR(const ValueWithSystematics<U>& other) const
{
  ValueWithSystematics<double> retVal = deltaR(value, other.Value());
  
  for(auto& kv: systematics)
  {
    auto tmp = other.Systematics();
    if(std::find(tmp.begin(), tmp.end(), kv.first) == tmp.end())
      retVal(kv.first) = deltaR(kv.second, other.Value());
  }

  for(auto& syst: other.Systematics())
  {
    const auto& tmp = other.Systematic(syst);
    if(systematics.count(syst) == 0)
      retVal(syst) = deltaR(value, tmp);
    else
      retVal(syst) = deltaR(systematics.at(syst), tmp);
  }
  
  return retVal;
}

template<class T>
ValueWithSystematics<TLorentzVector> ValueWithSystematics<T, typename std::enable_if<std::is_base_of<LorentzVectorF, T>::value>::type>::ToTLorentzVector() const
{
  ValueWithSystematics<TLorentzVector> retVal = TLorentzVector(value.Px(), value.Py(), value.Pz(), value.E());
  
  for(auto& kv: systematics)
    retVal(kv.first) = TLorentzVector(kv.second.Px(), kv.second.Py(), kv.second.Pz(), kv.second.E());
  
  return retVal;
}

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
  inline ValueWithSystematics<double>& GetDouble(std::string name);
  ValueWithSystematics<int>&    AddInt   (std::string name, int defaultVal);
  inline ValueWithSystematics<int>&    GetInt   (std::string name);
  ValueWithSystematics<bool>&   AddBool  (std::string name, bool defaultVal);
  inline ValueWithSystematics<bool>&   GetBool  (std::string name);
  
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

template<class T>
void EventInfo::OutputValueList(ofstream& file, const ValueWithSystematics<T>& val) const
{
  std::string metadata = val.GetMetadata("eventlist");
  if(metadata == "true")
  {
    const std::string widthStr = val.GetMetadata("eventlistWidth");
    int width = 15;
    if(widthStr != "")
    {
      std::stringstream tmp;
      tmp << widthStr;
      tmp >> width;
      if(width == 0)
        width = 15;
    }

    file << std::setw(width) << val.Value() << " | ";
  }

  return;
}

template<>
void EventInfo::OutputValueList(ofstream& file, const ValueWithSystematics<bool>& val) const
{
  std::string metadata = val.GetMetadata("eventlist");
  if(metadata == "true")
  {
    const std::string widthStr = val.GetMetadata("eventlistWidth");
    int width = 15;
    if(widthStr != "")
    {
      std::stringstream tmp;
      tmp << widthStr;
      tmp >> width;
      if(width == 0)
        width = 15;
    }

    file << std::setw(width) << (val.Value()?("True"):("False")) << " | ";
  }

  return;
}

EventInfo::EventInfo(): isLocked(false)
{
}

void EventInfo::Reset()
{
  if(isLocked)
  {
    for(auto& kv : eventDoubles)
      kv.second.Reset();
    for(auto& kv : eventInts)
      kv.second.Reset();
    for(auto& kv : eventBools)
      kv.second.Reset();
  }
  else
  {
    eventDoubles.clear();
    eventInts.clear();
    eventBools.clear();
  }
}

ValueWithSystematics<double>& EventInfo::AddDouble(std::string name, double defaultVal = 0.0)
{
  if(eventDoubles.count(name) == 0)
  {
    if(isLocked)
      throw AnalyserException("Tried to add more contents after locking the event content");
    eventDoubles[name] = ValueWithSystematics<double>(defaultVal);
    eventDoubles[name].DefaultValue() = defaultVal;
  }
  else
    std::cout << "The variable " << name << " already exists. No action taken." << std::endl;

  return eventDoubles.at(name);
}

inline ValueWithSystematics<double>& EventInfo::GetDouble(std::string name)
{
  if(eventDoubles.count(name) == 0)
    throw AnalyserException("Tried to access non-existing value: "+name);
  return eventDoubles.at(name);
}

ValueWithSystematics<int>&    EventInfo::AddInt   (std::string name, int defaultVal = 0)
{
  if(eventInts.count(name) == 0)
  {
    if(isLocked)
      throw AnalyserException("Tried to add more contents after locking the event content");
    eventInts[name] = ValueWithSystematics<int>(defaultVal);
    eventInts[name].DefaultValue() = defaultVal;
  }
  else
    std::cout << "The variable " << name << " already exists. No action taken." << std::endl;

  return eventInts.at(name);
}

inline ValueWithSystematics<int>&    EventInfo::GetInt   (std::string name)
{
  if(eventInts.count(name) == 0)
    throw AnalyserException("Tried to access non-existing value: "+name);
  return eventInts.at(name);
}

ValueWithSystematics<bool>&   EventInfo::AddBool  (std::string name, bool defaultVal = false)
{
  if(eventBools.count(name) == 0)
  {
    if(isLocked)
      throw AnalyserException("Tried to add more contents after locking the event content");
    eventBools[name] = ValueWithSystematics<bool>(defaultVal);
    eventBools[name].DefaultValue() = defaultVal;
  }
  else
    std::cout << "The variable " << name << " already exists. No action taken." << std::endl;

  return eventBools.at(name);
}

inline ValueWithSystematics<bool>&   EventInfo::GetBool  (std::string name)
{
  if(eventBools.count(name) == 0)
    throw AnalyserException("Tried to access non-existing value: "+name);
  return eventBools.at(name);
}

void EventInfo::SetSummaryTreeBranches(TTree* const tree)
{
  for(auto& kv: eventDoubles)
    AddBranch(tree, kv.second, "d_"+kv.first);
  for(auto& kv: eventInts)
    AddBranch(tree, kv.second, "i_"+kv.first);
  for(auto& kv: eventBools)
    AddBranch(tree, kv.second, "b_"+kv.first);

  return;
}

template<class T>
void EventInfo::AddBranch(TTree* const tree, ValueWithSystematics<T>& val, std::string name)
{
  std::string metadata = val.GetMetadata("eventtree");
  if(metadata == "true" || metadata == "")
  {
    tree->Branch(name.c_str(), &(val.Value()));
    
    for(auto& syst: val.Systematics())
    {
      tree->Branch((name + "_" + syst).c_str(), &(val.Systematic(syst)));
    }
  }

  return;
}

// TODO: make an option for the event list to be outputted as a tsv
void EventInfo::OutputEventListHeader(ofstream& file, const std::vector<std::string>& priority) const
{
  for(auto& entry: priority)
  {
    if(eventDoubles.count(entry) != 0)
      OutputValueListHeader(file, eventDoubles.at(entry), entry);
    if(eventInts.count(entry) != 0)
      OutputValueListHeader(file, eventInts.at(entry), entry);
    if(eventBools.count(entry) != 0)
      OutputValueListHeader(file, eventBools.at(entry), entry);
  }
  
  for(auto& kv: eventDoubles)
    if (std::find(priority.begin(), priority.end(), kv.first) == priority.end())
      OutputValueListHeader(file, kv.second, kv.first);
  for(auto& kv: eventInts)
    if (std::find(priority.begin(), priority.end(), kv.first) == priority.end())
      OutputValueListHeader(file, kv.second, kv.first);
  for(auto& kv: eventBools)
    if (std::find(priority.begin(), priority.end(), kv.first) == priority.end())
      OutputValueListHeader(file, kv.second, kv.first);

  file << "\n";
  return;
}

void EventInfo::OutputEventList(ofstream& file, const std::vector<std::string>& priority) const
{
  for(auto& entry: priority)
  {
    if(eventDoubles.count(entry) != 0)
      OutputValueList(file, eventDoubles.at(entry));
    if(eventInts.count(entry) != 0)
      OutputValueList(file, eventInts.at(entry));
    if(eventBools.count(entry) != 0)
      OutputValueList(file, eventBools.at(entry));
  }
  
  for(auto& kv: eventDoubles)
    if (std::find(priority.begin(), priority.end(), kv.first) == priority.end())
      OutputValueList(file, kv.second);
  for(auto& kv: eventInts)
    if (std::find(priority.begin(), priority.end(), kv.first) == priority.end())
      OutputValueList(file, kv.second);
  for(auto& kv: eventBools)
    if (std::find(priority.begin(), priority.end(), kv.first) == priority.end())
      OutputValueList(file, kv.second);

  file << "\n";
  return;
}

template<class T>
void EventInfo::OutputValueListHeader(ofstream& file, const ValueWithSystematics<T>& val, const std::string& name) const
{
  std::string metadata = val.GetMetadata("eventlist");
  if(metadata == "true")
  {
    const std::string widthStr = val.GetMetadata("eventlistWidth");
    int width = 15;
    if(widthStr != "")
    {
      std::stringstream tmp;
      tmp << widthStr;
      tmp >> width;
      if(width <= 0)
        width = 15;
    }

    file << std::setw(width) << name << " | ";
  }

  return;
}

class Analyser
{
public:
  Analyser(std::string cfgFile);
  virtual ~Analyser();

  virtual void Setup();
  virtual void LoopOverEvents();

  inline void SetEventLimit(size_t val) { limitEvents = val; };
  inline void SetDebugEvent(bool val)   { debugEvent = val; };
  inline void SetSkipEvents(int val)    { skipEvents = val; };

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

Analyser::Analyser(std::string cfgFile): limitEvents(0), debugEvent(false), skipEvents(0), isSetup(false), analyserCout(std::cout.rdbuf()), saveRedirect(false), keepAllEvents(false), mergeBoostedTaus(false)
{
  // Read the cfgFile
  cfgOptions = (edm::readPSetsFrom(cfgFile.c_str())->getParameter<edm::ParameterSet>("runProcess"));
}

Analyser::~Analyser()
{
}

void Analyser::LoadCfgOptions()
{
  isMC               = cfgOptions.getParameter<bool>("isMC");
  crossSection       = cfgOptions.getParameter<double>("xsec");
  crossSection_      = crossSection;
  fileList           = cfgOptions.getParameter<std::vector<std::string>>("input");
  baseDir            = cfgOptions.getParameter<std::string>("dirName");
  outDir             = cfgOptions.getParameter<std::string>("outdir");
  jecDir             = cfgOptions.getParameter<std::string>("jecDir");
  pdfDir             = cfgOptions.getParameter<std::string>("pdfDir");
  runSystematics     = cfgOptions.getParameter<bool>("runSystematics");
  saveSummaryTree    = cfgOptions.getParameter<bool>("saveSummaryTree");
  pileupDistribution = cfgOptions.getParameter<std::vector<double> >("datapileup");

  applyScaleFactors = true;
  debug = false;
  doDDBkg = false;
  outputEventList = false;

  if(cfgOptions.exists("applyScaleFactors"))
    applyScaleFactors = cfgOptions.getParameter<bool>("applyScaleFactors");
  if(cfgOptions.exists("debug"))
    debug             = cfgOptions.getParameter<bool>("debug");
  if(cfgOptions.exists("doDDBkg"))
    doDDBkg           = cfgOptions.getParameter<bool>("doDDBkg");
  if(cfgOptions.exists("outputEventList"))
    outputEventList   = cfgOptions.getParameter<bool>("outputEventList");
    
  if(!isMC && !doDDBkg)
    runSystematics = false;

  if(runSystematics && isMC)
  {
    // From: https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
    crossSection_("xsec_UP");
    crossSection_("xsec_DOWN");
    //crossSection_("PDF_UP");
    //crossSection_("PDF_DOWN");
    //crossSection_("topMass_UP");
    //crossSection_("topMass_DOWN");

    TString turl(fileList[0]);
    if(turl.Contains("MC8TeV_SingleT") || turl.Contains("MC8TeV_TTJetsMassiveBinDecay") || turl.Contains("MC8TeV_TTWJets"))
    {
      crossSection_("PDF_UP");
      crossSection_("PDF_DOWN");
      
      if(turl.Contains("MC8TeV_SingleT_t"))
      {
        if(turl.Contains("MC8TeV_SingleT_tW"))
        {
          crossSection_("xsec_UP")   += 0.3;
          crossSection_("xsec_DOWN") -= 0.3;
          crossSection_("PDF_UP")    += 0.7;
          crossSection_("PDF_DOWN")  -= 0.7;
        }
        else
        {
          crossSection_("xsec_UP")      += 1.64;
          crossSection_("xsec_DOWN")    -= 1.09;
          crossSection_("PDF_UP")       += 1.60;
          crossSection_("PDF_DOWN")     -= 1.60;
//          crossSection_("topMass_UP")   += 0.52;
//          crossSection_("topMass_DOWN") -= 0.52;
        }
      }
      if(turl.Contains("MC8TeV_SingleT_s"))
      {
        crossSection_("xsec_UP")   += 0.07;
        crossSection_("xsec_DOWN") -= 0.07;
        crossSection_("PDF_UP")    += 0.13;
        crossSection_("PDF_DOWN")  -= 0.13;
      }
      
      if(turl.Contains("MC8TeV_SingleTbar_t"))
      {
        if(turl.Contains("MC8TeV_SingleTbar_tW"))
        {
          crossSection_("xsec_UP")   += 0.3;
          crossSection_("xsec_DOWN") -= 0.3;
          crossSection_("PDF_UP")    += 0.7;
          crossSection_("PDF_DOWN")  -= 0.7;
        }
        else
        {
          crossSection_("xsec_UP")      += 0.92;
          crossSection_("xsec_DOWN")    -= 0.59;
          crossSection_("PDF_UP")       += 1.39;
          crossSection_("PDF_DOWN")     -= 1.39;
//          crossSection_("topMass_UP")   += 0.30;
//          crossSection_("topMass_DOWN") -= 0.30;
        }
      }
      if(turl.Contains("MC8TeV_SingleTbar_s"))
      {
        crossSection_("xsec_UP")   += 0.01;
        crossSection_("xsec_DOWN") -= 0.01;
        crossSection_("PDF_UP")    += 0.08;
        crossSection_("PDF_DOWN")  -= 0.08;
      }
      
      if(turl.Contains("MC8TeV_TTJetsMassiveBinDecay"))
      {
        crossSection_("xsec_UP")   += 6.23;
        crossSection_("xsec_DOWN") -= 8.41;
        crossSection_("PDF_UP")    += 11.43;
        crossSection_("PDF_DOWN")  -= 11.43;
      }
      
      if(turl.Contains("MC8TeV_TTWJets"))
      {
        crossSection_("xsec_UP")   += 0.067;
        crossSection_("xsec_DOWN") -= 0.067;
        crossSection_("PDF_UP")    += 0.03;
        crossSection_("PDF_DOWN")  -= 0.03;
      }
    }
    
    if(turl.Contains("MC8TeV_TTZJets"))
    {
      crossSection_("xsec_UP")   += 0.019;
      crossSection_("xsec_DOWN") -= 0.024;
    }
    if(turl.Contains("MC8TeV_WW"))
    {
      if(turl.Contains("MC8TeV_WWWJets"))
      {
        crossSection_("xsec_UP")   *= 1.047;
        crossSection_("xsec_DOWN") *= 0.961;
      }
      else
      {
        if(turl.Contains("MC8TeV_WWZJets"))
        {
          crossSection_("xsec_UP")   *= 1.056;
          crossSection_("xsec_DOWN") *= 0.954;
        }
        else
        {
          crossSection_("xsec_UP")   += 0.2079;
          crossSection_("xsec_DOWN") -= 0.2079;
          crossSection_("PDF_UP")    += 0.2394;
          crossSection_("PDF_DOWN")  -= 0.2394;
        }
      }
    }
    if(turl.Contains("MC8TeV_WZ"))
    {
      if(turl.Contains("MC8TeV_WZZJets"))
      {
        crossSection_("xsec_UP")   *= 1.06;
        crossSection_("xsec_DOWN") *= 0.951;
      }
      else
      {
        crossSection_("xsec_UP")   += 0.02938;
        crossSection_("xsec_DOWN") -= 0.02938;
        crossSection_("PDF_UP")    += 0.03072;
        crossSection_("PDF_DOWN")  -= 0.03072;
      }
    }
    if(turl.Contains("MC8TeV_ZZ"))
    {
      if(turl.Contains("MC8TeV_ZZZJets"))
      {
        crossSection_("xsec_UP")   *= 1.027;
        crossSection_("xsec_DOWN") *= 0.976;
      }
      else
      {
        crossSection_("xsec_UP")   += 0.0099;
        crossSection_("xsec_DOWN") -= 0.0099;
        crossSection_("PDF_UP")    += 0.0099;
        crossSection_("PDF_DOWN")  -= 0.0099;
      }
    }
    if(turl.Contains("MC8TeV_GJets"))
    {
      crossSection_("xsec_UP")   *= 2;
      crossSection_("xsec_DOWN") *= 0.5;
    }
    if(turl.Contains("MC8TeV_QCD"))
    {
      crossSection_("xsec_UP")   *= 2;
      crossSection_("xsec_DOWN") *= 0.5;
    }
    if(turl.Contains("MC8TeV_DY"))
    {
      if(turl.Contains("50toInf"))
      {
        crossSection_("xsec_UP")   += 17.7;
        crossSection_("xsec_DOWN") -= 10.8;
        crossSection_("PDF_UP")    += 116.4;
        crossSection_("PDF_DOWN")  -= 116.4;
      }
      else
      {
        crossSection_("xsec_UP")   += 17.7; // Bah?
        crossSection_("xsec_DOWN") -= 10.8;
        crossSection_("PDF_UP")    += 28;
        crossSection_("PDF_DOWN")  -= 28;
      }
    }
    if(turl.Contains("MC8TeV_WJets") || turl.Contains("MC8TeV_W1Jets") || turl.Contains("MC8TeV_W2Jets") || turl.Contains("MC8TeV_W3Jets") || turl.Contains("MC8TeV_W4Jets"))
    {
      crossSection_("xsec_UP")   += 237;
      crossSection_("xsec_DOWN") -= 119.1;
      crossSection_("PDF_UP")    += 1244.1;
      crossSection_("PDF_DOWN")  -= 1244.1;
    }

    crossSection_.Lock();
  }

  if(debug)
    std::cout << "Finished Analyser::LoadCfgOptions()" << std::endl;

  UserLoadCfgOptions();

  return;
}

void Analyser::Setup()
{
  LoadCfgOptions();

  // Create output directory if it doesn't exist
  gSystem->Exec(("mkdir -p " + outDir).c_str());

  std::string url = fileList[0];
  std::string outFileUrl(gSystem->BaseName(url.c_str()));
  while(outFileUrl.find(".root", 0) != std::string::npos)
    outFileUrl.replace(outFileUrl.find(".root", 0), 5, "");
  outFile = outDir + "/" + outFileUrl + ".root";
  TString turl(url);

  if(saveSummaryTree)
  {
    TDirectory* cwd = gDirectory;

    summaryOutFile = outFile;
    summaryOutFile.replace(summaryOutFile.find(".root", 0), 5, "_summary.root");

    summaryOutTFile = new TFile(summaryOutFile.c_str(), "RECREATE");
    summaryTree = new TTree("Events", "Events");
    summaryTree->SetDirectory(summaryOutTFile);  // This line is probably not needed

    cwd->cd();
  }

  if(outputEventList)
  {
    eventlistOutFile = outFile;
    eventlistOutFile.replace(eventlistOutFile.find(".root", 0), 5, "_eventlist.txt");
    std::cout << "Saving event list to " << eventlistOutFile << std::endl;
    eventListFile.open(eventlistOutFile);
  }
  
  if(isMC && runSystematics)
  {
    // PDF variations
    std::string baseFileName = gSystem->BaseName(url.c_str());
    while(baseFileName.find(".root", 0) != std::string::npos)
      baseFileName.replace(baseFileName.find(".root", 0), 5, "");
    baseFileName = pdfDir + "/" + baseFileName + "_pdf.root";
    
    pdfVariations.push_back(new PDFInfo(baseFileName.c_str(), "CT10.LHgrid"));
    pdfVariations.push_back(new PDFInfo(baseFileName.c_str(), "MSTW2008nlo68cl.LHgrid"));
    pdfVariations.push_back(new PDFInfo(baseFileName.c_str(), "NNPDF23_nlo_as_0119.LHgrid"));
  }
  
  isV0JetsMC = isMC && (turl.Contains("DYJetsToLL_50toInf") || turl.Contains("WJets"));

  if(debug)
    std::cout << "Finished Analyser::Setup()" << std::endl;

  UserSetup();
  EventContentSetup();

  isSetup = true;

  return;
}

void Analyser::LoopOverEvents()
{
  if(!isSetup)
    Setup();

  InitHistograms();

  if(debug)
    std::cout << "Preparing for event loop" << std::endl;
  fwlite::ChainEvent ev(fileList);
  const size_t totalEntries = ev.size();

  // MC normalization to 1/pb
  double nInitEvent = 1.;
  xsecWeight = 1.;
  if(isMC)
  {
    nInitEvent = static_cast<double>(utils::getMergeableCounterValue(fileList, "startCounter"));
    xsecWeight = crossSection_/nInitEvent;
  }

  // Jet Energy Scale and Uncertainties
  jecDir = gSystem->ExpandPathName(jecDir.c_str());
  jesCor = utils::cmssw::getJetCorrector(jecDir, isMC);
  totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt"));

  // Muon Energy Scale and Uncertainties
  muCor = getMuonCorrector(jecDir, fileList[0]);

  // Pileup Weighting: Based on vtx
  LumiWeights = NULL;
  PUNorm[0] = 1;
  PUNorm[1] = 1;
  PUNorm[2] = 1;
  if(isMC)
  {
    std::vector<double> dataPileupDistributionDouble = pileupDistribution;
    std::vector<float> dataPileupDistribution;
    for(unsigned int i = 0; i < dataPileupDistributionDouble.size(); ++i)
      dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);
    std::vector<float> mcPileupDistribution;

    utils::getMCPileupDistribution(ev, dataPileupDistribution.size(), mcPileupDistribution);
    while(mcPileupDistribution.size() < dataPileupDistribution.size()) mcPileupDistribution.push_back(0.0);
    while(mcPileupDistribution.size() > dataPileupDistribution.size()) dataPileupDistribution.push_back(0.0);

    LumiWeights= new edm::LumiReWeighting(mcPileupDistribution, dataPileupDistribution);
    PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution, 0.05);
    utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
  }

  gROOT->cd(); //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
  
  int step = int(totalEntries/50);

  // Redirect stdout and stderr. It can be chosen to redirect to a buffer or to /dev/null
  // use analyserCout instead of cout between the redirects
  std::streambuf *coutbuf = std::cout.rdbuf();
  std::streambuf *cerrbuf = std::cerr.rdbuf();  
  std::ofstream devnull("/dev/null");
  if(saveRedirect)
  {
    std::cout.rdbuf(analyserBuffer.rdbuf());
    std::cerr.rdbuf(analyserBuffer.rdbuf());
  }
  else
  {
    std::cout.rdbuf(devnull.rdbuf());
    std::cerr.rdbuf(devnull.rdbuf());
  }

  analyserCout << "       Progress Bar:0%      20%       40%       60%       80%      100%" << std::endl;
  analyserCout << "Scanning the ntuple:";

  bool doneFirstEvent = false;
  std::vector<std::string> priorityOutput;
  priorityOutput.push_back("RunNo");
  priorityOutput.push_back("LumiNo");
  priorityOutput.push_back("EventNo");
  priorityOutput.push_back("selected");
  priorityOutput.push_back("weight");
  priorityOutput.push_back("PUweight");
  priorityOutput.push_back("xsecweight");
  size_t nEventsOut = 0;
  // Loop on events
  for(size_t iev = 0; iev < totalEntries; ++iev)
  {
    if(iev%step == 0)
      analyserCout << "_" << std::flush;
    if(iev < skipEvents)
      continue;
    
    if(debugEvent)
      analyserCout << "## Event " << iev << std::endl;

    if(doneFirstEvent)
      eventContent.Reset();
    
    ev.to(int(iev));
    
    //Load information/collections from the event
    // Number of vertexes
    auto& nvtx = eventContent.GetInt("nvtx");
    fwlite::Handle<int> nvtxHandle;
    nvtxHandle.getByLabel(ev, "llvvObjectProducersUsed", "nvtx");
    if(nvtxHandle.isValid()) nvtx = *nvtxHandle;
    else continue; // TODO: Maybe remove this?
    
    // Collection of generated particles
    fwlite::Handle<llvvGenEvent> genEventHandle;
    genEventHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!genEventHandle.isValid())
    {
      std::cout << "llvvGenEvent Object NotFound" << std::endl;
      continue;
    }
    genEv = *genEventHandle;
    
    /**** Get LHE comments ****/
    LHEHandle.getByLabel(ev, "source");
    
    // Trigger Bits
    fwlite::Handle<std::vector<bool> > triggerBitsHandle;
    triggerBitsHandle.getByLabel(ev, "llvvObjectProducersUsed", "triggerBits");
    if(!triggerBitsHandle.isValid())
    {
      std::cout << "triggerBits Object NotFound" << std::endl;
      continue;
    }
    triggerBits = *triggerBitsHandle;
    
    // Rest of Gen Particles
    fwlite::Handle<llvvGenParticleCollection> genPartCollHandle;
    genPartCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!genPartCollHandle.isValid())
    {
      std::cout << "llvvGenParticleCollection Object NotFound" << std::endl;
      continue;
    }
    gen = *genPartCollHandle;
    
    // Rho
    fwlite::Handle<double> rhoHandle;
    rhoHandle.getByLabel(ev, "kt6PFJets", "rho");
    if(!rhoHandle.isValid())
    {
      std::cout << "rho Object NotFound" << std::endl;
      continue;
    }
    eventContent.GetDouble("rho") = *rhoHandle;

    // Rho25
    fwlite::Handle<double> rho25Handle;
    rho25Handle.getByLabel(ev, "kt6PFJetsCentral", "rho");
    if(!rho25Handle.isValid())
    {
      std::cout << "rho25 Object NotFound" << std::endl;
      continue;
    }
    eventContent.GetDouble("rho25") = *rho25Handle;

    // Collection of leptons
    fwlite::Handle<llvvLeptonCollection> leptonCollHandle;
    leptonCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!leptonCollHandle.isValid())
    {
      std::cout << "llvvLeptonCollection Object NotFound" << std::endl;
      continue;
    }
    leptons = *leptonCollHandle;
    
    // Electron Information Collection
    fwlite::Handle<llvvElectronInfoCollection> electronInfoCollHandle;
    electronInfoCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!electronInfoCollHandle.isValid())
    {
      std::cout << "llvvElectronInfoCollection Object NotFound" << std::endl;
      continue;
    }

    // Muon Information Collection
    fwlite::Handle<llvvMuonInfoCollection> muonInfoCollHandle;
    muonInfoCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!muonInfoCollHandle.isValid())
    {
      std::cout << "llvvMuonInfoCollection Object NotFound" << std::endl;
      continue;
    }
    
    // Tau Collection
    fwlite::Handle<llvvTauCollection> tauCollHandle;
    tauCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!tauCollHandle.isValid())
    {
      std::cout << "llvvTauCollection Object NotFound" << std::endl;
      continue;
    }
    taus = *tauCollHandle;
    
    // Boosted tau Collection
    fwlite::Handle<llvvTauCollection> boostedTauCollHandle;
    boostedTauCollHandle.getByLabel(ev, "llvvObjectProducersUsed", "boosted");
    if(!boostedTauCollHandle.isValid())
    {
      std::cout << "llvvTauCollection Boosted Object NotFound" << std::endl;
    }
    boostedTaus = *boostedTauCollHandle;
    if(mergeBoostedTaus)
      for(size_t i = 0; i < boostedTaus.size(); ++i)
        taus.push_back(boostedTaus[i]);

    // Jet Collection
    fwlite::Handle<llvvJetCollection> jetCollHandle;
    jetCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!jetCollHandle.isValid())
    {
      std::cout << "llvvJetCollection Object NotFound" << std::endl;
      continue;
    }
    jets_ = *jetCollHandle;
    jets.clear();
    for(auto i = jetCollHandle->begin(); i != jetCollHandle->end(); ++i)
      jets.push_back(llvvJetExt(*i));
    for(auto& jet: jets)
    {
      if(debugEvent)
        analyserCout << "  Starting computation of jet energy scale and resolution uncertainties" << std::endl;
      // Apply jet corrections
      double toRawSF = jet.torawsf;
      LorentzVector rawJet(jet*toRawSF);
      jesCor->setJetEta(rawJet.eta());
      jesCor->setJetPt(rawJet.pt());
      jesCor->setJetA(jet.area);
      jesCor->setRho(static_cast<double>(eventContent.GetDouble("rho")));
//      jesCor->setNPV(nvtx); ?

      if(debugEvent)
        analyserCout << "   Uncorrected jet pt: " << jet.pt() << "; jet eta: " << jet.eta() << std::endl;

      double newJECSF(jesCor->getCorrection());
      rawJet *= newJECSF;
      jet.SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
      jet.torawsf = 1./newJECSF;

      // Compute scale and resolution uncertainties
      if(isMC)
      {
        if(debugEvent)
          analyserCout << "  Smearing JER" << std::endl;
        std::vector<float> smearPt = utils::cmssw::smearJER(jet.pt(),jet.eta(),jet.genj.pt());
        jet.jer     = smearPt[0];
        jet.jerup   = smearPt[1];
        jet.jerdown = smearPt[2];

        if(debugEvent)
          analyserCout << "   Scaled jet (JES) pt: " << jet.pt() << "; jet eta: " << jet.eta() << std::endl;

        double newJERSF = jet.jer/jet.pt();
        rawJet *= newJERSF;
        jet.SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
        jet.torawsf *= 1./newJERSF;

        if(debugEvent)
          analyserCout << "  Smearing JES" << std::endl;
        if(debugEvent)
          analyserCout << "   Smeared jet (JES & JER) pt: " << jet.pt() << "; jet eta: " << jet.eta() << std::endl;
        smearPt = utils::cmssw::smearJES(jet.pt(),jet.eta(), totalJESUnc);
        jet.jesup   = smearPt[0];
        jet.jesdown = smearPt[1];
        if(debugEvent)
          analyserCout << "  Done Smearing" << std::endl;
      }
      else
      {
        jet.jer     = jet.pt();
        jet.jerup   = jet.pt();
        jet.jerdown = jet.pt();
        jet.jesup   = jet.pt();
        jet.jesdown = jet.pt();
      }
    }

    // MET Collection
    fwlite::Handle<llvvMet> metHandle;
    metHandle.getByLabel(ev, "llvvObjectProducersUsed", "pfMETPFlow");
//    metHandle.getByLabel(ev, "llvvObjectProducersUsed", "pfType1CorrectedMet");
    if(!metHandle.isValid())
    {
      std::cout << "llvvMet Object NotFound" << std::endl;
      continue;
    }
    //eventContent.GetDouble("MET") = *metHandle;
    metVec = *metHandle;

    // Trigger Prescales
    fwlite::Handle<std::vector<int> > triggerPrescalesHandle;
    triggerPrescalesHandle.getByLabel(ev, "llvvObjectProducersUsed", "triggerPrescales");
    if(!triggerPrescalesHandle.isValid())
    {
      std::cout << "triggerPrescales Object NotFound" << std::endl;
      continue;
    }
    triggerPrescales = *triggerPrescalesHandle;

    if(debugEvent)
      analyserCout << " Finished loading collections" << std::endl;

    eventContent.GetInt("RunNo")   = ev.eventAuxiliary().run();
    eventContent.GetInt("LumiNo")  = ev.eventAuxiliary().luminosityBlock();
    eventContent.GetInt("EventNo") = ev.eventAuxiliary().event();

    ProcessEvent(iev);
    eventContent.Lock();
    FillHistograms();

    if(!doneFirstEvent)
    {
      if(outputEventList)
        eventContent.OutputEventListHeader(eventListFile, priorityOutput);
      if(saveSummaryTree)
        eventContent.SetSummaryTreeBranches(summaryTree);
    }
    if(outputEventList)
      eventContent.OutputEventList(eventListFile, priorityOutput);
    if(saveSummaryTree && (keepAllEvents || static_cast<bool>(eventContent.GetBool("selected"))))
    {
      TDirectory* cwd = gDirectory;
      summaryOutTFile->cd();
      summaryTree->Fill();
      cwd->cd();
    }
    doneFirstEvent = true;

    if(limitEvents != 0)
    {
      if(nEventsOut >= limitEvents - 1)
        break;
      else
        ++nEventsOut;
    }
  }

  // Output temporary buffer and restore cout and cerr behaviour
  std::cout.rdbuf(coutbuf);
  std::cerr.rdbuf(cerrbuf);
  std::cout << std::endl;
  if(saveRedirect)
    std::cout << analyserBuffer.str();

  std::cout << "totalEntries: " << totalEntries << "; vs nInitEvent: " << nInitEvent << ";" << std::endl;

  std::cout << "Saving histograms in " << outFile << std::endl;
  TFile* outfile = new TFile(outFile.c_str(), "RECREATE");
  histMonitor.Write();
  outfile->Close();
  delete outfile;

  if(saveSummaryTree)
  {
    TDirectory* cwd = gDirectory;
    summaryOutTFile->cd();
    summaryTree->Write();
    summaryOutTFile->Close();
    delete summaryOutTFile;
    cwd->cd();
  }

  if(outputEventList)
  {
    eventListFile.close();
  }
  
  isSetup = false;

  return;
}

void Analyser::InitHistograms()
{
  histMonitor.addHistogram(new TH1D("nup", ";NUP;Events", 10, 0, 10));
  histMonitor.addHistogram(new TH1D("nvtx", ";Vertices;Events", 50, -0.5, 49.5));
  histMonitor.addHistogram(new TH1D("nvtxraw", ";Vertices;Events", 50, -0.5, 49.5));
  histMonitor.addHistogram(new TH1D("rho", ";#rho;Events", 25, 0, 25));
  histMonitor.addHistogram(new TH1D("rho25", ";#rho(#eta<2.5);Events", 25, 0, 25));
  
  // MET
  histMonitor.addHistogram(new TH1D("MET", ";MET [GeV];Events", 25, 0, 200));

  if(debug)
    std::cout << "Finished Analyser::InitHistograms()" << std::endl;
  
  UserInitHistograms();
  return;
}

void Analyser::FillHistograms()
{
  auto& weight = eventContent.GetDouble("weight").Value();
  auto& puWeight = eventContent.GetDouble("PUweight").Value();
  auto& selected = eventContent.GetBool("selected").Value();

  histMonitor.fillHisto("nup", "", genEv.nup, 1);
  if(selected)
  {
    histMonitor.fillHisto("nvtx", chTags, eventContent.GetInt("nvtx").Value(), weight);
    histMonitor.fillHisto("nvtxraw", chTags, eventContent.GetInt("nvtx").Value(), weight/puWeight);

    histMonitor.fillHisto("rho", chTags, eventContent.GetDouble("rho").Value(), weight);
    histMonitor.fillHisto("rho25", chTags, eventContent.GetDouble("rho25").Value(), weight);
  }
  
  UserFillHistograms();
}

void Analyser::EventContentSetup()
{
  if(debug)
    std::cout << "Running Analyser::EventContentSetup()" << std::endl;
  
  auto& runNumber = eventContent.AddInt("RunNo", 0);
  runNumber.AddMetadata("eventlist", "true");
  runNumber.AddMetadata("eventlistWidth", "10");
  
  auto& lumiNumber = eventContent.AddInt("LumiNo", 0);
  lumiNumber.AddMetadata("eventlist", "true");
  lumiNumber.AddMetadata("eventlistWidth", "10");
  
  auto& eventNumber = eventContent.AddInt("EventNo", 0);
  eventNumber.AddMetadata("eventlist", "true");
  eventNumber.AddMetadata("eventlistWidth", "10");
  
  auto& nvtx = eventContent.AddInt("nvtx", -1);
  
  auto& rho   = eventContent.AddDouble("rho", -1);
  auto& rho25 = eventContent.AddDouble("rho25", -1);
  
  auto& weight = eventContent.AddDouble("weight", 1);
  auto& PUweight = eventContent.AddDouble("PUweight", 1);
  auto& xsecweight = eventContent.AddDouble("xsecweight", 1);
  auto& xsec = eventContent.AddDouble("crossSection", 1);
  if(runSystematics && isMC)
  {
    xsec("xsec_UP");
    xsec("xsec_DOWN");
    xsec("PDF_UP");
    xsec("PDF_DOWN");
    xsec.Lock();
  }
  weight.AddMetadata("eventlist", "true");
  PUweight.AddMetadata("eventlist", "true");
  xsecweight.AddMetadata("eventlist", "true");

  auto& selected = eventContent.AddBool("selected", false);
  selected.AddMetadata("eventlist", "true");
  selected.AddMetadata("eventlistWidth", "8");
  
  auto& met = eventContent.AddDouble("MET", -20.0);
  met.AddMetadata("eventtree", "true"); // If this metadata is not defined, it is assumed to be true, only set it to false for variables not to be in the eventtree
  met.AddMetadata("eventlist", "true"); // If this metadata is not defined, it is assumed to be false. If true, the base variable will be output in the event list
  met.AddMetadata("eventlistWidth", "12"); // This metadata will only be considered if eventlist metadata is true. In that situation this field is used to define the width, in characters of this variable in the eventlist
  if(runSystematics && isMC)
  {
    met.Systematic("JES_UP"); // As an alternative you can also write:   met("JES_UP");
    met.Systematic("JES_DOWN");
    met.Systematic("JER_UP"); // As an alternative you can also write:   met("JES_UP");
    met.Systematic("JER_DOWN");
  }

  if(debug)
    std::cout << "Finished Analyser::EventContentSetup()" << std::endl;

  UserEventContentSetup();

//  eventContent.Lock(); // It should only be locked after the first iteration through the event loop
  return;
}

void Analyser::ProcessEvent(size_t iev)
{
  chTags.clear();
  chTags.push_back("all");

  if(isMC)
  {
    auto& PUweight = eventContent.GetDouble("PUweight");
    PUweight = LumiWeights->weight(genEv.ngenITpu) * PUNorm[0];
    if(runSystematics && isMC)
    {
      PUweight("PU_UP")   = PUweight.Value() * PuShifters[utils::cmssw::PUUP  ]->Eval(genEv.ngenITpu) * (PUNorm[2]/PUNorm[0]);
      PUweight("PU_DOWN") = PUweight.Value() * PuShifters[utils::cmssw::PUDOWN]->Eval(genEv.ngenITpu) * (PUNorm[1]/PUNorm[0]);
    }

    eventContent.GetDouble("crossSection") = crossSection_;
    eventContent.GetDouble("xsecweight") = xsecWeight;
  }

  MET = getMETvariations();
  eventContent.GetDouble("MET") = MET.Pt();

  UserProcessEvent(iev);

  ValueWithSystematics<double> pdfvar(1.0);
  
  if(isMC && runSystematics)
  {
    pdfvar("PDFVAR_UP");
    pdfvar("PDFVAR_DOWN");
    
    for(auto &entry: pdfVariations)
    {
      std::vector<double> weights = entry->getWeights(iev);
      for(auto& iWeight: weights)
      {
        if(iWeight > pdfvar("PDFVAR_UP"))
          pdfvar("PDFVAR_UP") = iWeight;
        if(iWeight < pdfvar("PDFVAR_DOWN"))
          pdfvar("PDFVAR_DOWN") = iWeight;
      }
    }
  }

  if(isMC)
    eventContent.GetDouble("weight") *= eventContent.GetDouble("PUweight") * eventContent.GetDouble("xsecweight");

  eventContent.GetDouble("weight") *= pdfvar;
  
  return;
}

template<class T>
void Analyser::loadSystematics(std::vector<std::string>& list, ValueWithSystematics<T> variable)
{
  for(auto& syst: variable.Systematics())
  {
    if(std::find(list.begin(), list.end(), syst) == list.end())
    {
      list.push_back(syst);
    }
  }
  
  return;
}

ValueWithSystematics<TLorentzVector> Analyser::getMETvariations()
{
  ValueWithSystematics<llvvMet> tmpVal(metVec);
  if(!isMC)
    return tmpVal.ToTLorentzVector();
  ValueWithSystematics<TLorentzVector> retVal(tmpVal.ToTLorentzVector());

  TLorentzVector nullVec(0, 0, 0, 0);

  ValueWithSystematics<TLorentzVector> clusteredFlux(nullVec);
  for(auto& jet: jets)
  {
    TLorentzVector tmp(jet.Px(), jet.Py(), jet.Pz(), jet.E());
    ValueWithSystematics<TLorentzVector> jetVec(tmp);
    
    ValueWithSystematics<double> scale(jet.Pt());
    if(runSystematics)
    {
      scale("JES_UP") = jet.jesup;
      scale("JES_DOWN") = jet.jesdown;
      scale("JER_UP") = jet.jerup;
      scale("JER_DOWN") = jet.jerdown;
    }
    scale /= jet.Pt();
    
    jetVec *= scale;
    clusteredFlux += jetVec;
  }

  ValueWithSystematics<TLorentzVector> leptonFlux(nullVec);
  for(auto& lep: leptons)
  {
    TLorentzVector tmp(lep.Px(), lep.Py(), lep.Pz(), lep.E());
    ValueWithSystematics<TLorentzVector> lepVec(tmp);
    
    ValueWithSystematics<double> scale(1);
    if(runSystematics)
    {
      if(abs(lep.id) == 13)
      {
        scale("LES_UP") = 1.01;
        scale("LES_DOWN") = 0.99;
      }
      else
      {
        if(std::abs(lep.electronInfoRef->sceta) < 1.442)
        {
          scale("LES_UP") = 1.02;
          scale("LES_DOWN") = 0.98;
        }
        else
        {
          scale("LES_UP") = 1.05;
          scale("LES_DOWN") = 0.95;
        }
      }
    }
    
    lepVec *= scale;
    leptonFlux += lepVec;
  }
  
  ValueWithSystematics<TLorentzVector> unclusteredFlux = -((retVal + clusteredFlux.Value() + leptonFlux.Value()).Value());
  if(runSystematics)
  {
    unclusteredFlux("UMET_UP") *= 1.1;
    unclusteredFlux("UMET_DOWN") *= 0.9;
  }

  retVal = (-unclusteredFlux) - clusteredFlux - leptonFlux;
  if(runSystematics)
  {
    retVal("LES_UP");
    retVal("LES_DOWN");
    retVal("UMET_UP");
    retVal("UMET_DOWN");
    retVal("JES_UP");
    retVal("JES_DOWN");
    retVal("JER_UP");
    retVal("JER_DOWN");
  }
  
  return retVal;
}

class StauAnalyser : public Analyser
{
public:
  StauAnalyser(std::string cfgFile);
  
  enum class IDType {LooseID, MediumID, TightID};
  enum class TAU_E_ID {antiELoose, antiEMedium, antiETight, antiEMva, antiEMva3Loose, antiEMva3Medium, antiEMva3Tight, antiEMva3VTight, antiEMva5Medium};

private:

protected:
  bool exclusiveRun;
  double stauMtoPlot;
  double neutralinoMtoPlot;
  bool doSVfit;
  bool keepOnlyPromptTaus;
  
  bool isStauStau;
  
  double sqrtS;
  double minElPt;
  double maxElEta;
  double ECALGap_MinEta;
  double ECALGap_MaxEta;
  double minMuPt;
  double maxMuEta;
  double minTauPt;
  double maxTauEta;
  double minJetPt;
  double maxJetEta;
  double maxElDz;
  double maxElD0;
  double maxElDzVeto;
  double maxElD0Veto;
  double maxMuDz;
  double maxMuD0;
  double maxMuDzVeto;
  double maxMuD0Veto;
  double elIso;
  double elIsoVeto;
  double muIso;
  double muIsoVeto;
  double minElPtVeto;
  double maxElEtaVeto;
  double minMuPtVeto;
  double maxMuEtaVeto;
  double maxTauDz;
  double genMatchRCone;
  
  ValueWithSystematics<TH1*> fakeRateHist;
  ValueWithSystematics<TH1*> promptRateHist;
  ValueWithSystematics<TH1*> xSecHist;

  virtual void UserLoadCfgOptions();
  virtual void UserSetup();
  virtual void UserProcessEvent(size_t iev);
  virtual void UserInitHistograms();
  virtual void UserEventContentSetup();
  virtual void UserFillHistograms();
  
  ValueWithSystematics<double> LeptonTauTriggerScaleFactor(llvvLepton& lepton, llvvTau& tau);
  ValueWithSystematics<double> StauCrossSec();
  double Efficiency(double m, double m0, double sigma, double alpha, double n, double norm);
  bool electronMVAID(double mva, llvvLepton& lepton, IDType id);
  ValueWithSystematics<double> leptonIdAndIsoScaleFactor(ValueWithSystematics<llvvLepton>& lepton);
  ValueWithSystematics<double> tauScaleFactor(ValueWithSystematics<llvvTau>& tau, TAU_E_ID eId);
  
  ValueWithSystematics<double> computeMT2(const ValueWithSystematics<llvvTau>& tau, const ValueWithSystematics<llvvLepton>& lep, const ValueWithSystematics<TLorentzVector>& met);

};

StauAnalyser::StauAnalyser(std::string cfgFile): Analyser(cfgFile)
{
}

void StauAnalyser::UserLoadCfgOptions()
{
  exclusiveRun = cfgOptions.getParameter<bool>("exclusiveRun");

  stauMtoPlot        =   120;
  neutralinoMtoPlot  =    20; // Default mass point to place in plots
  doSVfit            = false;
  keepOnlyPromptTaus = false;

  if(cfgOptions.exists("stauMtoPlot"))
    stauMtoPlot        = cfgOptions.getParameter<double>("stauMtoPlot");
  if(cfgOptions.exists("neutralinoMtoPlot"))
    stauMtoPlot        = cfgOptions.getParameter<double>("neutralinoMtoPlot");
  if(cfgOptions.exists("doSVfit"))
    doSVfit            = cfgOptions.getParameter<bool>("doSVfit");
  if(cfgOptions.exists("keepOnlyPromptTaus"))
    keepOnlyPromptTaus = cfgOptions.getParameter<bool>("keepOnlyPromptTaus");

  // Consider setting here the cut values etc, will have to be added to the cfg file
  sqrtS          =  8;      // Center of mass energy
  minElPt        = 24;      // Selected electron pT and eta
  maxElEta       =  2.1;
  ECALGap_MinEta =  1.4442; // ECAL gap parameters
  ECALGap_MaxEta =  1.5660;
  minMuPt        = 20;      // Selected muon pT and eta
  maxMuEta       =  2.1;
  minTauPt       = 20;      // Selected tau pT and eta (I was using 25)
  maxTauEta      =  2.3;
  minJetPt       = 30;
  maxJetEta      =  4.7;    // Selected jet eta

  maxElDz        =  0.1;
  maxElD0        =  0.045;
  maxElDzVeto    =  0.2;
  maxElD0Veto    =  0.045;

  maxMuDz        =  0.5;
  maxMuD0        =  0.2;
  maxMuDzVeto    =  0.2;  // TODO: Does this make sense? it should probably at least be equal to maxMuDz
  maxMuD0Veto    =  0.2;

  maxTauDz       =  0.5;
  
  elIso          =  0.1;
  elIsoVeto      =  0.3;
  muIso          =  0.1;
  muIsoVeto      =  0.3;
  
  minElPtVeto    = 10;
  maxElEtaVeto   =  2.3;
  minMuPtVeto    = 10;
  maxMuEtaVeto   =  2.4;
  
  genMatchRCone  =  0.3;

  if(debug)
    std::cout << "Finished StauAnalyser::LoadCfgOptions()" << std::endl;

  return;
}

void StauAnalyser::UserSetup()
{
  TDirectory* cwd = gDirectory;

  if(doDDBkg)
  {
    for(auto & file : fileList)
    {
      if(file.find("DD", 0) != std::string::npos)
        file.replace(file.find("DD", 0), 2, "");
    }
    for(auto & file : fileList)
    {
      std::cout << "New input file name: " << file << std::endl;
    }


    std::string RatesFileName = gSystem->ExpandPathName("$CMSSW_BASE/src/UserCode/llvv_fwk/data/TStauStau/rates.root");
    std::cout << "Trying to open rates file: " << RatesFileName << std::endl;
    TFile RatesFile(RatesFileName.c_str(), "READ");
    if(!RatesFile.IsOpen())
      throw AnalyserException("Unable to open rates file.");
    cwd->cd();

    fakeRateHist   = static_cast<TH1*>(RatesFile.Get("data-Zprompt/data-Zprompt_InvMET_OS_etaSelectedTau_FR")->Clone("fakeRate"));
//    promptRateHist = static_cast<TH1*>(RatesFile.Get("Z #rightarrow ll/Zrightarrowll_InvMET_OS_Prompt_etaSelectedTau_FR")->Clone("promptRate"));
    promptRateHist = static_cast<TH1*>(RatesFile.Get("Z #rightarrow ll/Zrightarrowll_InvMET_OS_etaSelectedTau_FR")->Clone("promptRate"));

    if(fakeRateHist.Value() == NULL)
    {
      throw AnalyserException("Unable to open fake rate histogram.");
    }
    if(promptRateHist.Value() == NULL)
    {
      throw AnalyserException("Unable to open prompt rate histogram.");
    }
  }

  TString turl(fileList[0]);
  isStauStau = isMC && turl.Contains("TStauStau");
  if(isStauStau)
  {
    std::string xSecFileName  = gSystem->ExpandPathName("$CMSSW_BASE/src/UserCode/llvv_fwk/data/TStauStau/StauCrossSections.root");
    std::cout << "Trying to open stau cross sections file: " << xSecFileName << std::endl;
    TFile xSecFile(xSecFileName.c_str(), "READ");
    if(!xSecFile.IsOpen())
      throw AnalyserException("Unable to open stau cross sections file.");
    cwd->cd();
    
    xSecHist = NULL;
    xSecHist("xsec_UP");
    xSecHist("xsec_DOWN");

    xSecHist.Value() = static_cast<TH1*>(xSecFile.Get("xsec")->Clone("xsec"));
    xSecHist("xsec_UP") = static_cast<TH1*>(xSecFile.Get("xsec_UP")->Clone("xsec_UP"));
    xSecHist("xsec_DOWN") = static_cast<TH1*>(xSecFile.Get("xsec_DOWN")->Clone("xsec_DOWN"));

    if(xSecHist.Value() == NULL || xSecHist("xsec_UP") == NULL || xSecHist("xsec_DOWN") == NULL)
    {
      throw AnalyserException("Unable to open the stau cross section histograms.");
    }
  }

  if(debug)
    std::cout << "Finished StauAnalyser::UserSetup()" << std::endl;

  return;
}

void StauAnalyser::UserProcessEvent(size_t iev)
{
  std::vector<std::string> tmpLoop;
  auto& dropEvent = eventContent.GetBool("dropEvent");
  dropEvent = false;
  /**** Remove double counting if running on exclusive samples ****/
  if(exclusiveRun && isV0JetsMC)
  {
    if(genEv.nup > 5)
    {
      dropEvent = true;
    }
  }

  /**** Ensure that for the TStauStau dataset the LHE event info with the generation comments is there, it is needed to know the generated masses ****/
  if(isStauStau)
  {
    if(!LHEHandle.isValid())
    {
      std::cout << "LHEEventProduct Object not Found for TStauStau dataset" << std::endl;
      dropEvent = true;
    }
    else
    {
      if(LHEHandle->comments_size() == 0)
      {
        std::cout << "LHEEventProduct Object Found but empty for TStauStau dataset" << std::endl;
        dropEvent = true;
      }
      else
      {
        for(auto comment = LHEHandle->comments_begin(); comment != LHEHandle->comments_end(); ++comment)
        {
          auto modelPos = comment->find("# model TStauStau_");
          if(modelPos != std::string::npos)
          {
            double stauMass = 0, neutralinoMass = 0;
            std::stringstream tmp;
            auto numPos = comment->find_first_of("1234567890", modelPos);

            tmp << comment->substr(numPos, comment->find("_", numPos)-numPos);
            tmp >> stauMass;
            tmp.clear();

            numPos = comment->find("_", numPos);
            numPos = comment->find_first_of("1234567890", numPos);
            tmp << comment->substr(numPos, comment->find("\n", numPos)-numPos);
            tmp >> neutralinoMass;

            eventContent.GetDouble("stauMass") = stauMass;
            eventContent.GetDouble("neutralinoMass") = neutralinoMass;

            break;
          }
        }
      }
    }
  }

  // Moving this to the end to avoid uninitialised systematics
/*  if(dropEvent.Value())
  {
    eventContent.GetBool("selected") = false;
    return;
  }// */

  if(isStauStau)
  {
    int nEvents = 10000;
    crossSection_ = StauCrossSec();
    eventContent.GetDouble("crossSection") = crossSection_;
    eventContent.GetDouble("xsecweight")  = crossSection_/nEvents;
  }

  //bool singleETrigger  = triggerBits[13]; // HLT_Ele27_WP80_v*
  //bool singleMuTrigger = triggerBits[15]; // HLT_IsoMu24_v*
  bool TauPlusE2012A  = triggerBits[18]; // HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*
  bool TauPlusMu2012A = triggerBits[22]; // HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v*
  bool TauPlusE2012B  = triggerBits[17]; // HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v*
  bool TauPlusMu2012B = triggerBits[21]; // HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*

  bool TauPlusETrigger = TauPlusE2012A || TauPlusE2012B;
  bool TauPlusMuTrigger = TauPlusMu2012A || TauPlusMu2012B;
  auto& triggeredOn = eventContent.GetBool("triggeredOn");
  triggeredOn = TauPlusETrigger || TauPlusMuTrigger;

  // Get trigger Scale Factor
  if(applyScaleFactors && isMC)
  {
    auto& triggerSF = eventContent.GetDouble("triggerSF");
    triggerSF = 1;
    if(debugEvent)
    {
      if(TauPlusETrigger)
        analyserCout << "  Triggered on by TauPlusE\n";
      if(TauPlusMuTrigger)
        analyserCout << "  Triggered on by TauPlusMu\n";

      if(triggeredOn.Value())
      {
        analyserCout << "  Looping on leptons:\n";
        for(auto& lep: leptons)
          analyserCout << "    Lepton (" << lep.id << ", pT=" << lep.pt() << ") trigger bits: " << bitset<8*sizeof(int)>(lep.Tbits) << "\n";

        for(auto& tau: taus)
          analyserCout << "    Tau (pT=" << tau.pt() << ") trigger bits: " << bitset<8*sizeof(int)>(tau.Tbits) << "\n";
      }
      analyserCout << std::flush;
    }

    if(TauPlusETrigger)
    {
      llvvTau* trigTau = NULL, *leadTau = NULL;
      llvvLepton* trigE = NULL, *leadE = NULL;

      //Sometimes the triggered lepton can not be found, so we use the leading lepton instead
      for(auto& lep: leptons)
      {
        if(lep.Tbits & (3 << 17))
        {
          if(trigE == NULL) trigE = &lep;
          else if(lep.pt() > trigE->pt()) trigE = &lep;
        }

        if(abs(lep.id) == 11)
        {
          if(leadE == NULL) leadE = &lep;
          else if(lep.pt() > leadE->pt()) leadE = &lep;
        }
      }
      if(trigE == NULL)
        trigE = leadE;

      //For the taus (at the moment) Tbits is filled randomnly, so we only use the leading tau
      for(auto& tau: taus)
      {
//        if(tau.Tbits & (3 << 17))
//        {
//          if(trigTau == NULL) trigTau = &tau;
//          else if(tau.pt() > trigTau->pt()) trigTau = &tau;
//        }

        if(leadTau == NULL) leadTau = &tau;
        else if(tau.pt() > leadTau->pt()) leadTau = &tau;
      }
      if(trigTau == NULL)
        trigTau = leadTau;

      if(trigTau != NULL && trigE != NULL)
      {
        triggerSF *= (LeptonTauTriggerScaleFactor(*trigE, *trigTau));
      }
      else
      {
        if(debugEvent)
        {
          if(trigE == NULL)
            analyserCout << " TauPlusE trigSF: Unable to find triggered electron" << std::endl;
          if(trigTau == NULL)
            analyserCout << " TauPlusE trigSF: Unable to find triggered tau" << std::endl;
        }
      }
    }

    if(TauPlusMuTrigger)
    {
      llvvTau* trigTau = NULL, *leadTau = NULL;
      llvvLepton* trigMu = NULL, *leadMu = NULL;

      //Sometimes the triggered lepton can not be found, so we use the leading lepton instead
      for(auto& lep: leptons)
      {
        if(lep.Tbits & (3 << 21))
        {
          if(trigMu == NULL) trigMu = &lep;
          else if(lep.pt() > trigMu->pt()) trigMu = &lep;
        }

        if(abs(lep.id) == 13)
        {
          if(leadMu == NULL) leadMu = &lep;
          else if(lep.pt() > leadMu->pt()) leadMu = &lep;
        }
      }
      if(trigMu == NULL)
        trigMu = leadMu;

      //For the taus (at the moment) Tbits is filled randomnly, so we only use the leading tau
      for(auto& tau: taus)
      {
//        if(tau.Tbits & (3 << 21))
//        {
//          if(trigTau == NULL) trigTau = &tau;
//          else if(tau.pt() > trigTau->pt()) trigTau = &tau;
//        }

        if(leadTau == NULL) leadTau = &tau;
        else if(tau.pt() > leadTau->pt()) leadTau = &tau;
      }
      if(trigTau == NULL)
        trigTau = leadTau;

      if(trigTau != NULL && trigMu != NULL)
      {
        triggerSF *= LeptonTauTriggerScaleFactor(*trigMu, *trigTau);
      }
      else
      {
        if(debugEvent)
        {
          if(trigMu == NULL)
            analyserCout << " TauPlusMu trigSF: Unable to find triggered muon" << std::endl;
          if(trigTau == NULL)
            analyserCout << " TauPlusMu trigSF: Unable to find triggered tau" << std::endl;
        }
      }
    }

    if(debugEvent)
      analyserCout << "  Computed trigger SF: " << triggerSF.Value() << std::endl;

    eventContent.GetDouble("weight") *= triggerSF;
  }

  // Get Leptons
  if(debugEvent)
  {
    analyserCout << " Finished computing PU weight and trigger scale factors" << std::endl;
    analyserCout << " Getting leptons" << std::endl;
  }
  ValueWithSystematics<std::vector<llvvLepton>> selLeptons;
  if(runSystematics && isMC)
  {
    selLeptons("LES_UP");
    selLeptons("LES_DOWN");
    selLeptons.Lock();
  }
  for(auto& lep: leptons)
  {
    int lepId = abs(lep.id);

    double sf = 0.01;
    if(lepId == 11)
    {
      if(std::abs(lep.electronInfoRef->sceta) < 1.442)
        sf = 0.02;
      else
        sf = 0.05;
    }

    if(lepId == 13 && muCor)
    {
      TLorentzVector p4(lep.px(), lep.py(), lep.pz(), lep.energy());
//Fixed?      muCor->applyPtCorrection(p4, (lep.id>0)?1:-1);
      muCor->applyPtCorrection(p4, (lepId>0)?1:-1);
      if(isMC)
//Fixed?        muCor->applyPtSmearing(p4, (lep.id>0)?1:-1, false);
        muCor->applyPtSmearing(p4, (lepId>0)?1:-1, false);
      lep.SetPxPyPzE(p4.Px(), p4.Py(), p4.Pz(), p4.Energy());
    }

    // Lepton Kinematics
    double eta = (lepId == 11)?(lep.electronInfoRef->sceta):(lep.eta());
    ValueWithSystematics<bool> keepKin(true), passKin(true);
    if(lepId == 11) // If Electron
    {
      if(std::abs(eta) > maxElEta)
        passKin = false;

      if(std::abs(eta) > maxElEtaVeto)
        keepKin = false;

      if(std::abs(eta) > ECALGap_MinEta && std::abs(eta) < ECALGap_MaxEta) // Remove electrons that fall in ECAL Gap
      {
        passKin = false;
        keepKin = false;
      }
      
      if(runSystematics && isMC)
      {
        passKin("LES_UP");
        passKin("LES_DOWN");
        keepKin("LES_UP");
        keepKin("LES_DOWN");
      }
      if(lep.pt() < minElPt)
        passKin.Value() = false;
      if(lep.pt() < minElPtVeto)
        keepKin.Value() = false;
      if(runSystematics && isMC)
      {
        if(lep.pt()*(1+sf) < minElPt)
          passKin("LES_UP") = false;
        if(lep.pt()*(1+sf) < minElPtVeto)
          keepKin("LES_UP") = false;
        if(lep.pt()*(1-sf) < minElPt)
          passKin("LES_DOWN") = false;
        if(lep.pt()*(1-sf) < minElPtVeto)
          keepKin("LES_DOWN") = false;
      }
    }
    else // If Muon
    {
      if(std::abs(eta) > maxMuEta)
        passKin = false;

      if(std::abs(eta) > maxMuEtaVeto)
        keepKin = false;

      if(runSystematics && isMC)
      {
        passKin("LES_UP");
        passKin("LES_DOWN");
        keepKin("LES_UP");
        keepKin("LES_DOWN");
      }
      if(lep.pt() < minMuPt)
        passKin.Value() = false;
      if(lep.pt() < minMuPtVeto)
        keepKin.Value() = false;
      if(runSystematics && isMC)
      {
        if(lep.pt()*(1+sf) < minMuPt)
          passKin("LES_UP") = false;
        if(lep.pt()*(1+sf) < minMuPtVeto)
          keepKin("LES_UP") = false;
        if(lep.pt()*(1-sf) < minMuPt)
          passKin("LES_DOWN") = false;
        if(lep.pt()*(1-sf) < minMuPtVeto)
          keepKin("LES_DOWN") = false;
      }
    }

    // Lepton ID
    bool passID = true, keepID = true;
    Int_t idbits = lep.idbits;
    if(lepId == 11)
    {
      // bool isTight = electronMVAID(lep.electronInfoRef->mvanontrigv0, lep, IDType::TightID);
      // bool isLoose = electronMVAID(lep.electronInfoRef->mvanontrigv0, lep, IDType::LooseID);
      // bool isLoose = ((idbits >> 4) & 0x1);
      // bool isTight = ((idbits >> 6) & 0x1);
      passID = electronMVAID(lep.electronInfoRef->mvanontrigv0, lep, IDType::LooseID);
      keepID = passID;
      if(lep.d0 > maxElD0)
        passID = false;
      if(lep.dZ > maxElDz)
        passID = false;
      if(lep.d0 > maxElD0Veto)
        keepID = false;
      if(lep.dZ > maxElDzVeto)
        keepID = false;

      if(lep.electronInfoRef->isConv)
      {
        passID = false;
        keepID = false;
      }
      if(lep.trkLostInnerHits > 0)
      {
        passID = false;
        keepID = false;
      }
    }
    else
    {
      // bool isLoose = ((idbits >> 8) & 0x1);
      // bool isTight = ((idbits >> 10) & 0x1);
      passID = ((idbits >> 10) & 0x1);
      keepID = ((idbits >> 8) & 0x1);
      if(lep.d0 > maxMuD0)
        passID = false;
      if(lep.d0 > maxMuD0Veto)
        keepID = false;
      if(lep.dZ > maxMuDz)
        passID = false;
      if(lep.dZ > maxMuDzVeto)
        keepID = false;

      if(passID)
        keepID = true;
    }

    // Lepton Isolation
    bool passIso = true, keepIso = true;
    double relIso = utils::cmssw::relIso(lep, eventContent.GetDouble("rho").Value());
    if(lepId == 11)
    {
      if(relIso > elIso)
        passIso = false;
      if(relIso > elIsoVeto)
        keepIso = false;
    }
    else
    {
      if(relIso > muIso)
        passIso = false;
      if(relIso > muIsoVeto)
        keepIso = false;
    }

    // Keep desired leptons
    if(static_cast<bool>(keepKin) && keepID && keepIso)
    {
      if(runSystematics && isMC)
      {
        if(keepKin.GetSystematicOrValue("LES_UP"))
          selLeptons.Systematic("LES_UP").push_back(lep*(1+sf));
        if(keepKin.GetSystematicOrValue("LES_DOWN"))
          selLeptons.Systematic("LES_DOWN").push_back(lep*(1-sf));
      }
      if(keepKin.Value())
        selLeptons.Value().push_back(lep);
    }

    if(!(triggeredOn.Value()))
      continue;

    ValueWithSystematics<double> weight = (eventContent.GetDouble("weight") * eventContent.GetDouble("PUweight") * eventContent.GetDouble("xsecweight"));
    histMonitor.fillHisto("leptonCutFlow", chTags, 0, weight.Value());
    if(passID)
    {
      histMonitor.fillHisto("leptonCutFlow", chTags, 1, weight.Value());
      if(passKin.Value())
      {
        histMonitor.fillHisto("leptonCutFlow", chTags, 2, weight.Value());
        if(passIso)
          histMonitor.fillHisto("leptonCutFlow", chTags, 3, weight.Value());
      }
    }
  }

  // Get Taus
  if(debugEvent)
    analyserCout << " Getting taus" << std::endl;
  ValueWithSystematics<std::vector<llvvTau>> selTaus;
  if(runSystematics && isMC)
  {
    selTaus("TES_UP");
    selTaus("TES_DOWN");
    for(auto& syst: selLeptons.Systematics())
      selTaus(syst);
    selTaus.Lock();
  }
  for(auto& tau: taus)
  {
    // Tau Kinematics
    ValueWithSystematics<bool> passKin = true;
    if(std::abs(tau.eta()) > maxTauEta)
      passKin = false;
    if(runSystematics && isMC)
    {
      passKin("TES_UP");
      passKin("TES_DOWN");
    }
    if(tau.pt() < minTauPt)
      passKin.Value() = false;
    if(runSystematics && isMC)
    {
      if(tau.pt()*1.03 < minTauPt)
        passKin("TES_UP") = false;
      if(tau.pt()*0.97 < minTauPt)
        passKin("TES_DOWN") = false;
    }

    // Tau overlap with selected leptons
    ValueWithSystematics<bool> passIso(true);
    tmpLoop.clear();
    tmpLoop.push_back("Value");
    if(runSystematics && isMC)
    {
      loadSystematics(tmpLoop, selLeptons);
      for(auto& syst: selLeptons.Systematics())
      {
        passIso(syst);
      }
    }
    for(auto& val: tmpLoop)
    {
      for(auto& lep: selLeptons.GetSystematicOrValue(val))
      {
        int lepId = abs(lep.id);
        if(lepId == 11)
        {
          if(lep.pt() < minElPt)
            continue;
          if(std::abs(lep.dZ) > maxElDz)
            continue;
          double eta = lep.electronInfoRef->sceta;
          if(std::abs(eta) > maxElEta)
            continue;
          double relIso = utils::cmssw::relIso(lep, eventContent.GetDouble("rho").Value());
          if(relIso > elIso)
            continue;
        }
        else
        {
          if(lep.pt() < minMuPt)
            continue;
          if(std::abs(lep.eta()) > maxMuEta)
            continue;
          double relIso = utils::cmssw::relIso(lep, eventContent.GetDouble("rho").Value());
          if(relIso > muIso)
            continue;
          Int_t idbits = lep.idbits;
          bool isTight = ((idbits >> 10) & 0x1);
          if(!isTight)
            continue;
        }

        if(deltaR(tau, lep) < 0.5)
        {
          passIso.GetSystematicOrValue(val) = false;
          break;
        }
      }
    }
    

    bool passQual = true;
    if(std::abs(tau.dZ) > maxTauDz)
      passQual = false;

    // Tau ID
    bool passID = true;
    if(!tau.passId(llvvTAUID::decayModeFinding)) passID = false;
    if(doDDBkg)
    {
      if(!tau.passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits)) passID = false;
    }
    else
    {
      if(!tau.passId(llvvTAUID::byTightCombinedIsolationDeltaBetaCorr3Hits)) passID = false;
    }
    if(!tau.passId(llvvTAUID::againstMuonTight3)) passID = false;
    if(!tau.passId(llvvTAUID::againstElectronMediumMVA5)) passID = false;

    if(passID && static_cast<bool>(passKin) && tau.isPF && static_cast<bool>(passIso) && passQual)
    {
      tmpLoop.clear();
      tmpLoop.push_back("Value");
      if(runSystematics && isMC)
      {
        loadSystematics(tmpLoop, passIso);
      }

      for(auto& val: tmpLoop)
      {
        if(val == "Value")
        {
          if(runSystematics && isMC)
          {
            if(passKin.GetSystematicOrValue("TES_UP") && passIso.GetSystematicOrValue(val))
              selTaus.Systematic("TES_UP").push_back(tau*1.03);
            if(passKin.GetSystematicOrValue("TES_DOWN") && passIso.GetSystematicOrValue(val))
              selTaus.Systematic("TES_DOWN").push_back(tau*0.97);
          }
          if(passKin.GetSystematicOrValue(val) && passIso.GetSystematicOrValue(val))
            selTaus.Value().push_back(tau);
        }
        else
        {
          if(passKin.GetSystematicOrValue(val) && passIso.GetSystematicOrValue(val))
            selTaus.Systematic(val).push_back(tau);
        }
      }
    }
    if(!(triggeredOn.Value()))
      continue;

    // Fill control histograms
    ValueWithSystematics<double> weightSys = (eventContent.GetDouble("weight") * eventContent.GetDouble("PUweight") * eventContent.GetDouble("xsecweight"));
    double weight = weightSys.Value();
    histMonitor.fillHisto("tauCutFlow", chTags, 0, weight);
    if(tau.isPF)
    {
      histMonitor.fillHisto("tauCutFlow", chTags, 1, weight);
      histMonitor.fillHisto("tauID", chTags, 0, weight);
      if(tau.passId(llvvTAUID::byTightCombinedIsolationDeltaBetaCorr3Hits))
      {
        histMonitor.fillHisto("tauID", chTags, 1, weight);
        if(tau.passId(llvvTAUID::decayModeFinding))
        {
          histMonitor.fillHisto("tauID", chTags, 2, weight);
          if(tau.passId(llvvTAUID::againstElectronMediumMVA5))
          {
            histMonitor.fillHisto("tauID", chTags, 3, weight);
            if(tau.passId(llvvTAUID::againstMuonTight3))
              histMonitor.fillHisto("tauID", chTags, 4, weight);
          }
        }
      }
      if(passID)
      {
        histMonitor.fillHisto("tauCutFlow", chTags, 2, weight);
        if(passQual)
        {
          histMonitor.fillHisto("tauCutFlow", chTags, 3, weight);
          if(passKin.Value())
          {
            histMonitor.fillHisto("tauCutFlow", chTags, 4, weight);
            if(passIso)
              histMonitor.fillHisto("tauCutFlow", chTags, 5, weight);
          }
        }
      }
    }
  }
  
  // Get Jets
  if(debugEvent)
    analyserCout << " Getting jets" << std::endl;
  ValueWithSystematics<llvvJetExtCollection> selJets;
  ValueWithSystematics<llvvJetExtCollection> selBJets;
  if(runSystematics && isMC)
  {
    selJets("JES_UP");
    selJets("JES_DOWN");
    selJets("JER_UP");
    selJets("JER_DOWN");

//    for(auto& syst: selTaus.Systematics())
//      selJets(syst);

    selJets.Lock();
    selBJets = selJets;
    selBJets.Lock();
  }
  if(debugEvent)
    analyserCout << "  There are " << jets.size() << " jets" << std::endl;
  for(auto& jet: jets)
  {
    if(debugEvent)
      analyserCout << "  Starting analysing a jet" << std::endl;

    // Jet ID
    bool passID = true;
    Int_t idbits = jet.idbits;
    bool passPFLoose = (idbits & 0x01);
    int fullPuId = (idbits >> 3) & 0x0f;
    bool passLooseFullPuId = ((fullPuId >> 2) & 0x01);
    passID = passLooseFullPuId;

    // Jet Kinematics
    ValueWithSystematics<bool> passKin = true;
    if(std::abs(jet.eta()) > maxJetEta)
      passKin = false;
    if(runSystematics && isMC)
    {
      passKin("JES_UP");
      passKin("JES_DOWN");
      passKin("JER_UP");
      passKin("JER_DOWN");

      if(jet.jerup <= minJetPt)
        passKin("JER_UP") = false;
      if(jet.jerdown <= minJetPt)
        passKin("JER_DOWN") = false;
      if(jet.jesup <= minJetPt)
        passKin("JES_UP") = false;
      if(jet.jesdown <= minJetPt)
        passKin("JES_DOWN") = false;
    }
    if(jet.pt() <= minJetPt)
      passKin.Value() = false;

    // B-jets
    bool isBJet = false;
//    bool hasBtagCorr = false;
    if(jet.csv > 0.679)
    {
      isBJet = true;
//      hasBtagCorr = true;
    }

    // Isolated tau
    ValueWithSystematics<bool> passIso(true);
    for(auto& tau: taus)
    {
      if(deltaR(tau, jet) < 0.4)
        passIso = false;
    }
/*    tmpLoop.clear();
    tmpLoop.push_back("Value");
    if(runSystematics && isMC)
    {
      loadSystematics(tmpLoop, selTaus);
      for(auto& syst: selTaus.Systematics())
      {
        passIso(syst);
      }
    }
    for(auto& val: tmpLoop)
    {
      for(auto& tau: selTaus.GetSystematicOrValue(val))
      {
        if(deltaR(tau, jet) < 0.4)
        {
          passIso.GetSystematicOrValue(val) = false;
          break;
        }
      }
    }// */
    // TODO: add systematics. Whenever they are added (if added) the loop before the jet loop should have the addition of the tau systematics uncommented

    if(passPFLoose && passID && static_cast<bool>(passKin) && static_cast<bool>(passIso))
    {
      tmpLoop.clear();
      tmpLoop.push_back("Value");
      if(runSystematics && isMC)
      {
        loadSystematics(tmpLoop, passIso);
      }

      for(auto& val: tmpLoop)
      {
        if(val == "Value")
        {
          if(runSystematics && isMC)
          {
            if(passKin.GetSystematicOrValue("JES_UP") && passIso.Value())
              selJets.Systematic("JES_UP").push_back(jet*(jet.jesup/jet.pt()));
            if(passKin.GetSystematicOrValue("JES_DOWN") && passIso.Value())
              selJets.Systematic("JES_DOWN").push_back(jet*(jet.jesdown/jet.pt()));
            if(passKin.GetSystematicOrValue("JER_UP") && passIso.Value())
              selJets.Systematic("JER_UP").push_back(jet*(jet.jerup/jet.pt()));
            if(passKin.GetSystematicOrValue("JER_DOWN") && passIso.Value())
              selJets.Systematic("JER_DOWN").push_back(jet*(jet.jerdown/jet.pt()));
          }
          if(passKin.Value() && passIso.Value())
            selJets.Value().push_back(jet);
        }
        else
        {
          if(passKin.GetSystematicOrValue(val) && passIso.GetSystematicOrValue(val))
            selJets.Systematic(val).push_back(jet);
        }
      }
    }
    if(passPFLoose && passID && static_cast<bool>(passKin) && isBJet && static_cast<bool>(passIso))
    {
      tmpLoop.clear();
      tmpLoop.push_back("Value");
      if(runSystematics && isMC)
      {
        loadSystematics(tmpLoop, passIso);
      }

      for(auto& val: tmpLoop)
      {
        if(val == "Value")
        {
          if(runSystematics && isMC)
          {
            if(passKin.GetSystematicOrValue("JES_UP") && passIso.Value())
              selBJets.Systematic("JES_UP").push_back(jet*(jet.jesup/jet.pt()));
            if(passKin.GetSystematicOrValue("JES_DOWN") && passIso.Value())
              selBJets.Systematic("JES_DOWN").push_back(jet*(jet.jesdown/jet.pt()));
            if(passKin.GetSystematicOrValue("JER_UP") && passIso.Value())
              selBJets.Systematic("JER_UP").push_back(jet*(jet.jerup/jet.pt()));
            if(passKin.GetSystematicOrValue("JER_DOWN") && passIso.Value())
              selBJets.Systematic("JER_DOWN").push_back(jet*(jet.jerdown/jet.pt()));
          }
          if(passKin.Value() && passIso.Value())
            selBJets.Value().push_back(jet);
        }
        else
        {
          if(passKin.GetSystematicOrValue(val) && passIso.GetSystematicOrValue(val))
            selBJets.Systematic(val).push_back(jet);
        }
      }
    }
    if(!(triggeredOn.Value()))
      continue;

    // Fill Jet control histograms
    ValueWithSystematics<double> weightSys = (eventContent.GetDouble("weight") * eventContent.GetDouble("PUweight") * eventContent.GetDouble("xsecweight"));
    double weight = weightSys.Value();
    histMonitor.fillHisto("jetCutFlow", chTags, 0, weight);
    if(passPFLoose)
    {
      histMonitor.fillHisto("jetCutFlow", chTags, 1, weight);
      if(passID)
      {
        histMonitor.fillHisto("jetCutFlow", chTags, 2, weight);
        if(passKin.Value())
        {
          histMonitor.fillHisto("jetCutFlow", chTags, 3, weight);
          if(passIso.Value())
          {
            histMonitor.fillHisto("jetCutFlow", chTags, 4, weight);
            if(isBJet)
            {
              histMonitor.fillHisto("jetCutFlow", chTags, 5, weight);
            }
          }
        }
      }
    }
  }


  if(debugEvent)
    analyserCout << " Sorting leptons, taus and jets" << std::endl;

  tmpLoop.clear();
  tmpLoop.push_back("Value");
  if(runSystematics && isMC)
  {
    loadSystematics(tmpLoop, selLeptons);
    loadSystematics(tmpLoop, selTaus);
    loadSystematics(tmpLoop, selJets);
  }

  for(auto& val: tmpLoop)
  {
    auto& lleptons = selLeptons.GetSystematicOrValue(val);
    if(lleptons.size() != 0)
      std::sort(lleptons.begin(), lleptons.end(), sort_llvvObjectByPt);

    auto& ltaus = selTaus.GetSystematicOrValue(val);
    if(ltaus.size() != 0)
      std::sort(ltaus.begin(), ltaus.end(), sort_llvvObjectByPt);

    auto& ljets = selJets.GetSystematicOrValue(val);
    if(ljets.size() != 0)
      std::sort(ljets.begin(), ljets.end(), sort_llvvObjectByPt);

//    if(val != "Value")
//      nBJet.Systematic(val) = selBJets.GetSystematicOrValue(val).size();
  }
  
  if(debugEvent)
    analyserCout << "Finished sorting" << std::endl;

  if(debugEvent)
  {
    analyserCout << " selBJets systematics:\n";
    for(auto& syst: selBJets.Systematics())
      analyserCout << "    " << syst << ": " << selBJets(syst).size() << "\n";

    analyserCout << " selJets systematics:\n";
    for(auto& syst: selJets.Systematics())
      analyserCout << "    " << syst << ": " << selJets(syst).size() << "\n";

    analyserCout << " selTaus systematics:\n";
    for(auto& syst: selTaus.Systematics())
      analyserCout << "    " << syst << ": " << selTaus(syst).size() << "\n";

    analyserCout << " selLeptons systematics:\n";
    for(auto& syst: selLeptons.Systematics())
      analyserCout << "    " << syst << ": " << selLeptons(syst).size() << "\n";
    analyserCout << std::endl;
  }

  auto& nBJet = eventContent.GetInt("nBJet");
  nBJet = selBJets.size();
  eventContent.GetInt("nJet") = selJets.size();
  eventContent.GetInt("nTau") = selTaus.size();
  eventContent.GetInt("nLep") = selLeptons.size();
  eventContent.GetInt("num") = jets.size();
  
  {
    tmpLoop.clear();
    tmpLoop.push_back("Value");
    if(runSystematics && isMC)
    {
      loadSystematics(tmpLoop, selLeptons);
    }

    ValueWithSystematics<int> nEl(0), nMu(0);
    for(auto& val: tmpLoop)
    {
      int nEl_s = 0, nMu_s = 0;
      for(auto& lep: selLeptons.GetSystematicOrValue(val))
      {
        if(abs(lep.id) == 11)
          nEl_s++;
        else
          nMu_s++;
      }
      
      if(val == "Value")
      {
        nEl.Value() = nEl_s;
        nMu.Value() = nMu_s;
      }
      else
      {
        nEl(val) = nEl_s;
        nMu(val) = nMu_s;
      }
    }

    eventContent.GetInt("nEl") = nEl;
    eventContent.GetInt("nMu") = nMu;
  }

  if(debugEvent)
    analyserCout << " Requiring an opposite sign pair" << std::endl;

  tmpLoop.clear();
  tmpLoop.push_back("Value");
  if(runSystematics && isMC)
  {
    loadSystematics(tmpLoop, selLeptons);
    loadSystematics(tmpLoop, selTaus);
  }
  
  // Opposite Sign requirements
  llvvLepton tmpLep;
  llvvTau tmpTau;
  ValueWithSystematics<llvvLepton> selectedLepton(tmpLep);
  ValueWithSystematics<llvvTau> selectedTau(tmpTau);
  auto& isOS = eventContent.GetBool("isOS");
  auto& isPromptLep = eventContent.GetBool("isPromptLep");
  auto& isPromptTau = eventContent.GetBool("isPromptTau");
  auto& isntMultilepton = eventContent.GetBool("isntMultilepton");
  auto& isETau = eventContent.GetBool("isETau");
  auto& isMuTau = eventContent.GetBool("isMuTau");
  for(auto& val: tmpLoop)
  {
    if(val != "Value")
    {
      isOS(val);
      selectedLepton(val);
      selectedTau(val);
      isPromptLep(val);
      isPromptTau(val);
      isntMultilepton(val);
      isETau(val);
      isMuTau(val);
    }
  }
  for(auto& val: tmpLoop)
  {
    double maxPtSum = 0;
    auto& leptons = selLeptons.GetSystematicOrValue(val);
    auto& taus = selTaus.GetSystematicOrValue(val);
    
    size_t lepIndex = 0;
    size_t tauIndex = 0;
    bool found = false;
    for(size_t j = 0; j < taus.size(); ++j) // The tau array is typically size 1 or 2, in this order it should normally be faster to loop (because of the break at the end of the inner loop)
    {
      for(size_t i = 0; i < leptons.size(); ++i)
      {
        if(abs(leptons[i].id) == 11) // Electron
        {
          if(leptons[i].pt() < minElPt)
            continue;
          if(std::abs(leptons[i].dZ) > maxElDz)
            continue;
          double eta = leptons[i].electronInfoRef->sceta;
          if(std::abs(eta) > maxElEta)
            continue;
          double relIso = utils::cmssw::relIso(leptons[i], eventContent.GetDouble("rho").Value());
          if(relIso > elIso)
            continue;
        }
        else
        {
          if(leptons[i].pt() < minMuPt)
            continue;
          if(std::abs(leptons[i].eta()) > maxMuEta)
            continue;
          double relIso = utils::cmssw::relIso(leptons[i], eventContent.GetDouble("rho").Value());
          if(relIso > muIso)
            continue;
          Int_t idbits = leptons[i].idbits;
          bool isTight = ((idbits >> 10) & 0x1);
          if(!isTight)
            continue;
        }
      
        double PtSum = leptons[i].pt() + taus[j].pt();
        if(PtSum > maxPtSum)
        {
          if(leptons[i].id * taus[j].id < 0) // If opposite sign
          {
            maxPtSum = PtSum;
            found = true;
            lepIndex = i;
            tauIndex = j;
          }
        }

        if(PtSum < maxPtSum)
          break;
      }
    }
    
    isOS.GetSystematicOrValue(val) = found;
    isETau.GetSystematicOrValue(val) = false;
    isMuTau.GetSystematicOrValue(val) = false;
    if(found)
    {
      selectedLepton.GetSystematicOrValue(val) = leptons[lepIndex];
      selectedTau.GetSystematicOrValue(val)    = taus[tauIndex];
      isPromptLep.GetSystematicOrValue(val) = false;
      isPromptTau.GetSystematicOrValue(val) = false;
      
      //Promptness
      if(isMC)
      {
        for(auto& genPart : gen)
        {
          if(genPart.status == 3)
          {
            if(genPart.id == leptons[lepIndex].id)
            {
              if(deltaR(leptons[lepIndex], genPart) < genMatchRCone)
              {
                isPromptLep.GetSystematicOrValue(val) = true;
              }
            }
            if(genPart.id == taus[tauIndex].id)
            {
              if(deltaR(taus[tauIndex], genPart) < genMatchRCone)
              {
                isPromptTau.GetSystematicOrValue(val) = true;
              }
            }
          }
        }
      }
    
      //Multilepton veto
      auto& isntMultilepton_ = isntMultilepton.GetSystematicOrValue(val);
      isntMultilepton_ = true;
      for(size_t i = 0; i < leptons.size(); ++i)
      {
        if(i == lepIndex)
          continue;
        if(abs(leptons[i].id) != 11 && std::abs(leptons[i].dZ) > maxMuDzVeto)
          continue;
        isntMultilepton_ = false;
        break;
      }
      
      //Channels
      if(abs(leptons[lepIndex].id) == 11)
      {
        if(val == "Value")
          chTags.push_back("etau");
        isETau.GetSystematicOrValue(val) = true;
      }
      else
      {
        if(val == "Value")
          chTags.push_back("mutau");
        isMuTau.GetSystematicOrValue(val) = true;
      }
      
      // TODO: add a check if all = etau+mutau, but only if val == "Value"
    }
  }
  
  //Lepton and tau SF
  if(debugEvent)
    analyserCout << " Getting lepton and tau scale factors" << std::endl;
  if(isMC && applyScaleFactors && static_cast<bool>(isOS))
  {
    auto& leptonSF = eventContent.GetDouble("leptonSF");
    auto& tauSF = eventContent.GetDouble("tauSF");

    if(debugEvent)
      analyserCout << "  Calling methods for the scale factors" << std::endl;
    leptonSF = leptonIdAndIsoScaleFactor(selectedLepton);
    if(debugEvent)
    {
      analyserCout << "  The leptonSF Systematics are:" << std::endl;
      analyserCout << "   - Value:" << leptonSF.Value() << std::endl;
      for(auto& syst: leptonSF.Systematics())
        analyserCout << "   - " << syst << ": " << leptonSF(syst) << std::endl;
    }
    tauSF = tauScaleFactor(selectedTau, TAU_E_ID::antiEMva5Medium);
    if(debugEvent)
    {
      analyserCout << "  The tauSF Systematics are:" << std::endl;
      analyserCout << "   - Value:" << tauSF.Value() << std::endl;
      for(auto& syst: tauSF.Systematics())
        analyserCout << "   - " << syst << ": " << tauSF(syst) << std::endl;
    }
    
    if(debugEvent)
      analyserCout << "  Multiplying scale factors" << std::endl;
    eventContent.GetDouble("weight") *= leptonSF * tauSF;
  }
  
  //Data driven stuff - TODO: make bins uncorrelated
  if(debugEvent)
    analyserCout << " Doing data-driven estimate (if applicable)" << std::endl;
  if(doDDBkg && static_cast<bool>(isOS))
  {
    auto& fakeRate   = eventContent.GetDouble("fakeRate");
    auto& promptRate = eventContent.GetDouble("promptRate");
    auto& DDweight = eventContent.GetDouble("DDweight");

    tmpLoop.clear();
    tmpLoop.push_back("Value");
    if(runSystematics)
    {
      if(isMC)
        loadSystematics(tmpLoop, selectedTau);
      loadSystematics(tmpLoop, fakeRateHist);
      loadSystematics(tmpLoop, promptRateHist);
    }

    for(auto& val: tmpLoop)
    {
      if(val != "Value")
      {
        fakeRate(val);
        promptRate(val);
        DDweight(val);
      }
      else
      {
        fakeRate("FR_UP");
        fakeRate("FR_DOWN");
        promptRate("PR_UP");
        promptRate("PR_DOWN");
        DDweight("FR_UP");
        DDweight("FR_DOWN");
        DDweight("PR_UP");
        DDweight("PR_DOWN");
      }

      auto& tau = selectedTau.GetSystematicOrValue(val);
      auto& FRhist = fakeRateHist.GetSystematicOrValue(val);
      auto& PRhist = promptRateHist.GetSystematicOrValue(val);
      auto& FR = fakeRate.GetSystematicOrValue(val);
      auto& PR = promptRate.GetSystematicOrValue(val);

      double eta = tau.eta();
      if(eta > FRhist->GetXaxis()->GetXmax())
        eta = FRhist->GetXaxis()->GetXmax();
      if(eta < FRhist->GetXaxis()->GetXmin())
        eta = FRhist->GetXaxis()->GetXmin();
      int bin = FRhist->FindBin(eta);
      std::vector<std::string> tmpLoop2;
      tmpLoop2.push_back(val);

      FR = FRhist->GetBinContent(bin);
      PR = PRhist->GetBinContent(bin);
      PR = 0.783028;
      if(val == "Value")
      {
        double tmp = FRhist->GetBinError(bin);
        fakeRate("FR_UP")     = FR + tmp;
        fakeRate("FR_DOWN")   = FR - tmp;
        tmp = PRhist->GetBinError(bin);
        promptRate("PR_UP")   = PR + 0.00379267;
        promptRate("PR_DOWN") = PR - 0.00379267;

        tmpLoop2.push_back("FR_UP");
        tmpLoop2.push_back("FR_DOWN");
        tmpLoop2.push_back("PR_UP");
        tmpLoop2.push_back("PR_DOWN");
      }

      for(auto& val2: tmpLoop2)
      {
        auto& weight = DDweight.GetSystematicOrValue(val2);
        auto& FR2 = fakeRate.GetSystematicOrValue(val2);
        auto& PR2 = promptRate.GetSystematicOrValue(val2);

        if(tau.passId(llvvTAUID::byTightCombinedIsolationDeltaBetaCorr3Hits))
        {
          weight = ((PR2 - 1) * FR2) / (PR2 - FR2);
        }
        else
        {
          weight = (PR2 * FR2) / (PR2 - FR2);
        }
      }
    }

    eventContent.GetDouble("weight") *= DDweight;
  }
  
  
  
  
  
  if(debugEvent)
    analyserCout << " Getting kinematic/topological/other variables" << std::endl;
  auto lep = selectedLepton.ToTLorentzVector();
  auto tau = selectedTau.ToTLorentzVector();
  auto met = MET;
  
  eventContent.GetDouble("deltaAlphaLepTau")      = lep.Angle(tau);
  eventContent.GetDouble("deltaRLepTau")          = selectedLepton.DeltaR<llvvTau>(selectedTau);
//  eventContent.GetDouble("deltaRLepTau")          = lep.DeltaR(tau);
  eventContent.GetDouble("deltaPhiLepTau")        = lep.DeltaPhi(tau);
  eventContent.GetDouble("deltaPhiLepTauMET")     = met.DeltaPhi(lep + tau);
  eventContent.GetDouble("minDeltaPhiMETJetPt40") = met.MinDeltaPhi<llvvJetExt>(selJets);
  eventContent.GetDouble("cosThetaLep")           = lep.CosTheta();
  eventContent.GetDouble("cosThetaTau")           = tau.CosTheta();
  eventContent.GetDouble("cosThetaMET")           = met.CosTheta();
  
  
  ValueWithSystematics<TLorentzVector> tauSystem = lep + tau;
  auto tauCS = tau;
  auto lepCS = lep;
  auto metCS = met;
  double energy = sqrtS * 500.; // sqrtS / 2 * 1000
  double mom = sqrt(energy*energy + 0.938*0.938);
  ValueWithSystematics<TLorentzVector> beam1(TLorentzVector(0, 0,  mom, energy));
  ValueWithSystematics<TLorentzVector> beam2(TLorentzVector(0, 0, -mom, energy));

  auto boost = tauSystem.BoostVector();
  boost = -boost;
  tauCS.Boost(boost);
  lepCS.Boost(boost);
  metCS.Boost(boost);
  beam1.Boost(boost);
  beam2.Boost(boost);
  
  ValueWithSystematics<TLorentzVector> SQA = beam1 - beam2;
  auto rotation = SQA.RotateTozz();
  
  SQA.Transform(rotation);
  tauCS.Transform(rotation);
  lepCS.Transform(rotation);
  metCS.Transform(rotation);
  beam1.Transform(rotation);
  beam2.Transform(rotation);
  
  eventContent.GetDouble("deltaAlphaLepTauCS")  = lepCS.Angle(tauCS);
  eventContent.GetDouble("deltaPhiLepTauCS")    = lepCS.DeltaPhi(tauCS);
  eventContent.GetDouble("deltaPhiLepTauMETCS") = metCS.DeltaPhi(lepCS + tauCS);
  eventContent.GetDouble("deltaPhiLepMETCS")    = metCS.DeltaPhi(lepCS);
  eventContent.GetDouble("cosThetaLepCS")       = lepCS.CosTheta();
  eventContent.GetDouble("cosThetaTauCS")       = tauCS.CosTheta();
  eventContent.GetDouble("cosThetaMETCS")       = metCS.CosTheta();
  
  eventContent.GetDouble("InvMass") = tauSystem.M();
  
  auto& isSVfit = eventContent.GetBool("isSVfit");
  if(doSVfit) // TODO: add svfit (remember to add the isOS and isntMultilepton checks, but then add the systematics to the eventContentSetup)
  {
    isSVfit = true;
  }
  else
  {
    isSVfit = true;
  }
  
  const ValueWithSystematics<double> unit(1);
  ValueWithSystematics<double> value;

  if(debugEvent)
    analyserCout << "  Getting MT, MT2 and deconstructed MT variables" << std::endl;
  auto& cosDeltaPhiLep = eventContent.GetDouble("cosPhiLep");
  cosDeltaPhiLep = (lep.DeltaPhi(met)).Cos();
  ValueWithSystematics<double> fac = met.Pt() * lep.Pt() * 2;
  eventContent.GetDouble("MTLep") = fac * (unit - cosDeltaPhiLep);
  eventContent.GetDouble("MTLep") = eventContent.GetDouble("MTLep").Sqrt();
  value = 80;
  eventContent.GetDouble("Q80Lep") = unit - (value*value) / fac;
  value = 100;
  eventContent.GetDouble("Q100Lep") = unit - (value*value) / fac;

  auto& cosDeltaPhiTau = eventContent.GetDouble("cosPhiTau");
  cosDeltaPhiTau = (tau.DeltaPhi(met)).Cos();
  fac = met.Pt() * tau.Pt() * 2;
  eventContent.GetDouble("MTTau") = fac * (unit - cosDeltaPhiTau);
  eventContent.GetDouble("MTTau") = eventContent.GetDouble("MTTau").Sqrt();
  value = 80;
  eventContent.GetDouble("Q80Tau") = unit - (value*value) / fac;
  value = 100;
  eventContent.GetDouble("Q100Tau") = unit - (value*value) / fac;

  eventContent.GetDouble("SumMT") = eventContent.GetDouble("MTTau") + eventContent.GetDouble("MTLep");

  eventContent.GetDouble("MT2") = computeMT2(selectedTau, selectedLepton, MET);

  if(debugEvent)
    analyserCout << " Saving some transformed variables in the event content" << std::endl;
  eventContent.GetDouble("LeptonPt") = lep.Pt();
  eventContent.GetDouble("TauPt") = tau.Pt();
  eventContent.GetDouble("METPlusLeptonPt") = met.Pt() + lep.Pt();
  eventContent.GetDouble("METPlusTauPt") = met.Pt() + tau.Pt();
  eventContent.GetDouble("absDeltaAlphaLepTau") = eventContent.GetDouble("deltaAlphaLepTau").abs();
  eventContent.GetDouble("absDeltaRLepTau") = eventContent.GetDouble("deltaRLepTau").abs();
  eventContent.GetDouble("absDeltaPhiLepTau") = eventContent.GetDouble("deltaPhiLepTau").abs();
  eventContent.GetDouble("absDeltaPhiLepTauMET") = eventContent.GetDouble("deltaPhiLepTauMET").abs();
  eventContent.GetDouble("absMinDeltaPhiMETJetPt40") = eventContent.GetDouble("minDeltaPhiMETJetPt40").abs();
  eventContent.GetDouble("absCosThetaLep") = eventContent.GetDouble("cosThetaLep").abs();
  eventContent.GetDouble("absCosThetaTau") = eventContent.GetDouble("cosThetaTau").abs();
  eventContent.GetDouble("absCosThetaMET") = eventContent.GetDouble("cosThetaMET").abs();
  eventContent.GetDouble("absDeltaAlphaLepTauCS") = eventContent.GetDouble("deltaAlphaLepTauCS").abs();
  eventContent.GetDouble("absDeltaPhiLepTauCS") = eventContent.GetDouble("deltaPhiLepTauCS").abs();
  eventContent.GetDouble("absDeltaPhiLepTauMETCS") = eventContent.GetDouble("deltaPhiLepTauMETCS").abs();
  eventContent.GetDouble("absDeltaPhiLepMETCS") = eventContent.GetDouble("deltaPhiLepMETCS").abs();
  eventContent.GetDouble("absCosThetaLepCS") = eventContent.GetDouble("cosThetaLepCS").abs();
  eventContent.GetDouble("absCosThetaTauCS") = eventContent.GetDouble("cosThetaTauCS").abs();
  eventContent.GetDouble("absCosThetaMETCS") = eventContent.GetDouble("cosThetaMETCS").abs();
  eventContent.GetDouble("absCosPhiLep") = eventContent.GetDouble("cosPhiLep").abs();
  eventContent.GetDouble("absCosPhiTau") = eventContent.GetDouble("cosPhiTau").abs();
  auto tmp = eventContent.GetDouble("InvMass");
  tmp -= 61;
  eventContent.GetDouble("absInvMassMinus61") = tmp.abs();
  tmp = eventContent.GetDouble("InvMass");
  tmp -= 60;
  eventContent.GetDouble("absInvMassMinus60") = tmp.abs();
  eventContent.GetDouble("EffMass") = met.Pt() + lep.Pt() + tau.Pt();





  if(debugEvent)
    analyserCout << " Is the event selected?" << std::endl;
  auto& selected = eventContent.GetBool("selected");
  selected = triggeredOn && isOS && isntMultilepton && (nBJet == 0) && (eventContent.GetDouble("MET") > 30);
  if(keepOnlyPromptTaus && isMC && !doDDBkg)
  {
    selected = selected && eventContent.GetBool("isPromptTau");
  }
  if(dropEvent.Value())
  {
    selected.Lock();
    selected = false;
  }
  if(debugEvent)
  {
    analyserCout << "  selected = ";
    if(selected.Value())
      analyserCout << "true";
    else
      analyserCout << "false";
    analyserCout << std::endl;
  }
}

void StauAnalyser::UserInitHistograms()
{
  histMonitor.addHistogram(new TH1D("nlep", ";# lep;Events", 10, 0, 10));
  histMonitor.addHistogram(new TH1D("nel", ";# e;Events", 10, 0, 10));
  histMonitor.addHistogram(new TH1D("nmu", ";# #mu;Events", 10, 0, 10));
  histMonitor.addHistogram(new TH1D("ntau", ";# #tau;Events", 10, 0, 10));
  histMonitor.addHistogram(new TH1D("njets", ";# jets;Events", 6, 0, 6));
  histMonitor.addHistogram(new TH1D("nbjets", ";# jets_{b};Events", 6, 0, 6));

  // Eventflow
  TH1D *eventflow = (TH1D*)histMonitor.addHistogram(new TH1D("eventflow", ";;Events", 8, 0, 8));
  eventflow->GetXaxis()->SetBinLabel(1, "HLT");
  eventflow->GetXaxis()->SetBinLabel(2, "MET > 30");
  eventflow->GetXaxis()->SetBinLabel(3, "1 lepton");
  eventflow->GetXaxis()->SetBinLabel(4, "1 #tau");
  eventflow->GetXaxis()->SetBinLabel(5, "B-veto");
  eventflow->GetXaxis()->SetBinLabel(6, "OS");
  eventflow->GetXaxis()->SetBinLabel(7, "lep veto");
  eventflow->GetXaxis()->SetBinLabel(8, "SVfit");

  // Leptons
  histMonitor.addHistogram(new TH1D("ptSelectedLep", ";p_{T}^{l};Events", 50, 0, 100));
  histMonitor.addHistogram(new TH1D("etaSelectedLep", ";#eta^{l};Events", 25, -2.6, 2.6));
  histMonitor.addHistogram(new TH1D("chargeSelectedLep", ";q^{l};Events", 5, -2, 2));
  TH1D *leptonCutFlow = (TH1D*)histMonitor.addHistogram(new TH1D("leptonCutFlow", ";;Leptons", 4, 0, 4));
  leptonCutFlow->GetXaxis()->SetBinLabel(1, "All");
  leptonCutFlow->GetXaxis()->SetBinLabel(2, "ID");
  leptonCutFlow->GetXaxis()->SetBinLabel(3, "Kin");
  leptonCutFlow->GetXaxis()->SetBinLabel(4, "Iso");

  // Taus
  histMonitor.addHistogram(new TH1D("ptSelectedTau", ";p_{T}^{#tau};Events", 50, 0, 100));
  histMonitor.addHistogram(new TH1D("ptSelectedTauExtended", ";p_{T}^{#tau};Events", 50, 0, 250));
  histMonitor.addHistogram(new TH1D("etaSelectedTau", ";#eta^{#tau};Events", 25, -2.6, 2.6));
  histMonitor.addHistogram(new TH1D("chargeSelectedTau", ";q^{#tau};Events", 5, -2, 2));
  histMonitor.addHistogram(new TH1D("dzSelectedTau", ";dz^{#tau};Events", 25, 0, 2));
  histMonitor.addHistogram(new TH1D("emfracSelectedTau", ";emf^{#tau};Events", 25, 0, 5));
  TH1D *tauCutFlow = (TH1D*)histMonitor.addHistogram(new TH1D("tauCutFlow", ";;#tau", 6, 0, 6));
  tauCutFlow->GetXaxis()->SetBinLabel(1, "All");
  tauCutFlow->GetXaxis()->SetBinLabel(2, "PF");
  tauCutFlow->GetXaxis()->SetBinLabel(3, "ID");
  tauCutFlow->GetXaxis()->SetBinLabel(4, "Quality");
  tauCutFlow->GetXaxis()->SetBinLabel(5, "Kin");
  tauCutFlow->GetXaxis()->SetBinLabel(6, "Iso");
  TH1D *tauID = (TH1D*)histMonitor.addHistogram(new TH1D("tauID", ";;#tau", 5, 0, 5));
  tauID->GetXaxis()->SetBinLabel(1, "All");
  tauID->GetXaxis()->SetBinLabel(2, "Medium comb iso");
  tauID->GetXaxis()->SetBinLabel(3, "Decay mode");
  tauID->GetXaxis()->SetBinLabel(4, "Not e");
  tauID->GetXaxis()->SetBinLabel(5, "Not #mu");

  // Jets
  histMonitor.addHistogram(new TH1D("jetleadpt", ";p_{T}^{jet};Events", 25, 0, 500));
  histMonitor.addHistogram(new TH1D("jetleadeta", ";#eta^{jet};Events", 50, -5, 5));
  histMonitor.addHistogram(new TH1D("jetcsv", ";csv;jets", 25, 0, 1));
  TH1D *jetCutFlow = (TH1D*)histMonitor.addHistogram(new TH1D("jetCutFlow", ";;jets", 6, 0, 6));
  jetCutFlow->GetXaxis()->SetBinLabel(1, "All");
  jetCutFlow->GetXaxis()->SetBinLabel(2, "PF Loose");
  jetCutFlow->GetXaxis()->SetBinLabel(3, "ID");
  jetCutFlow->GetXaxis()->SetBinLabel(4, "Kin");
  jetCutFlow->GetXaxis()->SetBinLabel(5, "Iso");
  jetCutFlow->GetXaxis()->SetBinLabel(6, "B-jet");

  // MT
  histMonitor.addHistogram(new TH1D("MT", ";MT [GeV];Events", 25, 0, 200));
  histMonitor.addHistogram(new TH1D("MTTau", ";MT(#tau) [GeV];Events", 25, 0, 200));
  histMonitor.addHistogram(new TH1D("SumMT", ";SumMT [GeV];Events", 25, 0, 200));

  // Deconstructed MT: https://indico.cern.ch/event/344807/
  histMonitor.addHistogram(new TH1D("Q80", ";Q_{80};Events", 30, -2, 1));
  histMonitor.addHistogram(new TH1D("Q100", ";Q_{100};Events", 30, -2, 1));
  histMonitor.addHistogram(new TH1D("cosPhi", ";cos#Phi;Events", 30, -1, 1));
  histMonitor.addHistogram(new TH1D("Q80Tau", ";Q_{80};Events", 30, -2, 1));
  histMonitor.addHistogram(new TH1D("Q100Tau", ";Q_{100};Events", 30, -2, 1));
  histMonitor.addHistogram(new TH1D("cosPhiTau", ";cos#Phi;Events", 30, -1, 1));

  // MT2
  histMonitor.addHistogram(new TH1D("MT2", ";M_{T2} [GeV];Events", 25, 0, 500));

  // SVFit Mass
  if(doSVfit)
    histMonitor.addHistogram(new TH1D("SVFitMass", ";M_{SVFit};Events", 50, 0, 500));

  // Invariant Mass
  histMonitor.addHistogram(new TH1D("InvMass", ";M_{l-#tau};Events", 50, 0, 500));

  // Angles
  histMonitor.addHistogram(new TH1D("deltaAlphaLepTau", ";#Delta#alpha_{l-#tau}(Lab);Events", 30, 0, TMath::Pi()));
  histMonitor.addHistogram(new TH1D("deltaRLepTau", ";#Delta R_{l-#tau}(Lab);Events", 40, 0, 8));
  histMonitor.addHistogram(new TH1D("deltaPhiLepTauMET", ";#Delta#phi_{l#tau-MET}(Lab);Events", 30, -TMath::Pi(), TMath::Pi()));
  histMonitor.addHistogram(new TH1D("deltaPhiLepTau", ";#Delta#phi_{l-#tau}(Lab);Events", 30, -TMath::Pi(), TMath::Pi()));
  histMonitor.addHistogram(new TH1D("cosThetaTau", ";cos#theta_{#tau}(Lab);Events", 30, -1, 1));
  histMonitor.addHistogram(new TH1D("cosThetaLep", ";cos#theta_{l}(Lab);Events", 30, -1, 1));
  histMonitor.addHistogram(new TH1D("deltaPhiLepMETCS", ";#Delta#phi_{l-MET}(CS);Events", 30, -TMath::Pi(), TMath::Pi()));
  histMonitor.addHistogram(new TH1D("cosThetaCS", ";cos#theta(CS);Events", 30, -1, 1));
  histMonitor.addHistogram(new TH1D("minDeltaPhiMETJetPt40", ";min(#Delta#phi_{MET-Jet40});Events", 20, -TMath::Pi(), TMath::Pi()));

  // 2D variables
  histMonitor.addHistogram(new TH2D("metVsPtl", ";p_{T}(l);MET", 50, 0, 100, 25, 0, 200));
  histMonitor.addHistogram(new TH2D("metVsPtTau", ";p_{T}(#tau);MET", 50, 0, 100, 25, 0, 200));
  histMonitor.addHistogram(new TH2D("metPtVsmetEt", ";met.Et();met.pt()", 25, 0, 200, 25, 0, 200));
  // Deconstructed MT 2D Plots:
  histMonitor.addHistogram(new TH2D("Q80VsCosPhi", ";cos#Phi;Q_{80}", 20, -1, 1, 20, -2, 1));
  histMonitor.addHistogram(new TH2D("Q100VsCosPhi", ";cos#Phi;Q_{100}", 20, -1, 1, 20, -2, 1));
  histMonitor.addHistogram(new TH2D("Q80VsCosPhiTau", ";cos#Phi;Q_{80}", 20, -1, 1, 20, -2, 1));
  histMonitor.addHistogram(new TH2D("Q100VsCosPhiTau", ";cos#Phi;Q_{100}", 20, -1, 1, 20, -2, 1));

  if(debug)
    std::cout << "Finished StauAnalyser::UserInitHistograms()" << std::endl;

  return;
}

void StauAnalyser::UserFillHistograms()
{
  auto& weight = eventContent.GetDouble("weight").Value();
//  auto& puWeight = eventContent.GetDouble("PUweight").Value();
  auto& selected = eventContent.GetBool("selected").Value();
  auto& dropEvent = eventContent.GetBool("dropEvent").Value();
  
  bool plotEvent = true;
  if(isStauStau)
  {
    plotEvent = false;
    if(eventContent.GetDouble("stauMass").Value() == stauMtoPlot && eventContent.GetDouble("neutralinoMass").Value() == neutralinoMtoPlot)
      plotEvent = true;
  }


  // Eventflow
  if(!dropEvent && plotEvent && eventContent.GetBool("triggeredOn").Value() && (!(keepOnlyPromptTaus && isMC && !doDDBkg) || eventContent.GetBool("isPromptTau").Value()))
  {
    histMonitor.fillHisto("eventflow", chTags, 0, weight);
    if(eventContent.GetDouble("MET").Value() > 30)
    {
      histMonitor.fillHisto("eventflow", chTags, 1, weight);
      if(eventContent.GetInt("nLep").Value() >= 1)
      {
        histMonitor.fillHisto("eventflow", chTags, 2, weight);
        if(eventContent.GetInt("nTau").Value() >= 1)
        {
          histMonitor.fillHisto("eventflow", chTags, 3, weight);
          if(eventContent.GetInt("nBJet").Value() == 0)
          {
            histMonitor.fillHisto("eventflow", chTags, 4, weight);
            if(eventContent.GetBool("isOS").Value())
            {
              histMonitor.fillHisto("eventflow", chTags, 5, weight);
              if(eventContent.GetBool("isntMultilepton").Value())
              {
                histMonitor.fillHisto("eventflow", chTags, 6, weight);
                if(eventContent.GetBool("isSVfit").Value())
                {
                  histMonitor.fillHisto("eventflow", chTags, 7, weight);
                }
              }
            }
          }
        }
      }
    }
  }


  if(selected && plotEvent)
  {
    histMonitor.fillHisto("nlep", chTags, eventContent.GetInt("nLep").Value(), weight);
    histMonitor.fillHisto("nel", chTags, eventContent.GetInt("nEl").Value(), weight);
    histMonitor.fillHisto("nmu", chTags, eventContent.GetInt("nMu").Value(), weight);
    histMonitor.fillHisto("ntau", chTags, eventContent.GetInt("nTau").Value(), weight);
    histMonitor.fillHisto("njets", chTags, eventContent.GetInt("nJet").Value(), weight);
    histMonitor.fillHisto("nbjets", chTags, eventContent.GetInt("nBJet").Value(), weight);

    histMonitor.fillHisto("MET", chTags, eventContent.GetDouble("MET").Value(), weight);
    
    // Selected Lepton
    histMonitor.fillHisto("ptSelectedLep", chTags, eventContent.GetDouble("LeptonPt").Value(), weight);
//    histMonitor.fillHisto("etaSelectedLep", chTags, eventContent.GetInt("nBJet").Value(), weight);
//    histMonitor.fillHisto("chargeSelectedLep", chTags, eventContent.GetInt("nBJet").Value(), weight);

    // Selected Tau
    histMonitor.fillHisto("ptSelectedTau", chTags, eventContent.GetDouble("TauPt").Value(), weight);
    histMonitor.fillHisto("ptSelectedTauExtended", chTags, eventContent.GetDouble("TauPt").Value(), weight);
//    histMonitor.fillHisto("etaSelectedTau", chTags, eventContent.GetDouble("TauPt").Value(), weight);
//    histMonitor.fillHisto("chargeSelectedTau", chTags, eventContent.GetDouble("TauPt").Value(), weight);
//    histMonitor.fillHisto("dzSelectedTau", chTags, eventContent.GetDouble("TauPt").Value(), weight);
//    histMonitor.fillHisto("emfracSelectedTau", chTags, eventContent.GetDouble("TauPt").Value(), weight);

    // Jets
//    histMonitor.fillHisto("jetleadpt", chTags, eventContent.GetDouble("TauPt").Value(), weight);
//    histMonitor.fillHisto("jetleadeta", chTags, eventContent.GetDouble("TauPt").Value(), weight);
//    histMonitor.fillHisto("jetcsv", chTags, eventContent.GetDouble("TauPt").Value(), weight);

    // MT
    histMonitor.fillHisto("MT", chTags, eventContent.GetDouble("MTLep").Value(), weight);
    histMonitor.fillHisto("MTTau", chTags, eventContent.GetDouble("MTTau").Value(), weight);
    histMonitor.fillHisto("SumMT", chTags, eventContent.GetDouble("SumMT").Value(), weight);

    // Deconstructed MT: https://indico.cern.ch/event/344807/
    histMonitor.fillHisto("Q80", chTags, eventContent.GetDouble("Q80Lep").Value(), weight);
    histMonitor.fillHisto("Q100", chTags, eventContent.GetDouble("Q100Lep").Value(), weight);
    histMonitor.fillHisto("cosPhi", chTags, eventContent.GetDouble("cosPhiLep").Value(), weight);
    histMonitor.fillHisto("Q80Tau", chTags, eventContent.GetDouble("Q80Tau").Value(), weight);
    histMonitor.fillHisto("Q100Tau", chTags, eventContent.GetDouble("Q100Tau").Value(), weight);
    histMonitor.fillHisto("cosPhiTau", chTags, eventContent.GetDouble("cosPhiTau").Value(), weight);

    // MT2
    histMonitor.fillHisto("MT2", chTags, eventContent.GetDouble("MT2").Value(), weight);

    // SVFit Mass
    if(doSVfit)
      histMonitor.fillHisto("SVFitMass", chTags, eventContent.GetDouble("SVFitMass").Value(), weight);

    // Invariant Mass
    histMonitor.fillHisto("InvMass", chTags, eventContent.GetDouble("InvMass").Value(), weight);
  
  
/*  eventContent.AddDouble("cosThetaMET", -999);
  eventContent.AddDouble("deltaAlphaLepTauCS", -999);
  eventContent.AddDouble("deltaPhiLepTauCS", -999);
  eventContent.AddDouble("deltaPhiLepTauMETCS", -999);
  eventContent.AddDouble("cosThetaMETCS", -999);

  eventContent.AddDouble("METPlusLeptonPt", -999);
  eventContent.AddDouble("METPlusTauPt", -999);
  eventContent.AddDouble("absDeltaAlphaLepTau", -999);
  eventContent.AddDouble("absDeltaRLepTau", -999);
  eventContent.AddDouble("absDeltaPhiLepTau", -999);
  eventContent.AddDouble("absDeltaPhiLepTauMET", -999);
  eventContent.AddDouble("absMinDeltaPhiMETJetPt40", -999);
  eventContent.AddDouble("absCosThetaLep", -999);
  eventContent.AddDouble("absCosThetaTau", -999);
  eventContent.AddDouble("absCosThetaMET", -999);
  eventContent.AddDouble("absDeltaAlphaLepTauCS", -999);
  eventContent.AddDouble("absDeltaPhiLepTauCS", -999);
  eventContent.AddDouble("absDeltaPhiLepTauMETCS", -999);
  eventContent.AddDouble("absDeltaPhiLepMETCS", -999);
  eventContent.AddDouble("absCosThetaLepCS", -999);
  eventContent.AddDouble("absCosThetaTauCS", -999);
  eventContent.AddDouble("absCosThetaMETCS", -999);
  eventContent.AddDouble("absCosPhiLep", -999);
  eventContent.AddDouble("absCosPhiTau", -999);
  eventContent.AddDouble("absInvMassMinus61", -999);
  eventContent.AddDouble("absInvMassMinus60", -999);
  eventContent.AddDouble("EffMass", -999);// */

    // Angles
    histMonitor.fillHisto("deltaAlphaLepTau", chTags, eventContent.GetDouble("deltaAlphaLepTau").Value(), weight);
    histMonitor.fillHisto("deltaRLepTau", chTags, eventContent.GetDouble("deltaRLepTau").Value(), weight);
    histMonitor.fillHisto("deltaPhiLepTauMET", chTags, eventContent.GetDouble("deltaPhiLepTauMET").Value(), weight);
    histMonitor.fillHisto("deltaPhiLepTau", chTags, eventContent.GetDouble("deltaPhiLepTau").Value(), weight);
    histMonitor.fillHisto("cosThetaTau", chTags, eventContent.GetDouble("cosThetaTau").Value(), weight);
    histMonitor.fillHisto("cosThetaLep", chTags, eventContent.GetDouble("cosThetaLep").Value(), weight);
    histMonitor.fillHisto("deltaPhiLepMETCS", chTags, eventContent.GetDouble("deltaPhiLepMETCS").Value(), weight);
    histMonitor.fillHisto("cosThetaCS", chTags, eventContent.GetDouble("cosThetaLepCS").Value(), weight);
    histMonitor.fillHisto("minDeltaPhiMETJetPt40", chTags, eventContent.GetDouble("minDeltaPhiMETJetPt40").Value(), weight);

    // 2D variables
    histMonitor.fillHisto("metVsPtl", chTags, eventContent.GetDouble("LeptonPt").Value(), eventContent.GetDouble("MET").Value(), weight);
    histMonitor.fillHisto("metVsPtTau", chTags, eventContent.GetDouble("TauPt").Value(), eventContent.GetDouble("MET").Value(), weight);
//    histMonitor.fillHisto("metPtVsmetEt", chTags, eventContent.GetDouble("LeptonPt").Value(), eventContent.GetDouble("MET").Value(), weight);
    // Deconstructed MT 2D Plots:
    histMonitor.fillHisto("Q80VsCosPhi", chTags, eventContent.GetDouble("cosPhiLep").Value(), eventContent.GetDouble("Q80Lep").Value(), weight);
    histMonitor.fillHisto("Q100VsCosPhi", chTags, eventContent.GetDouble("cosPhiLep").Value(), eventContent.GetDouble("Q100Lep").Value(), weight);
    histMonitor.fillHisto("Q80VsCosPhiTau", chTags, eventContent.GetDouble("cosPhiTau").Value(), eventContent.GetDouble("Q80Tau").Value(), weight);
    histMonitor.fillHisto("Q100VsCosPhiTau", chTags, eventContent.GetDouble("cosPhiTau").Value(), eventContent.GetDouble("Q100Tau").Value(), weight);
  }
}

void StauAnalyser::UserEventContentSetup()
{
  eventContent.AddDouble("stauMass", 0);
  eventContent.AddDouble("neutralinoMass", 0);

  eventContent.AddBool("triggeredOn", false);
  eventContent.AddBool("dropEvent", false);
  eventContent.AddBool("isSVfit", false);

  auto& triggerSF = eventContent.AddDouble("triggerSF", 1);
  triggerSF.AddMetadata("eventlist", "true");
  if(runSystematics && isMC)
  {
/*    triggerSF("etauTrig_UP");
    triggerSF("etauTrig_DOWN");
    triggerSF("mutauTrig_UP");
    triggerSF("mutauTrig_DOWN");// */
    triggerSF.Lock();
  }
  auto& leptonSF = eventContent.AddDouble("leptonSF", 1);
  leptonSF.AddMetadata("eventlist", "true");
  if(runSystematics && isMC)
  {
    leptonSF.Systematic("elID_UP");
    leptonSF.Systematic("elID_DOWN");
    leptonSF.Systematic("muID_UP");
    leptonSF.Systematic("muID_DOWN");
    leptonSF.Systematic("elISO_UP");
    leptonSF.Systematic("elISO_DOWN");
    leptonSF.Systematic("muISO_UP");
    leptonSF.Systematic("muISO_DOWN");
    leptonSF.Systematic("TES_UP");
    leptonSF.Systematic("TES_DOWN");
    leptonSF.Systematic("LES_UP");
    leptonSF.Systematic("LES_DOWN");

    leptonSF.Lock();
  }
  auto& tauSF = eventContent.AddDouble("tauSF", 1);
  tauSF.AddMetadata("eventlist", "true");
  if(runSystematics && isMC)
  {
    tauSF.Systematic("tauID_UP");
    tauSF.Systematic("tauID_DOWN");
    tauSF.Systematic("tauFromESF_UP");
    tauSF.Systematic("tauFromESF_DOWN");
    tauSF.Systematic("tauFromMu_UP");
    tauSF.Systematic("tauFromMu_DOWN");
    tauSF.Systematic("tauFromJet_UP");
    tauSF.Systematic("tauFromJet_DOWN");
    tauSF.Systematic("TES_UP");
    tauSF.Systematic("TES_DOWN");
    tauSF.Systematic("LES_UP");
    tauSF.Systematic("LES_DOWN");

    tauSF.Lock();
  }
  auto& weight = eventContent.GetDouble("weight");
  if(runSystematics && isMC)
  {
    weight.Systematic("elID_UP");
    weight.Systematic("elID_DOWN");
    weight.Systematic("muID_UP");
    weight.Systematic("muID_DOWN");
    weight.Systematic("elISO_UP");
    weight.Systematic("elISO_DOWN");
    weight.Systematic("muISO_UP");
    weight.Systematic("muISO_DOWN");
    weight.Systematic("tauID_UP");
    weight.Systematic("tauID_DOWN");
    weight.Systematic("tauFromESF_UP");
    weight.Systematic("tauFromESF_DOWN");
    weight.Systematic("tauFromMu_UP");
    weight.Systematic("tauFromMu_DOWN");
    weight.Systematic("tauFromJet_UP");
    weight.Systematic("tauFromJet_DOWN");
    weight.Systematic("TES_UP");
    weight.Systematic("TES_DOWN");
    weight.Systematic("LES_UP");
    weight.Systematic("LES_DOWN");
  }
  if(runSystematics && doDDBkg)
  {
    weight("FR_UP");
    weight("FR_DOWN");
    weight("PR_UP");
    weight("PR_DOWN");
  }
  auto& DDweight   = eventContent.AddDouble("DDweight", 1);
  auto& fakeRate   = eventContent.AddDouble("fakeRate", 1);
  auto& promptRate = eventContent.AddDouble("promptRate", 1);
  if(runSystematics && doDDBkg)
  {
    fakeRate("FR_UP");
    fakeRate("FR_DOWN");
    promptRate("PR_UP");
    promptRate("PR_DOWN");
    DDweight("FR_UP");
    DDweight("FR_DOWN");
    DDweight("PR_UP");
    DDweight("PR_DOWN");
    if(isMC)
    {
      fakeRate("TES_UP");
      fakeRate("TES_DOWN");
      promptRate("TES_UP");
      promptRate("TES_DOWN");
      DDweight("TES_UP");
      DDweight("TES_DOWN");
    }

    fakeRate.Lock();
    promptRate.Lock();
    DDweight.AddMetadata("eventlist", "true");
    DDweight.Lock();
  }

  auto& num = eventContent.AddInt("num", 0);
  num.AddMetadata("eventlist", "true");
  num.AddMetadata("eventtree", "false");

  auto& nLep = eventContent.AddInt("nLep", 0);
  nLep.AddMetadata("eventlist", "true");
//  nLep.AddMetadata("eventtree", "false");

  auto& nEl = eventContent.AddInt("nEl", 0);
  nEl.AddMetadata("eventlist", "true");
//  nEl.AddMetadata("eventtree", "false");

  auto& nMu = eventContent.AddInt("nMu", 0);
  nMu.AddMetadata("eventlist", "true");
//  nMu.AddMetadata("eventtree", "false");

  auto& nTau = eventContent.AddInt("nTau", 0);
  nTau.AddMetadata("eventlist", "true");
//  nTau.AddMetadata("eventtree", "false");

  auto& nJet = eventContent.AddInt("nJet", 0);
  nJet.AddMetadata("eventlist", "true");
//  nJet.AddMetadata("eventtree", "false");

  auto& nBJet = eventContent.AddInt("nBJet", 0);
  nBJet.AddMetadata("eventlist", "true");
  eventContent.AddBool("isOS", false);
  eventContent.AddBool("isPromptLep", false);
  eventContent.AddBool("isPromptTau", false);
  eventContent.AddBool("isntMultilepton", false);
  eventContent.AddBool("isETau", false);
  eventContent.AddBool("isMuTau", false);
  
  
  eventContent.AddDouble("deltaAlphaLepTau", -999);
  eventContent.AddDouble("deltaRLepTau", -999);
  eventContent.AddDouble("deltaPhiLepTau", -999);
  eventContent.AddDouble("deltaPhiLepTauMET", -999);
  eventContent.AddDouble("minDeltaPhiMETJetPt40", -999);
  eventContent.AddDouble("cosThetaLep", -999);
  eventContent.AddDouble("cosThetaTau", -999);
  eventContent.AddDouble("cosThetaMET", -999);
  eventContent.AddDouble("deltaAlphaLepTauCS", -999);
  eventContent.AddDouble("deltaPhiLepTauCS", -999);
  eventContent.AddDouble("deltaPhiLepTauMETCS", -999);
  eventContent.AddDouble("deltaPhiLepMETCS", -999);
  eventContent.AddDouble("cosThetaLepCS", -999);
  eventContent.AddDouble("cosThetaTauCS", -999);
  eventContent.AddDouble("cosThetaMETCS", -999);
  eventContent.AddDouble("InvMass", -999);
  eventContent.AddDouble("cosPhiLep", -999);
  eventContent.AddDouble("cosPhiTau", -999);
  auto& MTLep = eventContent.AddDouble("MTLep", -999);
  auto& MTTau = eventContent.AddDouble("MTTau", -999);
  eventContent.AddDouble("Q80Lep", -999);
  eventContent.AddDouble("Q80Tau", -999);
  eventContent.AddDouble("Q100Lep", -999);
  eventContent.AddDouble("Q100Tau", -999);
  auto& SumMT = eventContent.AddDouble("SumMT", -999);
  auto& MT2 = eventContent.AddDouble("MT2", -999);

  MTLep.AddMetadata("eventlist", "true");
  MTTau.AddMetadata("eventlist", "true");
  SumMT.AddMetadata("eventlist", "true");
  MT2.AddMetadata("eventlist", "true");

  eventContent.AddDouble("LeptonPt", -999);
  eventContent.AddDouble("TauPt", -999);
  eventContent.AddDouble("METPlusLeptonPt", -999);
  eventContent.AddDouble("METPlusTauPt", -999);
  eventContent.AddDouble("absDeltaAlphaLepTau", -999);
  eventContent.AddDouble("absDeltaRLepTau", -999);
  eventContent.AddDouble("absDeltaPhiLepTau", -999);
  eventContent.AddDouble("absDeltaPhiLepTauMET", -999);
  eventContent.AddDouble("absMinDeltaPhiMETJetPt40", -999);
  eventContent.AddDouble("absCosThetaLep", -999);
  eventContent.AddDouble("absCosThetaTau", -999);
  eventContent.AddDouble("absCosThetaMET", -999);
  eventContent.AddDouble("absDeltaAlphaLepTauCS", -999);
  eventContent.AddDouble("absDeltaPhiLepTauCS", -999);
  eventContent.AddDouble("absDeltaPhiLepTauMETCS", -999);
  eventContent.AddDouble("absDeltaPhiLepMETCS", -999);
  eventContent.AddDouble("absCosThetaLepCS", -999);
  eventContent.AddDouble("absCosThetaTauCS", -999);
  eventContent.AddDouble("absCosThetaMETCS", -999);
  eventContent.AddDouble("absCosPhiLep", -999);
  eventContent.AddDouble("absCosPhiTau", -999);
  eventContent.AddDouble("absInvMassMinus61", -999);
  eventContent.AddDouble("absInvMassMinus60", -999);
  eventContent.AddDouble("EffMass", -999);
  if(doSVfit)
    eventContent.AddDouble("SVFitMass", -999);

  if(debug)
    std::cout << "Finished StauAnalyser::UserEventContentSetup()" << std::endl;

  return;
}

ValueWithSystematics<double> StauAnalyser::LeptonTauTriggerScaleFactor(llvvLepton& lepton, llvvTau& tau)
{
  ValueWithSystematics<double> scaleFactor(1);
  double m0[2], sigma[2], alpha[2], n[2], norm[2]; // Index 0 - Data; Index 1 - MC
  double pt, eta;
  
  if(abs(lepton.id) == 11) // eTau channel
  {
    // Electron leg
    eta = lepton.electronInfoRef->sceta;
    pt  = lepton.pt();
    if(std::abs(eta) < 1.479) // In barrel
    {
      m0[0]    = 22.9704;
      m0[1]    = 21.7243;
      sigma[0] = 1.0258;
      sigma[1] = 0.619015;
      alpha[0] = 1.26889;
      alpha[1] = 0.739301;
      n[0]     = 1.31024;
      n[1]     = 1.34903;
      norm[0]  = 1.06409;
      norm[1]  = 1.02594;
    }
    else // In endcap
    {
      m0[0] = 21.9816;
      m0[1] = 22.1217;
      sigma[0] = 1.40993;
      sigma[1] = 1.34054;
      alpha[0] = 0.978597;
      alpha[1] = 1.8885;
      n[0] = 2.33144;
      n[1] = 1.01855;
      norm[0] = 0.937552;
      norm[1] = 4.7241;
    }

    double electronSF = 1;
    if(pt >= 20) // Do not apply for electrons with pt below threshold
    {
      double electronDataEff = Efficiency(pt, m0[0], sigma[0], alpha[0], n[0], norm[0]);
      double electronMCEff   = Efficiency(pt, m0[1], sigma[1], alpha[1], n[1], norm[1]);
      electronSF = electronDataEff/electronMCEff;
    }

    // Tau leg
    eta = tau.eta();
    pt  = tau.pt();
    if(std::abs(eta) < 1.5) // In barrel
    {
      m0[0]    = 18.538229;
      m0[1]    = 18.605055;
      sigma[0] = 0.651562;
      sigma[1] = 0.264062;
      alpha[0] = 0.324869;
      alpha[1] = 0.139561;
      n[0]     = 13.099048;
      n[1]     = 4.792849;
      norm[0]  = 0.902365;
      norm[1]  = 0.915035;
    }
    else // In endcap
    {
      m0[0]    = 18.756548;
      m0[1]    = 18.557810;
      sigma[0] = 0.230732;
      sigma[1] = 0.280908;
      alpha[0] = 0.142859;
      alpha[1] = 0.119282;
      n[0]     = 3.358497;
      n[1]     = 17.749043;
      norm[0]  = 0.851919;
      norm[1]  = 0.865756;
    }

    double tauSF = 1;
    if(pt >= 20)
    {
      double tauDataEff = Efficiency(pt, m0[0], sigma[0], alpha[0], n[0], norm[0]);
      double tauMCEff   = Efficiency(pt, m0[1], sigma[1], alpha[1], n[1], norm[1]);
      tauSF = tauDataEff/tauMCEff;
    }

    scaleFactor = electronSF * tauSF;
  }
  else // muTau channel
  {
    // Muon leg
    eta = lepton.eta();
    pt  = lepton.pt();
    if(eta < -1.2)
    {
      m0[0]    = 15.9977;
      m0[1]    = 16.0051;
      sigma[0] = 7.64004e-05;
      sigma[1] = 2.45144e-05;
      alpha[0] = 6.4951e-08;
      alpha[1] = 4.3335e-09;
      n[0]     = 1.57403;
      n[1]     = 1.66134;
      norm[0]  = 0.865325;
      norm[1]  = 0.87045;
    }
    else if(eta < -0.8)
    {
      m0[0]    = 17.3974;
      m0[1]    = 17.3135;
      sigma[0] = 0.804001;
      sigma[1] = 0.747636;
      alpha[0] = 1.47145;
      alpha[1] = 1.21803;
      n[0]     = 1.24295;
      n[1]     = 1.40611;
      norm[0]  = 0.928198;
      norm[1]  = 0.934983;
    }
    else if(eta < 0)
    {
      m0[0]    = 16.4307;
      m0[1]    = 15.9556;
      sigma[0] = 0.226312;
      sigma[1] = 0.0236127;
      alpha[0] = 0.265553;
      alpha[1] = 0.00589832;
      n[0]     = 1.55756;
      n[1]     = 1.75409;
      norm[0]  = 0.974462;
      norm[1]  = 0.981338;
    }
    else if(eta < 0.8)
    {
      m0[0]    = 17.313;
      m0[1]    = 15.9289;
      sigma[0] = 0.662731;
      sigma[1] = 0.0271317;
      alpha[0] = 1.3412;
      alpha[1] = 0.00448573;
      n[0]     = 1.05778;
      n[1]     = 1.92101;
      norm[0]  = 1.26624;
      norm[1]  = 0.978625;
    }
    else if(eta < 1.2)
    {
      m0[0]    = 16.9966;
      m0[1]    = 16.5678;
      sigma[0] = 0.550532;
      sigma[1] = 0.328333;
      alpha[0] = 0.807863;
      alpha[1] = 0.354533;
      n[0]     = 1.55402;
      n[1]     = 1.67085;
      norm[0]  = 0.885134;
      norm[1]  = 0.916992;
    }
    else
    {
      m0[0]    = 15.9962;
      m0[1]    = 15.997;
      sigma[0] = 0.000106195;
      sigma[1] = 7.90069e-05;
      alpha[0] = 4.95058e-08;
      alpha[1] = 4.40036e-08;
      n[0]     = 1.9991;
      n[1]     = 1.66272;
      norm[0]  = 0.851294;
      norm[1]  = 0.884502;
    }

    double muonSF = 1;
    if(pt >= 17)
    {
      double muonDataEff = Efficiency(pt, m0[0], sigma[0], alpha[0], n[0], norm[0]);
      double muonMCEff   = Efficiency(pt, m0[1], sigma[1], alpha[1], n[1], norm[1]);
      muonSF = muonDataEff/muonMCEff;
    }

    // Tau leg
    eta = tau.eta();
    pt  = tau.pt();
    if(std::abs(eta) < 1.5) // In barrel
    {
      m0[0]    = 18.604910;
      m0[1]    = 18.532997;
      sigma[0] = 0.276042;
      sigma[1] = 1.027880;
      alpha[0] = 0.137039;
      alpha[1] = 2.262950;
      n[0]     = 2.698437;
      n[1]     = 1.003322;
      norm[0]  = 0.940721;
      norm[1]  = 5.297292;
    }
    else // In endcap
    {
      m0[0]    = 18.701715;
      m0[1]    = 18.212782;
      sigma[0] = 0.216523;
      sigma[1] = 0.338119;
      alpha[0] = 0.148111;
      alpha[1] = 0.122828;
      n[0]     = 2.245081;
      n[1]     = 12.577926;
      norm[0]  = 0.895320;
      norm[1]  = 0.893975;
    }

    double tauSF = 1;
    if(pt >= 20)
    {
      double tauDataEff = Efficiency(pt, m0[0], sigma[0], alpha[0], n[0], norm[0]);
      double tauMCEff   = Efficiency(pt, m0[1], sigma[1], alpha[1], n[1], norm[1]);
      tauSF = tauDataEff/tauMCEff;
    }

    scaleFactor = muonSF * tauSF;
  }
  
  return scaleFactor;
}

// Following function from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2012#ETau_MuTau_trigger_turn_on_Joshu
// it parametrizes a trigger efficiency turn on curve. m is the pT of the object
double StauAnalyser::Efficiency(double m, double m0, double sigma, double alpha, double n, double norm)
{
  const double sqrtPiOver2 = 1.2533141373;
  const double sqrt2 = 1.4142135624;
  double sig = fabs((double) sigma);
  double t = (m - m0)/sig;
  if(alpha < 0)
    t = -t;
  double absAlpha = fabs(alpha/sig);
  double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
  double b = absAlpha - n/absAlpha;
  double ApproxErf;
  double arg = absAlpha / sqrt2;
  if (arg > 5.) ApproxErf = 1;
  else if (arg < -5.) ApproxErf = -1;
  else ApproxErf = TMath::Erf(arg);
  double leftArea = (1 + ApproxErf) * sqrtPiOver2;
  double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
  double area = leftArea + rightArea;
  if( t <= absAlpha )
  {
    arg = t / sqrt2;
    if(arg > 5.) ApproxErf = 1;
    else if (arg < -5.) ApproxErf = -1;
    else ApproxErf = TMath::Erf(arg);
    return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
  }
  else
  {
    return norm * (leftArea + a * (1/TMath::Power(t-b,n-1) - 1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area;
  }
}

ValueWithSystematics<double> StauAnalyser::StauCrossSec()
{
  // Cross sections from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections8TeVsleptonsleptonCMS
  // For MStau < 100, points taken from  http://arxiv.org/abs/1204.2379
  //   Mstau == 100 => 0.1
  //   Mstau == 125 => 0.05
  //   Mstau == 145 => 0.03
  //   Mstau == 195 => 0.01
  //   Mstau == 240 => 0.005
  //   Mstau == 275 => 0.003
  //   Mstau == 300 => 0.002
  //   Mstau == 360 => 0.001
  //   Mstau == 425 => 0.0005

  ValueWithSystematics<double> stauMass = eventContent.GetDouble("stauMass");
  double stauM = stauMass.Value();
  ValueWithSystematics<double> retVal;
  
  // New code
  Int_t bin = xSecHist.Value()->FindBin(stauM);
  retVal = xSecHist.Value()->GetBinContent(bin);
  if(runSystematics)
  {
    retVal("xsec_UP") = xSecHist("xsec_UP")->GetBinContent(bin);
    retVal("xsec_DOWN") = xSecHist("xsec_DOWN")->GetBinContent(bin);
  }
  return retVal;

  // Old code
  double a = 0.2979;
  double b = 17.626;
  double c = 67.632;
  double d = 3.463;
  retVal = (a / (1 + std::pow((stauM - b) / c, d)));
  retVal("xsec_UP") *= 1.03;
  retVal("xsec_DOWN") *= 0.97;
  return retVal;
}

bool StauAnalyser::electronMVAID(double mva, llvvLepton& lepton, IDType id)
{
  if(id == IDType::MediumID)
    id = IDType::TightID;

  double eta = lepton.electronInfoRef->sceta;
  bool pass = false;

  switch(id)
  {
  case IDType::LooseID:
    if(lepton.pt() < 20)
    {
      if(std::abs(eta) < 0.8)
      {
        if(mva > 0.925)
          pass = true;
      }
      else
      {
        if(std::abs(eta) < 1.479)
        {
          if(mva > 0.915)
            pass = true;
        }
        else
        {
          if(mva > 0.965)
            pass = true;
        }
      }
    }
    else
    {
      if(std::abs(eta) < 0.8)
      {
        if(mva > 0.905)
          pass = true;
      }
      else
      {
        if(std::abs(eta) < 1.479)
        {
          if(mva > 0.955)
            pass = true;
        }
        else
        {
          if(mva > 0.975)
            pass = true;
        }
      }
    }
    break;
  case IDType::TightID:
  default:
    if(lepton.pt() >= 20)
    {
      if(std::abs(eta) < 0.8)
      {
        if(mva > 0.925)
          pass = true;
      }
      else
      {
        if(std::abs(eta) < 1.479)
        {
          if(mva > 0.975)
            pass = true;
        }
        else
        {
          if(mva > 0.985)
            pass = true;
        }
      }
    }
    break;
  }

  return pass;
}

ValueWithSystematics<double> StauAnalyser::leptonIdAndIsoScaleFactor(ValueWithSystematics<llvvLepton>& selLepton)
{
  ValueWithSystematics<double> scaleFactor(1);

  ValueWithSystematics<double> idSF(1);
  ValueWithSystematics<double> isoSF(1);
  auto systematics = selLepton.Systematics();
  if(runSystematics && isMC)
  {
/*    scaleFactor.Systematic("elID_UP");
    scaleFactor.Systematic("elID_DOWN");
    scaleFactor.Systematic("muID_UP");
    scaleFactor.Systematic("muID_DOWN");
    scaleFactor.Systematic("elISO_UP");
    scaleFactor.Systematic("elISO_DOWN");
    scaleFactor.Systematic("muISO_UP");
    scaleFactor.Systematic("muISO_DOWN");
    
    for(auto& syst: systematics)
      scaleFactor.Systematic(syst);
    
    scaleFactor.Lock();// */
  }
  else
    systematics.clear();

  systematics.push_back("Value");
  
  for(auto& val: systematics)
  {
    auto& lepton = selLepton.GetSystematicOrValue(val);
  
    if(abs(lepton.id) == 11)
    {
      double pt = lepton.pt();
      double eta = lepton.electronInfoRef->sceta;
      if(std::abs(eta) < 1.479) // Electron in barrel
      {
        if(pt < 30)
        {
          if(val == "Value")
          {
            idSF.Value() = 0.8999;
            idSF("elID_UP")   += 0.0018;
            idSF("elID_DOWN") -= 0.0018;
            isoSF.Value() = 0.9417;
            isoSF("elISO_UP")   += 0.0019;
            isoSF("elISO_DOWN") -= 0.0019;
          }
          else
          {
            idSF(val) = 0.8999;
            isoSF(val) = 0.9417;
          }
        }
        else
        {
          if(val == "Value")
          {
            idSF.Value() = 0.9486;
            idSF("elID_UP")   += 0.0003;
            idSF("elID_DOWN") -= 0.0003;
            isoSF.Value() = 0.9804;
            isoSF("elISO_UP")   += 0.0003;
            isoSF("elISO_DOWN") -= 0.0003;
          }
          else
          {
            idSF(val) = 0.9486;
            isoSF(val) = 0.9804;
          }
        }
      }
      else // Electron in endcap
      {
        if(pt < 30)
        {
          if(val == "Value")
          {
            idSF.Value() = 0.7945;
            idSF("elID_UP")   += 0.0055;
            idSF("elID_DOWN") -= 0.0055;
            isoSF.Value() = 0.9471;
            isoSF("elISO_UP")   += 0.0037;
            isoSF("elISO_DOWN") -= 0.0037;
          }
          else
          {
            idSF(val) = 0.7945;
            isoSF(val) = 0.9471;
          }
        }
        else
        {
          if(val == "Value")
          {
            idSF.Value() = 0.8866;
            idSF("elID_UP")   += 0.0001;
            idSF("elID_DOWN") -= 0.0001;
            isoSF.Value() = 0.9900;
            isoSF("elISO_UP")   += 0.0002;
            isoSF("elISO_DOWN") -= 0.0002;
          }
          else
          {
            idSF(val) = 0.8866;
            isoSF(val) = 0.9900;
          }
        }
      }
    }
    else
    {
      double eta = lepton.eta();
      double pt  = lepton.pt();
      if(std::abs(eta) < 0.8) // Barrel muons
      {
        if(pt < 30)
        {
          if(val == "Value")
          {
            idSF.Value() = 0.9818;
            idSF("muID_UP")   += 0.0005;
            idSF("muID_DOWN") -= 0.0005;
            isoSF.Value() = 0.9494;
            isoSF("muISO_UP")   += 0.0015;
            isoSF("muISO_DOWN") -= 0.0015;
          }
          else
          {
            idSF(val) = 0.9818;
            isoSF(val) = 0.9494;
          }
        }
        else
        {
          if(val == "Value")
          {
            idSF.Value() = 0.9852;
            idSF("muID_UP")   += 0.0001;
            idSF("muID_DOWN") -= 0.0001;
            isoSF.Value() = 0.9883;
            isoSF("muISO_UP")   += 0.0003;
            isoSF("muISO_DOWN") -= 0.0003;
          }
          else
          {
            idSF(val) = 0.9852;
            isoSF(val) = 0.9883;
          }
        }
      }
      else
      {
        if(std::abs(eta) < 1.2) // Transition muons
        {
          if(pt < 30)
          {
            if(val == "Value")
            {
              idSF.Value() = 0.9829;
              idSF("muID_UP")   += 0.0009;
              idSF("muID_DOWN") -= 0.0009;
              isoSF.Value() = 0.9835;
              isoSF("muISO_UP")   += 0.0020;
              isoSF("muISO_DOWN") -= 0.0020;
            }
            else
            {
              idSF(val) = 0.9829;
              isoSF(val) = 0.9835;
            }
          }
          else
          {
            if(val == "Value")
            {
              idSF.Value() = 0.9852;
              idSF("muID_UP")   += 0.0002;
              idSF("muID_DOWN") -= 0.0002;
              isoSF.Value() = 0.9937;
              isoSF("muISO_UP")   += 0.0004;
              isoSF("muISO_DOWN") -= 0.0004;
            }
            else
            {
              idSF(val) = 0.9852;
              isoSF(val) = 0.9937;
            }
          }
        }
        else
        {
          if(pt < 30)
          {
            if(val == "Value")
            {
              idSF.Value() = 0.9869;
              idSF("muID_UP")   += 0.0007;
              idSF("muID_DOWN") -= 0.0007;
              isoSF.Value() = 0.9923;
              isoSF("muISO_UP")   += 0.0013;
              isoSF("muISO_DOWN") -= 0.0013;
            }
            else
            {
              idSF(val) = 0.9869;
              isoSF(val) = 0.9923;
            }
          }
          else
          {
            if(val == "Value")
            {
              idSF.Value() = 0.9884;
              idSF("muID_UP")   += 0.0001;
              idSF("muID_DOWN") -= 0.0001;
              isoSF.Value() = 0.9996;
              isoSF("muISO_UP")   += 0.0005;
              isoSF("muISO_DOWN") -= 0.0005;
            }
            else
            {
              idSF(val) = 0.9884;
              isoSF(val) = 0.9996;
            }
          }
        }
      }
    }
  }
  
  if(runSystematics && isMC)
  {
    scaleFactor = idSF * isoSF;
    scaleFactor("elID_UP");
    scaleFactor("elID_DOWN");
    scaleFactor("elISO_UP");
    scaleFactor("elISO_DOWN");
    scaleFactor("muID_UP");
    scaleFactor("muID_DOWN");
    scaleFactor("muISO_UP");
    scaleFactor("muISO_DOWN");// */
  }
  else
  {
    double tmp = idSF.Value() * isoSF.Value();
    scaleFactor.Reset();
    scaleFactor = tmp;
  }
  
  if(debugEvent)
  {
    analyserCout << "  ID SF:" << idSF.Value() << std::endl;
    for(auto& syst: idSF.Systematics())
      analyserCout << "   - " << syst << ": " << idSF(syst) << std::endl;
    analyserCout << "  ISO SF:" << isoSF.Value() << std::endl;
    for(auto& syst: isoSF.Systematics())
      analyserCout << "   - " << syst << ": " << isoSF(syst) << std::endl;
    analyserCout << "  ScaleFactor:" << scaleFactor.Value() << std::endl;
    for(auto& syst: scaleFactor.Systematics())
      analyserCout << "   - " << syst << ": " << scaleFactor(syst) << std::endl;
  }

  return scaleFactor;
}

ValueWithSystematics<double> StauAnalyser::tauScaleFactor(ValueWithSystematics<llvvTau>& selTau, TAU_E_ID eId)
{
  ValueWithSystematics<double> scaleFactor(1);
  auto systematics = selTau.Systematics();
  if(!(runSystematics&& isMC))
    systematics.clear();
  systematics.push_back("Value");

  // No correction necessary for tau ID
  // No correction necessary for normalization of Jet->Tau fake (if doing shape analysis, should investigate pt dependence of this)
  // No correction necessary for mu->Tau fake if using tight muon discriminator
  // Hadronic tau energy scale, no correction necessary
  // Tau charge misidentification rate, no correction necessary
  // This leaves only e->Tau fake, which must be corrected according to the anti-e discriminator used
  
  for(auto& val: systematics)
  {
    auto& tau = selTau.GetSystematicOrValue(val);
    bool isElectronFakingTau = false;
    bool isMuonFakingTau = false;
    bool isPromptTau = false;
  
    for(auto& genPart: gen)
    {
      if(genPart.status == 3 && abs(genPart.id) == 11)
      {
        if(deltaR(tau, genPart) < genMatchRCone)
        {
          isElectronFakingTau = true;
        }
      }
      if(genPart.status == 3 && abs(genPart.id) == 13)
      {
        if(deltaR(tau, genPart) < genMatchRCone)
        {
          isMuonFakingTau = true;
        }
      }
      if(genPart.status == 3 && abs(genPart.id) == 15)
      {
        if(deltaR(tau, genPart) < genMatchRCone)
        {
          isPromptTau = true;
        }
      }
    }
  
    if(isElectronFakingTau)
    {
      // Getting the correct endcap and barrel scale factors depending on the desired ID
      double barrelSF = 1, barrelShift = 0;
      double endcapSF = 1, endcapShift = 0;
      
      switch(eId)
      {
      case TAU_E_ID::antiELoose:
        barrelShift = barrelSF*0.05;
        endcapShift = endcapSF*0.1;
        break;
      case TAU_E_ID::antiEMedium:
        barrelSF = 0.95;
        endcapSF = 0.75;
        barrelShift = barrelSF*0.1;
        endcapShift = endcapSF*0.15;
        break;
      case TAU_E_ID::antiETight:
        barrelSF = 0.90;
        endcapSF = 0.70;
        barrelShift = barrelSF*0.15;
        endcapShift = endcapSF*0.2;
        break;
      case TAU_E_ID::antiEMva:
        barrelSF = 0.85;
        endcapSF = 0.65;
        barrelShift = barrelSF*0.2;
        endcapShift = endcapSF*0.25;
        break;
      case TAU_E_ID::antiEMva3Loose:
        barrelSF = 1.4; // +- 0.3
        endcapSF = 0.8; // +- 0.3
        barrelShift = 0.3;
        endcapShift = 0.3;
        break;
      case TAU_E_ID::antiEMva3Medium:
        barrelSF = 1.6; // +- 0.3
        endcapSF = 0.8; // +- 0.3
        barrelShift = 0.3;
        endcapShift = 0.3;
        break;
      case TAU_E_ID::antiEMva3Tight:
        barrelSF = 2.0; // +- 0.4
        endcapSF = 1.2; // +- 0.4
        barrelShift = 0.4;
        endcapShift = 0.4;
        break;
      case TAU_E_ID::antiEMva3VTight:
        barrelSF = 2.4; // +- 0.5
        endcapSF = 1.2; // +- 0.5
        barrelShift = 0.5;
        endcapShift = 0.5;
        break;
      case TAU_E_ID::antiEMva5Medium: // 1.6 +/- 0.3 for the barrel (std::abs(tauEta) < 1.5) and 1.1 +/- 0.3 for the endcap.
      default:
        barrelSF = 1.6;
        endcapSF = 1.1;
        barrelShift = 0.3;
        endcapShift = 0.3;
        break;
      }
      
      double SF = 1, SFshift = 0; 
      if(std::abs(tau.eta()) < 1.5)
      {
        SF = barrelSF;
        SFshift = barrelShift;
      }
      else
      {
        SF = endcapSF;
        SFshift = endcapShift;
      }
      
      if(val == "Value")
      {
        scaleFactor.Value() = SF;
        if(runSystematics && isMC)
        {
          scaleFactor("tauFromESF_UP") = SF + SFshift;
          scaleFactor("tauFromESF_DOWN") = SF - SFshift;
        }
      }
      else
      {
        scaleFactor(val) = SF;
      }
    }
    
    if(val == "Value" && runSystematics && isMC)
    {
      scaleFactor("tauID_UP") = scaleFactor.Value()*1.06;
      scaleFactor("tauID_DOWN") = scaleFactor.Value()*0.94;
      if(isMuonFakingTau)
      {
        scaleFactor("tauFromMu_UP") = scaleFactor.Value()*1.3;
        scaleFactor("tauFromMu_DOWN") = scaleFactor.Value()*0.7;
      }
      if(!isMuonFakingTau && !isElectronFakingTau && !isPromptTau)
      {
        scaleFactor("tauFromJet_UP") = scaleFactor.Value()*1.2;
        scaleFactor("tauFromJet_DOWN") = scaleFactor.Value()*0.8;
      }
    }
  }
  
  
  if(runSystematics && isMC)
  {
    scaleFactor.Systematic("tauID_UP");
    scaleFactor.Systematic("tauID_DOWN");
    scaleFactor.Systematic("tauFromESF_UP");
    scaleFactor.Systematic("tauFromESF_DOWN");
    scaleFactor.Systematic("tauFromMu_UP");
    scaleFactor.Systematic("tauFromMu_DOWN");
    scaleFactor.Systematic("tauFromJet_UP");
    scaleFactor.Systematic("tauFromJet_DOWN");
//    scaleFactor.Systematic("tauISOUP");
//    scaleFactor.Systematic("tauISODOWN");
    
    for(auto& syst: systematics)
      if(syst != "Value")
        scaleFactor.Systematic(syst);
    
    scaleFactor.Lock();
  }
  
  return scaleFactor;
}

ValueWithSystematics<double> StauAnalyser::computeMT2(const ValueWithSystematics<llvvTau>& tau, const ValueWithSystematics<llvvLepton>& lep, const ValueWithSystematics<TLorentzVector>& met)
{
  ValueWithSystematics<double> retVal;
  std::vector<std::string> tmpLoop;
  tmpLoop.push_back("Value");
  if(runSystematics && isMC)
  {
    loadSystematics(tmpLoop, tau);
    loadSystematics(tmpLoop, lep);
    loadSystematics(tmpLoop, met);
  }
  
  for(auto& val: tmpLoop)
  {
    if(val != "Value")
      retVal(val);
    
    auto& lep_ = lep.GetSystematicOrValue(val);
    auto& tau_ = tau.GetSystematicOrValue(val);
    auto& met_ = met.GetSystematicOrValue(val);

    double pa[3];
    double pb[3];
    double pmiss[3];
    double mn;
    
    pa[0] = lep_.M();
    pa[1] = lep_.px();
    pa[2] = lep_.py();
    pb[0] = tau_.M();
    pb[1] = tau_.px();
    pb[2] = tau_.py();
    pmiss[0] = 0;
    pmiss[1] = met_.Px();
    pmiss[2] = met_.Py();
    mn = 0;
    
    mt2_bisect::mt2 mt2_event;
    mt2_event.set_momenta(pa,pb,pmiss);
    mt2_event.set_mn(mn);

    retVal.GetSystematicOrValue(val) = mt2_event.get_mt2();
  }
  
  return retVal;
}



enum ID_Type {LooseID, MediumID, TightID};
enum TAU_E_ID {antiELoose, antiEMedium, antiETight, antiEMva, antiEMva3Loose, antiEMva3Medium, antiEMva3Tight, antiEMva3VTight, antiEMva5Medium};

double stauCrossSec(double stauM, double neutM);
bool electronMVAID(double mva, llvvLepton& lepton, ID_Type id);
double tauSF(llvvTau& tau, llvvGenParticleCollection& genPartColl, TAU_E_ID eId);
double leptonIdAndIsoScaleFactor(llvvLepton& lepton);
double leptonTauTriggerScaleFactor(llvvLepton& lepton, llvvTau& tau);
double efficiency(double m, double m0, double sigma, double alpha, double n, double norm);

/*****************************************************************************/
/* Return Codes:                                                             */
/*   0 - Everything OK                                                       */
/*   1 - Missing parameters_cfg.py configuration file                        */
/*   2 - Something failed                                                    */
/*****************************************************************************/
int main(int argc, char* argv[])
{
  if(argc < 2)
    std::cout << "Usage: " << argv[0] << " parameters_cfg.py" << std::endl, exit(1);

  size_t limit = 0;
  bool keepAllEvents = false;
  bool debugEvent = false;
  int skipEvents = 0;
  bool doOld = false;

  int fileIndex = 1;
  if(argc > 2)
  {
    std::stringstream parser;

    for(int i = 1; i < argc; ++i)
    {
      if(argv[i][0] != '-')
      {
        fileIndex = i;
        break;
      }

      std::string arg = argv[i];
      if(arg.find("--limit") != std::string::npos)
      {
        char first = argv[i+1][0];
        if(!isdigit(first))
          continue;

        parser << argv[i+1];
        parser >> limit;

        ++i;
        continue;
      }
      
      if(arg.find("--old") != std::string::npos)
      {
        doOld = true;
        continue;
      }
      
      if(arg.find("--keepAllEvents") != std::string::npos)
      {
        keepAllEvents = true;
        continue;
      }
      
      if(arg.find("--debugEvent") != std::string::npos)
      {
        debugEvent = true;
        continue;
      }

      if(arg.find("--skipEvents") != std::string::npos)
      {
        char first = argv[i+1][0];
        if(!isdigit(first))
          continue;

        parser << argv[i+1];
        parser >> skipEvents;

        ++i;
        continue;
      }
    }
  }

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();
  
  
  if(!doOld)
  {
    StauAnalyser testing(argv[fileIndex]);
    if(limit != 0)
      testing.SetEventLimit(limit);
    if(keepAllEvents)
      testing.KeepAllEvents();
    if(debugEvent)
      testing.SetDebugEvent(debugEvent);
    if(skipEvents != 0)
      testing.SetSkipEvents(skipEvents);
    testing.LoopOverEvents();
  
    return 0;// */
  }
  

  // Read parameters from the configuration file
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[fileIndex])->getParameter<edm::ParameterSet>("runProcess");

  bool isMC = runProcess.getParameter<bool>("isMC");
  double xsec = runProcess.getParameter<double>("xsec");
  std::vector<std::string> urls = runProcess.getParameter<std::vector<std::string> >("input");
  std::string baseDir = runProcess.getParameter<std::string>("dirName");
  std::string outdir = runProcess.getParameter<std::string>("outdir");
  std::string jecDir = runProcess.getParameter<std::string>("jecDir");
  bool runSystematics = runProcess.getParameter<bool>("runSystematics");
  bool saveSummaryTree = runProcess.getParameter<bool>("saveSummaryTree");
  bool exclusiveRun = runProcess.getParameter<bool>("exclusiveRun");
  double stauMtoPlot = 120;
  double neutralinoMtoPlot = 20;  // TStauStau mass point to plot by default
  if(runProcess.exists("stauMtoPlot"))
    stauMtoPlot = runProcess.getParameter<double>("stauMtoPlot");
  if(runProcess.exists("neutralinoMtoPlot"))
    neutralinoMtoPlot = runProcess.getParameter<double>("neutralinoMtoPlot");
  bool doSVfit = false;
  if(runProcess.exists("doSVfit"))
    doSVfit = runProcess.getParameter<bool>("doSVfit");
  bool applyScaleFactors = false;
  if(runProcess.exists("applyScaleFactors"))
    applyScaleFactors = runProcess.getParameter<bool>("applyScaleFactors");
  bool debug = false;
  if(runProcess.exists("debug"))
    debug = runProcess.getParameter<bool>("debug");
  bool doTightTauID = false;
  if(runProcess.exists("doTightTauID"))
    doTightTauID = runProcess.getParameter<bool>("doTightTauID");
  bool doDDBkg = false;
  if(runProcess.exists("doDDBkg"))
    doDDBkg = runProcess.getParameter<bool>("doDDBkg");
  bool outputEventList = false;
  if(runProcess.exists("outputEventList"))
    outputEventList   = runProcess.getParameter<bool>("outputEventList");

  if(debug)
    std::cout << "Finished loading config file" << std::endl;

  // Hardcoded Values
  double sqrtS          =  8;      // Center of mass energy
  double minElPt        = 24;      // Selected electron pT and eta
  double maxElEta       =  2.1;
  double ECALGap_MinEta =  1.4442; // ECAL gap parameters
  double ECALGap_MaxEta =  1.5660;
  double minMuPt        = 20;      // Selected muon pT and eta
  double maxMuEta       =  2.1;
  double minTauPt       = 20;      // Selected tau pT and eta (I was using 25)
  double maxTauEta      =  2.3;
  double maxJetEta      =  4.7;    // Selected jet eta

  char runPeriod = 'X';

  // Setting up -------------------------------------------------------------------------
  if(debug)
    std::cout << "Setting up" << std::endl;
  gSystem->Exec(("mkdir -p " + outdir).c_str());
  std::string url = urls[0];
  std::string outFileUrl(gSystem->BaseName(url.c_str()));
  while(outFileUrl.find(".root", 0) != std::string::npos)
    outFileUrl.replace(outFileUrl.find(".root", 0), 5, "");
  if(doDDBkg)
  {
    for(auto & url : urls)
    {
      if(url.find("DD", 0) != std::string::npos)
        url.replace(url.find("DD", 0), 2, "");
    }
    for(auto & url : urls)
    {
      std::cout << "New url: " << url << std::endl;
    }
  }
  std::string outUrl = outdir;
  outUrl += "/";
  outUrl += outFileUrl + ".root";

  TString turl(url);
  bool isV0JetsMC(isMC && (turl.Contains("DYJetsToLL_50toInf") || turl.Contains("WJets")));
  bool isStauStau(isMC && turl.Contains("TStauStau"));
  if(!isMC)
  {
    if(turl.Contains("2012A"))
      runPeriod = 'A';
    if(turl.Contains("2012B"))
      runPeriod = 'B';
    if(turl.Contains("2012C"))
      runPeriod = 'C';
    if(turl.Contains("2012D"))
      runPeriod = 'D';
  }

//  TH2D *etauFR = NULL;
//  TH2D *mutauFR = NULL;
//  TH2D *etauPR = NULL;
//  TH2D *mutauPR = NULL;
  TH1 *etauFR = NULL;
  TH1 *mutauFR = NULL;
  TH1 *etauPR = NULL;
  TH1 *mutauPR = NULL;
  if(doDDBkg)
  {
    TDirectory* cwd = gDirectory;

    std::string RatesFileName = gSystem->ExpandPathName("$CMSSW_BASE/src/UserCode/llvv_fwk/data/TStauStau/rates.root");
    std::cout << "Trying to open Rates file: " << RatesFileName << std::endl;
    TFile RatesFile(RatesFileName.c_str(), "READ");
    cwd->cd();

    etauFR  = static_cast<TH1*>(RatesFile.Get("data-Zprompt/data-Zprompt_InvMET_OS_etaSelectedTau_FR")->Clone("etauFR"));
    mutauFR = static_cast<TH1*>(etauFR->Clone("mutauFR"));

    //PRompt rate should probably be pt dependent
    etauPR  = static_cast<TH1*>(RatesFile.Get("Z #rightarrow ll/Zrightarrowll_InvMET_OS_etaSelectedTau")->Clone("etauPR"));
    mutauPR = static_cast<TH1*>(etauPR->Clone("mutauPR"));

/*    if(isMC)
    {
//      etauFR  = static_cast<TH2D*>(FRFile.Get("W + Jets/WJets_leadingE_varptetaSelectedTau_FR")->Clone("etauFR"));
//      mutauFR = static_cast<TH2D*>(FRFile.Get("W + Jets/WJets_leadingMu_varptetaSelectedTau_FR")->Clone("mutauFR"));
      etauFR  = static_cast<TH1*>(FRFile.Get("W + Jets/WJets_leadingE_varetaSelectedTau_FR")->Clone("etauFR"));
      mutauFR = static_cast<TH1*>(FRFile.Get("W + Jets/WJets_leadingMu_varetaSelectedTau_FR")->Clone("mutauFR"));
    }
    else
    {
//      etauFR  = static_cast<TH2D*>(FRFile.Get("data/data_leadingE_varptetaSelectedTau_FR")->Clone("etauFR"));
//      mutauFR = static_cast<TH2D*>(FRFile.Get("data/data_leadingMu_varptetaSelectedTau_FR")->Clone("mutauFR"));
      etauFR  = static_cast<TH1*>(FRFile.Get("data/data_leadingE_varetaSelectedTau_FR")->Clone("etauFR"));
      mutauFR = static_cast<TH1*>(FRFile.Get("data/data_leadingMu_varetaSelectedTau_FR")->Clone("mutauFR"));
    }
//    etauPR  = static_cast<TH2D*>(PRFile.Get("Z #rightarrow ll/Zrightarrowll_leadingE_varptetaSelectedTau_FR")->Clone("etauPR"));
//    mutauPR = static_cast<TH2D*>(PRFile.Get("Z #rightarrow ll/Zrightarrowll_leadingMu_varptetaSelectedTau_FR")->Clone("mutauPR"));
    etauPR  = static_cast<TH1*>(PRFile.Get("Z #rightarrow ll/Zrightarrowll_leadingE_varetaSelectedTau_FR")->Clone("etauPR"));
    mutauPR = static_cast<TH1*>(PRFile.Get("Z #rightarrow ll/Zrightarrowll_leadingMu_varetaSelectedTau_FR")->Clone("mutauPR")); // */

    if(etauFR == NULL || mutauFR == NULL)
    {
      std::cout << "Unable to open fake rate histograms. Stopping execution." << std::endl;
      return 2;
    }
    if(etauPR == NULL || mutauPR == NULL)
    {
      std::cout << "Unabe to open prompt rate histograms. Stopping execution." << std::endl;
      return 2;
    }
  }


  TTree* summaryTree = NULL;
  TFile* summaryOutFile = NULL;
  if(saveSummaryTree)
  {
    TDirectory* cwd = gDirectory;

    std::string summaryOutUrl = outUrl;
    summaryOutUrl.replace(summaryOutUrl.find(".root", 0), 5, "_summary.root");
    std::cout << "Saving summary results in " << summaryOutUrl << std::endl;
    summaryOutFile = new TFile(summaryOutUrl.c_str(), "RECREATE");

    summaryTree = new TTree("Events", "Events");

    summaryTree->SetDirectory(summaryOutFile);

    cwd->cd();
  }
  
  ofstream eventListFile;
  if(outputEventList)
  {
    std::string eventlistOutUrl = outUrl;
    eventlistOutUrl.replace(eventlistOutUrl.find(".root", 0), 5, "_eventlist.txt");
    std::cout << "Saving event list to " << eventlistOutUrl << std::endl;
    eventListFile.open(eventlistOutUrl);
  }




  /***************************************************************************/
  /*                         Initializing Histograms                         */
  /***************************************************************************/
  if(debug)
    std::cout << "Initializing histograms" << std::endl;
  SmartSelectionMonitor mon;
  TH1D *eventflow = (TH1D*)mon.addHistogram(new TH1D("eventflow", ";;Events", 8, 0, 8));
  eventflow->GetXaxis()->SetBinLabel(1, "HLT");
  eventflow->GetXaxis()->SetBinLabel(2, "MET > 30");
  eventflow->GetXaxis()->SetBinLabel(3, "> 1l");
  eventflow->GetXaxis()->SetBinLabel(4, "> 1#tau");
  eventflow->GetXaxis()->SetBinLabel(5, "B-veto");
  eventflow->GetXaxis()->SetBinLabel(6, "OS");
  eventflow->GetXaxis()->SetBinLabel(7, "lep veto");
  eventflow->GetXaxis()->SetBinLabel(8, "SVfit");

  mon.addHistogram(new TH1D("nup", ";NUP;Events", 10, 0, 10));

  // Pile Up
//  mon.addHistogram(new TH1D("nvtxAll", ";Vertices;Events", 50, -0.5, 49.5));
  mon.addHistogram(new TH1D("nvtx", ";Vertices;Events", 50, -0.5, 49.5));
  mon.addHistogram(new TH1D("nvtxraw", ";Vertices;Events", 50, -0.5, 49.5));
  mon.addHistogram(new TH1D("rho", ";#rho;Events", 25, 0, 25));
  mon.addHistogram(new TH1D("rho25", ";#rho(#eta<2.5);Events", 25, 0, 25));


  // Leptons
  mon.addHistogram(new TH1D("nlep", ";nlep;Events", 10, 0, 10));
  mon.addHistogram(new TH1D("ptSelectedLep", ";p_{T}^{l};Events", 50, 0, 100));
  mon.addHistogram(new TH1D("etaSelectedLep", ";#eta^{l};Events", 25, -2.6, 2.6));
  mon.addHistogram(new TH1D("chargeSelectedLep", ";q^{l};Events", 5, -2, 2));
  TH1D *leptonCutFlow = (TH1D*)mon.addHistogram(new TH1D("leptonCutFlow", ";;Leptons", 4, 0, 4));
  leptonCutFlow->GetXaxis()->SetBinLabel(1, "All");
  leptonCutFlow->GetXaxis()->SetBinLabel(2, "ID");
  leptonCutFlow->GetXaxis()->SetBinLabel(3, "Kin");
  leptonCutFlow->GetXaxis()->SetBinLabel(4, "Iso");

  // Lepton Isolation
  mon.addHistogram(new TH1D("isomu", "RelIso(#mu);;Leptons", 100, -0.5, 9.5));
  mon.addHistogram(new TH1D("isoele", "RelIso(ele);;Leptons", 100, -0.5, 9.5));

  // Taus
  mon.addHistogram(new TH1D("ntaus", ";ntaus;Events", 10, 0, 10));
  mon.addHistogram(new TH1D("ptSelectedTau", ";p_{T}^{#tau};Events", 50, 0, 100));
  mon.addHistogram(new TH1D("ptSelectedTauExtended", ";p_{T}^{#tau};Events", 50, 0, 250));
  mon.addHistogram(new TH1D("etaSelectedTau", ";#eta^{#tau};Events", 25, -2.6, 2.6));
  mon.addHistogram(new TH1D("chargeSelectedTau", ";q^{#tau};Events", 5, -2, 2));
  mon.addHistogram(new TH1D("dzSelectedTau", ";dz^{#tau};Events", 25, 0, 2));
  mon.addHistogram(new TH1D("emfracSelectedTau", ";emf^{#tau};Events", 25, 0, 5));
  TH1D *tauCutFlow = (TH1D*)mon.addHistogram(new TH1D("tauCutFlow", ";;#tau", 6, 0, 6));
  tauCutFlow->GetXaxis()->SetBinLabel(1, "All");
  tauCutFlow->GetXaxis()->SetBinLabel(2, "PF");
  tauCutFlow->GetXaxis()->SetBinLabel(3, "ID");
  tauCutFlow->GetXaxis()->SetBinLabel(4, "Quality");
  tauCutFlow->GetXaxis()->SetBinLabel(5, "Kin");
  tauCutFlow->GetXaxis()->SetBinLabel(6, "Iso");
  TH1D *tauID = (TH1D*)mon.addHistogram(new TH1D("tauID", ";;#tau", 5, 0, 5));
  tauID->GetXaxis()->SetBinLabel(1, "All");
  tauID->GetXaxis()->SetBinLabel(2, "Medium comb iso");
  tauID->GetXaxis()->SetBinLabel(3, "Decay mode");
  tauID->GetXaxis()->SetBinLabel(4, "Not e");
  tauID->GetXaxis()->SetBinLabel(5, "Not #mu");

  // Jets
  mon.addHistogram(new TH1D("njets", ";njets;Events", 6, 0, 6));
  mon.addHistogram(new TH1D("nbjets", ";njets;Events", 6, 0, 6));
  mon.addHistogram(new TH1D("jetleadpt", ";p_{T}^{jet};Events", 25, 0, 500));
  mon.addHistogram(new TH1D("jetleadeta", ";#eta^{jet};Events", 50, -5, 5));
  mon.addHistogram(new TH1D("jetcsv", ";csv;jets", 25, 0, 1));
  TH1D *jetCutFlow = (TH1D*)mon.addHistogram(new TH1D("jetCutFlow", ";;jets", 6, 0, 6));
  jetCutFlow->GetXaxis()->SetBinLabel(1, "All");
  jetCutFlow->GetXaxis()->SetBinLabel(2, "PF Loose");
  jetCutFlow->GetXaxis()->SetBinLabel(3, "ID");
  jetCutFlow->GetXaxis()->SetBinLabel(4, "Kin");
  jetCutFlow->GetXaxis()->SetBinLabel(5, "Iso");
  jetCutFlow->GetXaxis()->SetBinLabel(6, "B-jet");

  // MET
  mon.addHistogram(new TH1D("MET", ";MET [GeV];Events", 20, 0, 200));

  // MT
  mon.addHistogram(new TH1D("MT", ";MT [GeV];Events", 25, 0, 200));
  mon.addHistogram(new TH1D("MTTau", ";MT(#tau) [GeV];Events", 25, 0, 200));
  mon.addHistogram(new TH1D("SumMT", ";SumMT [GeV];Events", 25, 0, 200));

  // Deconstructed MT: https://indico.cern.ch/event/344807/
  mon.addHistogram(new TH1D("Q80", ";Q_{80};Events", 30, -2, 1));
  mon.addHistogram(new TH1D("Q100", ";Q_{100};Events", 30, -2, 1));
  mon.addHistogram(new TH1D("cosPhi", ";cos#Phi;Events", 30, -1, 1));
  mon.addHistogram(new TH1D("Q80Tau", ";Q_{80};Events", 30, -2, 1));
  mon.addHistogram(new TH1D("Q100Tau", ";Q_{100};Events", 30, -2, 1));
  mon.addHistogram(new TH1D("cosPhiTau", ";cos#Phi;Events", 30, -1, 1));

  // MT2
  mon.addHistogram(new TH1D("MT2", ";M_{T2} [GeV];Events", 25, 0, 500));

  // SVFit Mass
  if(doSVfit)
    mon.addHistogram(new TH1D("SVFitMass", ";M_{SVFit};Events", 50, 0, 500));

  // Invariant Mass
  mon.addHistogram(new TH1D("InvMass", ";M_{l-#tau};Events", 50, 0, 500));

  // Angles
  mon.addHistogram(new TH1D("deltaAlphaLepTau", ";#Delta#alpha_{l-#tau}(Lab);Events", 30, 0, TMath::Pi()));
  mon.addHistogram(new TH1D("deltaRLepTau", ";#Delta R_{l-#tau}(Lab);Events", 40, 0, 8));
  mon.addHistogram(new TH1D("deltaPhiLepTauMET", ";#Delta#phi_{l#tau-MET}(Lab);Events", 30, -TMath::Pi(), TMath::Pi()));
  mon.addHistogram(new TH1D("deltaPhiLepTau", ";#Delta#phi_{l-#tau}(Lab);Events", 30, -TMath::Pi(), TMath::Pi()));
  mon.addHistogram(new TH1D("cosThetaTau", ";cos#theta_{#tau}(Lab);Events", 30, -1, 1));
  mon.addHistogram(new TH1D("cosThetaLep", ";cos#theta_{l}(Lab);Events", 30, -1, 1));
  mon.addHistogram(new TH1D("deltaPhiLepMETCS", ";#Delta#phi_{l-MET}(CS);Events", 30, -TMath::Pi(), TMath::Pi()));
  mon.addHistogram(new TH1D("cosThetaCS", ";cos#theta(CS);Events", 30, -1, 1));
  mon.addHistogram(new TH1D("minDeltaPhiMETJetPt40", ";min(#Delta#phi_{MET-Jet40});Events", 20, -TMath::Pi(), TMath::Pi()));

  // 2D variables
  mon.addHistogram(new TH2D("metVsPtl", ";p_{T}(l);MET", 50, 0, 100, 25, 0, 200));
  mon.addHistogram(new TH2D("metVsPtTau", ";p_{T}(#tau);MET", 50, 0, 100, 25, 0, 200));
  mon.addHistogram(new TH2D("metPtVsmetEt", ";met.Et();met.pt()", 25, 0, 200, 25, 0, 200));
  //  Deconstructed MT 2D Plots:
  mon.addHistogram(new TH2D("Q80VsCosPhi", ";cos#Phi;Q_{80}", 20, -1, 1, 20, -2, 1));
  mon.addHistogram(new TH2D("Q100VsCosPhi", ";cos#Phi;Q_{100}", 20, -1, 1, 20, -2, 1));
  mon.addHistogram(new TH2D("Q80VsCosPhiTau", ";cos#Phi;Q_{80}", 20, -1, 1, 20, -2, 1));
  mon.addHistogram(new TH2D("Q100VsCosPhiTau", ";cos#Phi;Q_{100}", 20, -1, 1, 20, -2, 1));


  if(outputEventList)
  {
    eventListFile << std::setw(EVENTLISTWIDTH) << "Run #" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "Lumi #" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "Event #" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "selected" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "weight" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "PUweight" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "xsecweight" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "triggerSF" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "lepSF" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "tauSF" << "|";
    if(doDDBkg)
      eventListFile << std::setw(EVENTLISTWIDTH) << "DDweight" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "MET" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "num" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "nLep" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "nTau" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "nJet" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "nBJet" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "MTLep" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "MTTau" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "SumMT" << "|";
    eventListFile << std::setw(EVENTLISTWIDTH) << "MT2" << "|";
    eventListFile << std::endl;
  }



  /***************************************************************************/
  /*                          Prepare for Event Loop                         */
  /***************************************************************************/
  if(debug)
    std::cout << "Preparing for event loop" << std::endl;
  fwlite::ChainEvent ev(urls);
  const size_t totalEntries = ev.size();

  // MC normalization to 1/pb
  double nInitEvent = 1.;
  if(isMC)
    nInitEvent = (double) utils::getMergeableCounterValue(urls, "startCounter");
  double xsecWeight = xsec/nInitEvent;
  if(!isMC)
    xsecWeight = 1.;

  // Jet Energy Scale and Uncertainties
  jecDir = gSystem->ExpandPathName(jecDir.c_str());
  FactorizedJetCorrector *jesCor = utils::cmssw::getJetCorrector(jecDir, isMC);
  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt"));

  // Muon Energy Scale and Uncertainties
  MuScleFitCorrector* muCor = getMuonCorrector(jecDir, url);

  // Pileup Weighting: Based on vtx (For now?)
  edm::LumiReWeighting* LumiWeights = NULL;
  utils::cmssw::PuShifter_t PuShifters;
  double PUNorm[] = {1, 1, 1};
  if(isMC)
  {
    std::vector<double> dataPileupDistributionDouble = runProcess.getParameter<std::vector<double> >("datapileup");
    std::vector<float> dataPileupDistribution;
    for(unsigned int i = 0; i < dataPileupDistributionDouble.size(); ++i)
      dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);
    std::vector<float> mcPileupDistribution;

    utils::getMCPileupDistribution(ev, dataPileupDistribution.size(), mcPileupDistribution);
    while(mcPileupDistribution.size() < dataPileupDistribution.size()) mcPileupDistribution.push_back(0.0);
    while(mcPileupDistribution.size() > dataPileupDistribution.size()) dataPileupDistribution.push_back(0.0);

    LumiWeights= new edm::LumiReWeighting(mcPileupDistribution, dataPileupDistribution);
    PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution, 0.05);
    utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
  }

  gROOT->cd(); //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

  DuplicatesChecker duplicatesChecker;
//  int nDuplicates(0);
  int step = int(totalEntries/50);

  // Redirect stdout and stderr to a temporary buffer, then output buffer after event loop
  std::ostream myCout(std::cout.rdbuf());
  std::streambuf *coutbuf = std::cout.rdbuf();
  std::streambuf *cerrbuf = std::cerr.rdbuf();

//  std::stringstream buffer;
//  std::cout.rdbuf(buffer.rdbuf());
//  std::cerr.rdbuf(buffer.rdbuf());
  std::ofstream devnull("/dev/null");
  std::cout.rdbuf(devnull.rdbuf());
  std::cerr.rdbuf(devnull.rdbuf());

  // Variables used in loop
  if(debug)
    myCout << "  Declaring all variables used in loop" << std::endl;
  int nvtx = 0;
  bool selected = false;
  bool isetau   = false;
  bool ismutau  = false;
  bool istautau = false;
  bool isloose  = false;
  bool istight  = false;
  std::vector<TString> chTags;
  std::vector<bool> triggerBits;
  bool triggeredOn = false;
  llvvMet met;
  double rho = 0;
  double rho25 = 0;
  double crossSection = xsec;
  double weight = 1.;
  double weight_plus = 1.;
  double weight_minus = 1.;
  double puWeight = 1.;
  double triggerSF = 1.;
  double leptonIdIsoSF = 1.;
  double tauSF = 1.;
  llvvLeptonCollection selLeptons;
  llvvJetExtCollection selJets;
  llvvJetCollection selJetsOut;
  llvvJetCollection selBJets;
  llvvTauCollection selTaus;
  int nJets = 0;
  int nBJets = 0;
  int tauIndex = -1, leptonIndex = -1;
  bool isOS = false;
  bool isMultilepton = false;
  bool isSVfit = true;
  double mass = -1;
  double invMass = -1;
  double mt = -1;
  double mtTau = -1;
  double sumMt = -1;
  double Q80 = 2;
  double Q100 = 2;
  double cosPhi = -10;
  double Q80Tau = 2;
  double Q100Tau = 2;
  double cosPhiTau = -10;
  double mt2 = -1;
  double stauMass = 0;
  double neutralinoMass = 0;
  double deltaAlphaLepTau = 0;
  double deltaRLepTau = 0;
  double deltaPhiLepTauMET = 0;
  double deltaPhiLepTau = 0;
  double cosThetaTau = 0;
  double cosThetaLep = 0;
  double deltaPhiLepMETCS = 0;
  double cosThetaCS = 0;
  double minDeltaPhiMETJetPt40 = 0;
  double tauLeadPt = 0;
  double lepLeadPt = 0;
  double maxPtSum = 0;
  bool isPromptLep = false;
  bool isPromptTau = false;
//  int nTauJets = 0;

  // Prepare summary tree
  if(saveSummaryTree)
  {
    if(debug)
      std::cout << "  Defining all branches in output root file" << std::endl;

    TDirectory* cwd = gDirectory;
    summaryOutFile->cd();

    // Dataset specific variables
    summaryTree->Branch("isMC", &isMC);
    summaryTree->Branch("xSecWeight", &xsecWeight);

    // Event specific variables
    summaryTree->Branch("selected", &selected);
    summaryTree->Branch("isetau",   &isetau);
    summaryTree->Branch("ismutau",  &ismutau);
    summaryTree->Branch("istautau", &istautau);
    summaryTree->Branch("isloose",  &isloose);
    summaryTree->Branch("istight",  &istight);
    summaryTree->Branch("chTags", &chTags);
    summaryTree->Branch("nvtx", &nvtx);
    summaryTree->Branch("triggerBits", &triggerBits);
    summaryTree->Branch("triggeredOn", &triggeredOn);
    summaryTree->Branch("rho", &rho);
    summaryTree->Branch("rho25", &rho25);
    summaryTree->Branch("met", &met);
    summaryTree->Branch("crossSection", &crossSection);
    summaryTree->Branch("weight", &weight);
    summaryTree->Branch("weight_plus", &weight_plus);
    summaryTree->Branch("weight_minus", &weight_minus);
    summaryTree->Branch("puWeight", &puWeight);
    summaryTree->Branch("triggerSF", &triggerSF);
    summaryTree->Branch("leptonIdIsoSF", &leptonIdIsoSF);
    summaryTree->Branch("tauSF", &tauSF);
    summaryTree->Branch("selLeptons", &selLeptons);
    summaryTree->Branch("selTaus", &selTaus);
//    summaryTree->Branch("selJets", &selJetsOut);
//    summaryTree->Branch("selBJets", &selBJets);
    summaryTree->Branch("nJets", &nJets);
    summaryTree->Branch("nBJets", &nBJets);
    summaryTree->Branch("isOS", &isOS);
    summaryTree->Branch("isMultilepton", &isMultilepton);
    summaryTree->Branch("isSVfit", &isSVfit);
    summaryTree->Branch("tauIndex", &tauIndex);
    summaryTree->Branch("leptonIndex", &leptonIndex);
    if(doSVfit)
    {
      summaryTree->Branch("SVFitMass", &mass);
    }
    summaryTree->Branch("InvariantMass", &invMass);
    summaryTree->Branch("MT", &mt);
    summaryTree->Branch("MTTau", &mtTau);
    summaryTree->Branch("SumMT", &sumMt);
    summaryTree->Branch("Q80", &Q80);
    summaryTree->Branch("Q100", &Q100);
    summaryTree->Branch("cosPhi", &cosPhi);
    summaryTree->Branch("Q80Tau", &Q80Tau);
    summaryTree->Branch("Q100Tau", &Q100Tau);
    summaryTree->Branch("cosPhiTau", &cosPhiTau);
    summaryTree->Branch("MT2", &mt2);
    summaryTree->Branch("stauMass", &stauMass);
    summaryTree->Branch("neutralinoMass", &neutralinoMass);
    summaryTree->Branch("deltaAlphaLepTau", &deltaAlphaLepTau);
    summaryTree->Branch("deltaRLepTau", &deltaRLepTau);
    summaryTree->Branch("deltaPhiLepTauMET", &deltaPhiLepTauMET);
    summaryTree->Branch("deltaPhiLepTau", &deltaPhiLepTau);
    summaryTree->Branch("cosThetaTau", &cosThetaTau);
    summaryTree->Branch("cosThetaLep", &cosThetaLep);
    summaryTree->Branch("deltaPhiLepMETCS", &deltaPhiLepMETCS);
    summaryTree->Branch("cosThetaCS", &cosThetaCS);
    summaryTree->Branch("minDeltaPhiMETJetPt40", &minDeltaPhiMETJetPt40);
    summaryTree->Branch("tauLeadPt", &tauLeadPt);
    summaryTree->Branch("lepLeadPt", &lepLeadPt);
    summaryTree->Branch("isPromptLep", &isPromptLep);
    summaryTree->Branch("isPromptTau", &isPromptTau);

    cwd->cd();
  }



  /***************************************************************************/
  /*                                Event Loop                               */
  /***************************************************************************/
  myCout << "       Progress Bar:0%      20%       40%       60%       80%      100%" << std::endl;
  myCout << "Scanning the ntuple:";

  // Loop on events
  for(size_t iev = 0; iev < totalEntries; ++iev)
  {
    #if defined(DEBUG_EVENT)
    if(iev < skipEvents)
      continue;
    if(debugEvent)
      myCout << "## Event " << iev << std::endl;
    #endif
    if(iev%step == 0)
      myCout << "_" << std::flush;

    // Init variables
    selected = false;
    isetau   = false;
    ismutau  = false;
    istautau = false;
    isloose  = false;
    istight  = false;
    deltaAlphaLepTau = 0;
    deltaRLepTau = 0;
    deltaPhiLepTauMET = 0;
    deltaPhiLepTau = 0;
    cosThetaTau = 0;
    cosThetaLep = 0;
    deltaPhiLepMETCS = 0;
    cosThetaCS = 0;
    minDeltaPhiMETJetPt40 = 0;
    nvtx = 0;
    weight = 1.;
    weight_plus = 1.;
    weight_minus = 1.;
    puWeight = 1.;
    triggerSF = 1.;
    leptonIdIsoSF = 1.;
    tauSF = 1.;
    chTags.clear();
    selLeptons.clear();
    nJets = 0;
    nBJets = 0;
    selJets.clear();
    selJetsOut.clear();
    selBJets.clear();
    selTaus.clear();
    tauIndex = -1, leptonIndex = -1;
    isOS = false;
    isMultilepton = false;
    isSVfit = false;
    mass = -1;
    invMass = -1;
    mt = -1;
    mtTau = -1;
    sumMt = -1;
    Q80 = 2;
    Q100 = 2;
    cosPhi = -10;
    Q80Tau = 2;
    Q100Tau = 2;
    cosPhiTau = -10;
    mt2 = -1;
    stauMass = -1;
    neutralinoMass = -1;
    tauLeadPt = 0;
    lepLeadPt = 0;
    maxPtSum = 0;
    isPromptLep = false;
    isPromptTau = false;
//    nTauJets = 0;

    // Prepare tags to fill the histograms
    chTags.push_back("all");

    // Load the event content from tree
    ev.to(int(iev));


    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Loading collections" << std::endl;
    #endif

    /**** Get information/collections from the event ****/
    // Number of vertexes
    fwlite::Handle<int> nvtxHandle;
    nvtxHandle.getByLabel(ev, "llvvObjectProducersUsed", "nvtx");
    if(nvtxHandle.isValid()) nvtx = *nvtxHandle;

    // Collection of generated particles
    fwlite::Handle<llvvGenEvent> genEventHandle;
    genEventHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!genEventHandle.isValid())
    {
      std::cout << "llvvGenEvent Object NotFound" << std::endl;
      continue;
    }
    llvvGenEvent genEv = *genEventHandle;
    /**** Remove double counting if running on exclusive samples ****/
    if(exclusiveRun && isV0JetsMC)
    {
      if(genEv.nup > 5) // Drop V+{1,2,3,4}Jets from VJets samples to avoid double counting (but keep V+0Jets) [V = W,Z]
        continue;
    }

    /**** Get LHE comments with mass info ****/
    if(isStauStau)
    {
      fwlite::Handle<LHEEventProduct> LHEHandle;
      LHEHandle.getByLabel(ev, "source");
      if(!LHEHandle.isValid())
      {
        std::cout << "LHEEventProduct Object not Found" << std::endl;
        continue;
      }
      if(LHEHandle->comments_size() == 0)
        continue;

      for(auto comment = LHEHandle->comments_begin(); comment != LHEHandle->comments_end(); ++comment)
      {
        auto modelPos = comment->find("# model TStauStau_");
        if(modelPos != std::string::npos)
        {
          std::stringstream tmp;
          auto numPos = comment->find_first_of("1234567890", modelPos);

          tmp << comment->substr(numPos, comment->find("_", numPos)-numPos);
          tmp >> stauMass;
          tmp.clear();

          numPos = comment->find("_", numPos);
          numPos = comment->find_first_of("1234567890", numPos);
          tmp << comment->substr(numPos, comment->find("\n", numPos)-numPos);
          tmp >> neutralinoMass;

          break;
        }
      }
    }

    // Trigger Bits
    fwlite::Handle<std::vector<bool> > triggerBitsHandle;
    triggerBitsHandle.getByLabel(ev, "llvvObjectProducersUsed", "triggerBits");
    if(!triggerBitsHandle.isValid())
    {
      std::cout << "triggerBits Object NotFound" << std::endl;
      continue;
    }
    triggerBits = *triggerBitsHandle;

    //bool singleETrigger = triggerBits[13]; // HLT_Ele27_WP80_v*
    //bool singleMuTrigger = triggerBits[15]; // HLT_IsoMu24_v*
    bool TauPlusE2012A = triggerBits[18]; // HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*
    bool TauPlusMu2012A = triggerBits[22]; // HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v*
    //bool TauPlusE2012A = triggerBits[19]; // HLT_Ele22_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v
    //bool TauPlusMu2012A = triggerBits[23]; // HLT_IsoMu20_eta2p1_LooseIsoPFTau20_v*
    bool TauPlusE2012B = triggerBits[17]; // HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v*
    bool TauPlusMu2012B = triggerBits[21]; // HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*
    //bool TauPlusE2012D = triggerBits[16]; // HLT_Ele13_eta2p1_WP90Rho_LooseIsoPFTau20_v*
    //bool TauPlusMu2012D = triggerBits[20]; // HLT_IsoMu8_eta2p1_LooseIsoPFTau20_v*
    bool TauPlusETrigger = TauPlusE2012A || TauPlusE2012B;
    bool TauPlusMuTrigger = TauPlusMu2012A || TauPlusMu2012B;

    triggeredOn = TauPlusETrigger || TauPlusMuTrigger;
//    if(TauPlusETrigger)
//      chTags.push_back("TauPlusE");
//    if(TauPlusMuTrigger)
//      chTags.push_back("TauPlusMu");

    // Rest of Gen Particles
    fwlite::Handle<llvvGenParticleCollection> genPartCollHandle;
    genPartCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!genPartCollHandle.isValid())
    {
      std::cout << "llvvGenParticleCollection Object NotFound" << std::endl;
      continue;
    }
    llvvGenParticleCollection gen = *genPartCollHandle;

    // Collection of leptons
    fwlite::Handle<llvvLeptonCollection> leptonCollHandle;
    leptonCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!leptonCollHandle.isValid())
    {
      std::cout << "llvvLeptonCollection Object NotFound" << std::endl;
      continue;
    }
    llvvLeptonCollection leptons = *leptonCollHandle;

    // Electron Information Collection
    fwlite::Handle<llvvElectronInfoCollection> electronInfoCollHandle;
    electronInfoCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!electronInfoCollHandle.isValid())
    {
      std::cout << "llvvElectronInfoCollection Object NotFound" << std::endl;
      continue;
    }

    // Muon Information Collection
    fwlite::Handle<llvvMuonInfoCollection> muonInfoCollHandle;
    muonInfoCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!muonInfoCollHandle.isValid())
    {
      std::cout << "llvvMuonInfoCollection Object NotFound" << std::endl;
      continue;
    }

    // Tau Collection
    fwlite::Handle<llvvTauCollection> tauCollHandle;
    tauCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!tauCollHandle.isValid())
    {
      std::cout << "llvvTauCollection Object NotFound" << std::endl;
      continue;
    }
    llvvTauCollection taus = *tauCollHandle;

    // Boosted tau Collection
    fwlite::Handle<llvvTauCollection> boostedTauCollHandle;
    boostedTauCollHandle.getByLabel(ev, "llvvObjectProducersUsed", "boosted");
    if(!boostedTauCollHandle.isValid())
    {
      std::cout << "llvvTauCollection Boosted Object NotFound" << std::endl;
    }
    else
    {
      llvvTauCollection boostedTaus = *boostedTauCollHandle;
      //for(size_t i = 0; i < boostedTaus.size(); ++i)
      //  taus.push_back(boostedTaus[i]);
    }

    // Jet Collection
    fwlite::Handle<llvvJetCollection> jetCollHandle;
    jetCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!jetCollHandle.isValid())
    {
      std::cout << "llvvJetCollection Object NotFound" << std::endl;
      continue;
    }
    llvvJetCollection jets_ = *jetCollHandle;
    llvvJetExtCollection jets;
    for(auto i = jetCollHandle->begin(); i != jetCollHandle->end(); ++i)
      jets.push_back(llvvJetExt(*i));

    // MET Collection
    fwlite::Handle<llvvMet> metHandle;
    metHandle.getByLabel(ev, "llvvObjectProducersUsed", "pfMETPFlow");
//    metHandle.getByLabel(ev, "llvvObjectProducersUsed", "pfType1CorrectedMet");
    if(!metHandle.isValid())
    {
      std::cout << "llvvMet Object NotFound" << std::endl;
      continue;
    }
    met = *metHandle;

    // Trigger Prescales
    fwlite::Handle<std::vector<int> > triggerPrescalesHandle;
    triggerPrescalesHandle.getByLabel(ev, "llvvObjectProducersUsed", "triggerPrescales");
    if(!triggerPrescalesHandle.isValid())
    {
      std::cout << "triggerPrescales Object NotFound" << std::endl;
      continue;
    }
    std::vector<int> triggerPrescales = *triggerPrescalesHandle;

    // Rho
    fwlite::Handle<double> rhoHandle;
    rhoHandle.getByLabel(ev, "kt6PFJets", "rho");
    if(!rhoHandle.isValid())
    {
      std::cout << "rho Object NotFound" << std::endl;
      continue;
    }
    rho = *rhoHandle;

    // Rho25
    fwlite::Handle<double> rho25Handle;
    rho25Handle.getByLabel(ev, "kt6PFJetsCentral", "rho");
    if(!rho25Handle.isValid())
    {
      std::cout << "rho25 Object NotFound" << std::endl;
      continue;
    }
    rho25 = *rho25Handle;



    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Finished loading collections, now computing PU weight and trigger scale factor" << std::endl;
    #endif

    // Pileup Weight
    if(isMC)
    {
      if(isStauStau)
      {
        int nEvents = 10000;
        double xsec = stauCrossSec(stauMass, neutralinoMass);
        crossSection = xsec;
        xsecWeight  = xsec/nEvents;
      }
      puWeight     = LumiWeights->weight(genEv.ngenITpu) * PUNorm[0];
      weight       = xsecWeight*puWeight;
      weight_plus  = PuShifters[utils::cmssw::PUUP ]->Eval(genEv.ngenITpu) * (PUNorm[2]/PUNorm[0]);
      weight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(genEv.ngenITpu) * (PUNorm[1]/PUNorm[0]);
    }

    // Get trigger Scale Factor
    triggerSF = 1;
    if(isMC)
    {
      #if defined(DEBUG_EVENT)
      if(debugEvent)
      {
        myCout << " Event";
        if(TauPlusETrigger)
          myCout << ", it is a TauPlusE event";
        if(TauPlusMuTrigger)
          myCout << ", it is a TauPlusMu event";
        myCout << std::endl;

        if(triggeredOn)
        {
          myCout << "  Looping on leptons:" << std::endl;
          for(auto lep = leptons.begin(); lep != leptons.end(); ++lep)
            myCout << "    Lepton (" << lep->id << ", pT=" << lep->pt() << ") trigger bits: " << bitset<8*sizeof(int)>(lep->Tbits) << std::endl;

          myCout << "  Looping on taus:" << std::endl;
          for(auto tau = taus.begin(); tau != taus.end(); ++tau)
            myCout << "    Tau (pT=" << tau->pt() << ") trigger bits: " << bitset<8*sizeof(int)>(tau->Tbits) << std::endl;
        }
      }
      #endif

      if(TauPlusETrigger)
      {
        llvvTau* trigTau = NULL, *leadTau = NULL;
        llvvLepton* trigE = NULL, *leadE = NULL;

        // This is working, but sometimes it can't find the triggered lepton
        //so, when it can't be found, we select the leading pt lepton of the right type
        for(auto lep = leptons.begin(); lep != leptons.end(); ++lep)
        {
          if(lep->Tbits & (3 << 17))
          {
            if(trigE == NULL)
              trigE = &(*lep);
            else
              if(lep->pt() > trigE->pt())
                trigE = &(*lep);
          }

          if(abs(lep->id) == 11)
          {
            if(leadE == NULL)
              leadE = &(*lep);
            else
              if(lep->pt() > leadE->pt())
                leadE = &(*lep);
          }
        }
        if(trigE == NULL)
          trigE = leadE;

        // Tau trigger matching has not yet been enabled in the nTuple production (Tbits is filled with random data)
        //so we use the leading pt pt
        for(auto tau = taus.begin(); tau != taus.end(); ++tau)
        {
          //if(tau->Tbits & (3 << 17))
          //{
          //  trigTau = &(*tau);
          //  break;
          //}

          if(leadTau == NULL)
            leadTau = &(*tau);
          else
            if(tau->pt() > leadTau->pt())
              leadTau = &(*tau);
        }
        if(trigTau == NULL)
          trigTau = leadTau;

        if(trigTau != NULL && trigE != NULL)
        {
          triggerSF *= leptonTauTriggerScaleFactor(*trigE, *trigTau);
        }
        else
        {
          #if defined(DEBUG_EVENT)
          if(debugEvent)
          {
            if(trigE == NULL)
              myCout << " TauPlusE trigSF: Unable to find triggered electron" << std::endl;
            if(trigTau == NULL)
              myCout << " TauPlusE trigSF: Unable to find triggered tau" << std::endl;
          }
          #endif
        }
      }

      if(TauPlusMuTrigger)
      {
        llvvTau* trigTau = NULL, *leadTau = NULL;
        llvvLepton* trigMu = NULL, *leadMu = NULL;

        // This is working, but sometimes it can't find the triggered lepton
        //so, when it can't be found, we select the leading pt lepton of the right type
        for(auto lep = leptons.begin(); lep != leptons.end(); ++lep)
        {
          if(lep->Tbits & (3 << 21))
          {
            if(trigMu == NULL)
              trigMu = &(*lep);
            else
              if(lep->pt() > trigMu->pt())
                trigMu = &(*lep);
          }

          if(abs(lep->id) == 13)
          {
            if(leadMu == NULL)
              leadMu = &(*lep);
            else
              if(lep->pt() > leadMu->pt())
                leadMu = &(*lep);
          }
        }
        if(trigMu == NULL)
          trigMu = leadMu;

        // Tau trigger matching has not yet been enabled in the nTuple production (Tbits is filled with random data)
        //so we use the leading pt pt
        for(auto tau = taus.begin(); tau != taus.end(); ++tau)
        {
          //if(tau->Tbits & (3 << 21))
          //{
          //  trigTau = &(*tau);
          //  break;
          //}

          if(leadTau == NULL)
            leadTau = &(*tau);
          else
            if(tau->pt() > leadTau->pt())
              leadTau = &(*tau);
        }
        if(trigTau == NULL)
          trigTau = leadTau;

        if(trigTau != NULL && trigMu != NULL)
        {
          triggerSF *= leptonTauTriggerScaleFactor(*trigMu, *trigTau);
        }
        else
        {
          #if defined(DEBUG_EVENT)
          if(debugEvent)
          {
            if(trigMu == NULL)
              myCout << " TauPlusMu trigSF: Unable to find triggered muon" << std::endl;
            if(trigTau == NULL)
              myCout << " TauPlusMu trigSF: Unable to find triggered tau" << std::endl;
          }
          #endif
        }
      }

      #if defined(DEBUG_EVENT)
      if(debugEvent)
        myCout << "  Computed trigger SF: " << triggerSF << std::endl;
      #endif
    }
    if(applyScaleFactors && isMC)
      weight *= triggerSF;



    #if defined(DEBUG_EVENT)
    if(debugEvent)
    {
      myCout << " Finished computing PU weight and trigger scale factors" << std::endl;
      myCout << " Getting leading lepton" << std::endl;
    }
    #endif
    // Get Leading Lepton
    for(size_t i = 0; i < leptons.size(); ++i)
    {
      int lepId = abs(leptons[i].id);

      if(lepId == 13 && muCor)
      {
        TLorentzVector p4(leptons[i].px(), leptons[i].py(), leptons[i].pz(), leptons[i].energy());
        muCor->applyPtCorrection(p4, (lepId>0)?1:-1);
        if(isMC)
          muCor->applyPtSmearing(p4, (lepId>0)?1:-1, false);
        leptons[i].SetPxPyPzE(p4.Px(), p4.Py(), p4.Pz(), p4.Energy());
      }

      // Lepton Kinematics
      bool passKin = true;
      bool keepKin = true; // We want to keep leptons with a looser selection than the one for the ID so we can then remove events with multiple leptons
      double eta = (lepId == 11)?(leptons[i].electronInfoRef->sceta):(leptons[i].eta());
      if(lepId == 11) // If Electron
      {
        if(leptons[i].pt() < minElPt)
          passKin = false;
        if(abs(eta) > maxElEta)
          passKin = false;

        if(leptons[i].pt() < 10)
          keepKin = false;
        if(abs(eta) > 2.3)
          keepKin = false;

        if(abs(eta) > ECALGap_MinEta && abs(eta) < ECALGap_MaxEta) // Remove electrons that fall in ECAL Gap
        {
          passKin = false;
          keepKin = false;
        }
      }
      else // If Muon
      {
        if(leptons[i].pt() < minMuPt)
          passKin = false;
        if(abs(eta) > maxMuEta)
          passKin = false;

        if(leptons[i].pt() < 10)
          keepKin = false;
        if(abs(eta) > 2.4)
          keepKin = false;
      }

      // Lepton ID
      bool passID = true, keepID = true;
      Int_t idbits = leptons[i].idbits;
      if(lepId == 11)
      {
        // bool isTight = electronMVAID(leptons[i].electronInfoRef->mvanontrigv0, leptons[i], TightID);
        // bool isLoose = electronMVAID(leptons[i].electronInfoRef->mvanontrigv0, leptons[i], LooseID);
        // bool isLoose = ((idbits >> 4) & 0x1);
        // bool isTight = ((idbits >> 6) & 0x1);
        passID = electronMVAID(leptons[i].electronInfoRef->mvanontrigv0, leptons[i], LooseID);
        keepID = passID;
        if(leptons[i].d0 > 0.045)
          passID = false;
        if(leptons[i].dZ > 0.1)
          passID = false;
        if(leptons[i].d0 > 0.045)
          keepID = false;
        if(leptons[i].dZ > 0.2)
          keepID = false;

        if(leptons[i].electronInfoRef->isConv)
        {
          passID = false;
          keepID = false;
        }
        if(leptons[i].trkLostInnerHits > 0)
        {
          passID = false;
          keepID = false;
        }
      }
      else
      {
        // bool isLoose = ((idbits >> 8) & 0x1);
        // bool isTight = ((idbits >> 10) & 0x1);
        passID = ((idbits >> 10) & 0x1);
        keepID = ((idbits >> 8) & 0x1);
        if(leptons[i].d0 > 0.2)
        {
          passID = false;
          keepID = false;
        }
        if(leptons[i].dZ > 0.5)
        {
          passID = false;
//        if(leptons[i].dZ > 0.2)
          keepID = false;
        }
      }

      // Lepton Isolation
      bool passIso = true, keepIso = true;
      double relIso = utils::cmssw::relIso(leptons[i], rho);
      if(lepId == 11)
      {
        if(relIso > 0.1)
          passIso = false;
        if(relIso > 0.3)
          keepIso = false;
      }
      else
      {
        if(relIso > 0.1)
          passIso = false;
        if(relIso > 0.3)
          keepIso = false;
      }

      // Keep desired leptons
      if(keepKin && keepID && keepIso)
        selLeptons.push_back(leptons[i]);
      if(!triggeredOn)
        continue;

      mon.fillHisto("leptonCutFlow", chTags, 0, weight);
      if(passID)
      {
        mon.fillHisto("leptonCutFlow", chTags, 1, weight);
        if(passKin)
        {
          mon.fillHisto("leptonCutFlow", chTags, 2, weight);
          if(passIso)
            mon.fillHisto("leptonCutFlow", chTags, 3, weight);
        }
      }
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Getting taus" << std::endl;
    #endif
    // Get taus
    for(size_t i = 0; i < taus.size(); ++i)
    {
      llvvTau& tau = taus[i];

      // Tau Kinematics
      bool passKin = true;
      if(tau.pt() < minTauPt)
        passKin = false;
      if(abs(tau.eta()) > maxTauEta)
        passKin = false;

      // Tau overlap with selected leptons
      bool passIso = true;
      for(size_t lep = 0; lep < selLeptons.size(); ++lep)
      {
        int lepId = abs(selLeptons[lep].id);
        if(lepId == 11)
        {
          if(selLeptons[lep].pt() < minElPt)
            continue;
          if(abs(selLeptons[lep].dZ) > 0.1)
            continue;
          double eta = selLeptons[lep].electronInfoRef->sceta;
          if(abs(eta) > maxElEta)
            continue;
          double relIso = utils::cmssw::relIso(selLeptons[lep], rho);
          if(relIso > 0.1)
            continue;
        }
        else
        {
          if(selLeptons[lep].pt() < minMuPt)
            continue;
          if(abs(selLeptons[lep].eta()) > maxMuEta)
            continue;
          double relIso = utils::cmssw::relIso(selLeptons[lep], rho);
          if(relIso > 0.1)
            continue;
          Int_t idbits = selLeptons[lep].idbits;
          bool isTight = ((idbits >> 10) & 0x1);
          if(!isTight)
            continue;
        }
//        if(deltaR(tau, selLeptons[lep]) < 0.1)
        if(deltaR(tau, selLeptons[lep]) < 0.5)
        {
          passIso = false;
          break;
        }
      }

      bool passQual = true;
      if(abs(tau.dZ) > 0.5)
        passQual = false;

      // Tau ID
      bool passID = true;
      if(!tau.passId(llvvTAUID::decayModeFinding)) passID = false;
      if(doDDBkg)
      {
        if(!tau.passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits)) passID = false;
      }
      else
      {
        if(!doTightTauID)
        {
          if(!tau.passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits)) passID = false;
        }
        else
        {
          if(!tau.passId(llvvTAUID::byTightCombinedIsolationDeltaBetaCorr3Hits)) passID = false;
        }
      }
      if(!tau.passId(llvvTAUID::againstMuonTight3)) passID = false;
//      if(!tau.passId(llvvTAUID::againstMuonTight2)) passID = false;
      if(!tau.passId(llvvTAUID::againstElectronMediumMVA5)) passID = false;

      if(passID && passKin && tau.isPF && passIso && passQual)
        selTaus.push_back(tau);
      if(!triggeredOn)
        continue;

      // Fill control histograms
      mon.fillHisto("tauCutFlow", chTags, 0, weight);
      if(tau.isPF)
      {
        mon.fillHisto("tauCutFlow", chTags, 1, weight);
        mon.fillHisto("tauID", chTags, 0, weight);
        if((doTightTauID && tau.passId(llvvTAUID::byTightCombinedIsolationDeltaBetaCorr3Hits)) || (!doTightTauID && tau.passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits)))
        {
          mon.fillHisto("tauID", chTags, 1, weight);
          if(tau.passId(llvvTAUID::decayModeFinding))
          {
            mon.fillHisto("tauID", chTags, 2, weight);
            if(tau.passId(llvvTAUID::againstElectronMediumMVA5))
            {
              mon.fillHisto("tauID", chTags, 3, weight);
              if(tau.passId(llvvTAUID::againstMuonTight3))
                mon.fillHisto("tauID", chTags, 4, weight);
            }
          }
        }
        if(passID)
        {
          mon.fillHisto("tauCutFlow", chTags, 2, weight);
          if(passQual)
          {
            mon.fillHisto("tauCutFlow", chTags, 3, weight);
            if(passKin)
            {
              mon.fillHisto("tauCutFlow", chTags, 4, weight);
              if(passIso)
                mon.fillHisto("tauCutFlow", chTags, 5, weight);
            }
          }
        }
      }
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Getting jets" << std::endl;
    #endif
    // Get Jets
    for(size_t i = 0; i < jets.size(); ++i)
    {
      // Apply jet corrections
      double toRawSF = jets[i].torawsf;
      LorentzVector rawJet(jets[i]*toRawSF);
      jesCor->setJetEta(rawJet.eta());
      jesCor->setJetPt(rawJet.pt());
      jesCor->setJetA(jets[i].area);
      jesCor->setRho(rho);

      if(debugEvent)
        myCout << "   Uncorrected jet pt: " << jets[i].pt() << "; jet eta: " << jets[i].eta() << std::endl;
      double newJECSF(jesCor->getCorrection());
      jets[i].SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
      jets[i] *= newJECSF;
      jets[i].torawsf = 1./newJECSF;

      // Jet ID
      bool passID = true;
      Int_t idbits = jets[i].idbits;
      bool passPFLoose = (idbits & 0x01);
      int fullPuId = (idbits >> 3) & 0x0f;
      bool passLooseFullPuId = ((fullPuId >> 2) & 0x01);
      passID = passLooseFullPuId;

      bool passIso = true;
      for(auto& tau: taus)
      {
        if(deltaR(tau, jets[i]) < 0.4)
          passIso = false;
      }

      // Jet Kinematics
      bool passKin = true;
      if(abs(jets[i].eta()) > maxJetEta)
        passKin = false;
      if(jets[i].pt() <= 30)  // TODO: remove hardcoded value
        passKin = false;

      // B-jets
      bool isBJet = false;
//      bool hasBtagCorr = false;
      if(jets[i].csv > 0.679)
      {
        isBJet = true;
//        hasBtagCorr = true;
      }

      if(isMC)
      {
      }

      // Compute scale and resolution uncertainties
      if(isMC)
      {
        std::vector<float> smearPt = utils::cmssw::smearJER(jets[i].pt(),jets[i].eta(),jets[i].genj.pt());
        jets[i].jer = smearPt[0];
        jets[i].jerup = smearPt[1];
        jets[i].jerdown = smearPt[2];
        if(debugEvent)
          myCout << "   Scaled jet (JES) pt: " << jets[i].pt() << "; jet eta: " << jets[i].eta() << std::endl;
        smearPt = utils::cmssw::smearJES(jets[i].pt(),jets[i].eta(), totalJESUnc);
        jets[i].jesup = smearPt[0];
        jets[i].jesdown = smearPt[1];
      }
      else
      {
        jets[i].jer = jets[i].pt();
        jets[i].jerup = jets[i].pt();
        jets[i].jerdown = jets[i].pt();
        jets[i].jesup = jets[i].pt();
        jets[i].jesdown = jets[i].pt();
      }

      if(passPFLoose && passID && passKin && passIso)
      {
        selJets.push_back(jets[i]);
        selJetsOut.push_back(jets_[i]);
      }
      if(passPFLoose && passID && passKin && isBJet && passIso)
        selBJets.push_back(jets_[i]);
      if(!triggeredOn)
        continue;

      // Fill Jet control histograms
      mon.fillHisto("jetCutFlow", chTags, 0, weight);
      if(passPFLoose)
      {
        mon.fillHisto("jetCutFlow", chTags, 1, weight);
        if(passID)
        {
          mon.fillHisto("jetCutFlow", chTags, 2, weight);
          if(passKin)
          {
            mon.fillHisto("jetCutFlow", chTags, 3, weight);
            if(passIso)
            {
              mon.fillHisto("jetCutFlow", chTags, 4, weight);
              if(isBJet)
              {
                mon.fillHisto("jetCutFlow", chTags, 5, weight);
              }
            }
          }
        }
      }
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Sorting leptons, taus and jets" << std::endl;
    #endif
    if(selLeptons.size() != 0)
    {
      std::sort(selLeptons.begin(), selLeptons.end(), sort_llvvObjectByPt);

//      if(abs(selLeptons[0].id) == 11)
//        chTags.push_back("leadingE");
//      else
//        chTags.push_back("leadingMu");
    }

    if(selTaus.size() != 0)
      std::sort(selTaus.begin(), selTaus.end(), sort_llvvObjectByPt);

    nBJets = selBJets.size();
    nJets = selJets.size();
    if(nJets != 0)
      std::sort(selJets.begin(), selJets.end(), sort_llvvObjectByPt);



    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Requiring an opposite sign pair" << std::endl;
    #endif
    // Opposite Sign requirements
    maxPtSum = 0;
    tauIndex = -1;
    leptonIndex = -1;
    for(size_t i = 0; i < selLeptons.size(); ++i)
    {
      if(abs(selLeptons[i].id) == 11) // Electron
      {
        if(selLeptons[i].pt() < minElPt)
          continue;
        if(abs(selLeptons[i].dZ) > 0.1)
          continue;
        double eta = selLeptons[i].electronInfoRef->sceta;
        if(abs(eta) > maxElEta)
          continue;
        double relIso = utils::cmssw::relIso(selLeptons[i], rho);
        if(relIso > 0.1)
          continue;
      }
      else
      {
        if(selLeptons[i].pt() < minMuPt)
          continue;
        if(abs(selLeptons[i].eta()) > maxMuEta)
          continue;
        double relIso = utils::cmssw::relIso(selLeptons[i], rho);
        if(relIso > 0.1)
          continue;
        Int_t idbits = selLeptons[i].idbits;
        bool isTight = ((idbits >> 10) & 0x1);
        if(!isTight)
          continue;
      }

      for(size_t j = 0; j < selTaus.size(); ++j)
      {
        double PtSum = selLeptons[i].pt() + selTaus[j].pt();
        if(selLeptons[i].id * selTaus[j].id < 0)
        {
          if(PtSum > maxPtSum)
          {
            maxPtSum = PtSum;
            tauIndex = j;
            leptonIndex = i;
            isOS = true;
            tauLeadPt = selTaus[tauIndex].pt();
            lepLeadPt = selLeptons[leptonIndex].pt();
//            triggerSF = leptonTauTriggerScaleFactor(selLeptons[leptonIndex], selTaus[tauIndex]);
            leptonIdIsoSF = leptonIdAndIsoScaleFactor(selLeptons[leptonIndex]);
            tauSF = ::tauSF(selTaus[tauIndex], gen, antiEMva5Medium);
          }
        }
        if(PtSum < maxPtSum) // Skip a few iterations if it is not expected that we will find a better candidate
          break;
      }
    }

    if(isOS && isMC)
    {
      for(auto & genPart : gen)
      {
        if(genPart.status == 3)
        {
          if(genPart.id == selLeptons[leptonIndex].id)
          {
            if(deltaR(selLeptons[leptonIndex], genPart) < 0.3)
              isPromptLep = true;
          }
          if(genPart.id == selTaus[tauIndex].id)
          {
            if(deltaR(selTaus[tauIndex], genPart) < 0.3)
              isPromptTau = true;
          }
        }
      }
    }
    if(!isMC)
    {
//      triggerSF = 1;
      leptonIdIsoSF = 1;
      tauSF = 1;
    }
    if(isOS && applyScaleFactors)
    {
//      weight *= triggerSF;
      weight *= leptonIdIsoSF;
      weight *= tauSF;
    }

    if(isOS && doDDBkg)
    {
      double taupt  = selTaus[tauIndex].pt();
      double taueta = selTaus[tauIndex].eta();
      bool isetau = abs(selLeptons[leptonIndex].id) == 11;

      TH1* FRhist = mutauFR;
      TH1* PRhist = mutauPR;
      if(isetau)
      {
        FRhist = etauFR;
        PRhist = etauPR;
      }


      Int_t bin = -1;
      if(FRhist->InheritsFrom("TH2"))
      {
        if(taupt < FRhist->GetXaxis()->GetXmin())
          taupt = FRhist->GetXaxis()->GetXmin();
        if(taupt > FRhist->GetXaxis()->GetXmax())
          taupt = FRhist->GetXaxis()->GetXmax();
        if(taueta < FRhist->GetYaxis()->GetXmin())
          taueta = FRhist->GetYaxis()->GetXmin();
        if(taueta > FRhist->GetYaxis()->GetXmax())
          taueta = FRhist->GetYaxis()->GetXmax();

        bin = FRhist->FindBin(taupt, taueta);
      }
      else
      {
        if(taueta < FRhist->GetXaxis()->GetXmin())
          taueta = FRhist->GetXaxis()->GetXmin();
        if(taueta > FRhist->GetXaxis()->GetXmax())
          taueta = FRhist->GetXaxis()->GetXmax();

        bin = FRhist->FindBin(taueta);
      }

      // Fake Rate
      double fakeRate = FRhist->GetBinContent(bin);
//      double promptRate = PRhist->GetBinContent(bin);// */
//      double fakeRate = 0.497205; //\pm0.00413627
      double promptRate = 0.783028; //\pm0.00379267

      if(selTaus[tauIndex].passId(llvvTAUID::byTightCombinedIsolationDeltaBetaCorr3Hits))
      {
        weight *= ((promptRate - 1) * fakeRate) / (promptRate - fakeRate);
      }
      else
      {
        weight *= (promptRate * fakeRate) / (promptRate - fakeRate);
      }
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Rejecting event if multilepton" << std::endl;
    #endif
    // Reject events with more leptons
    isMultilepton = false;
    if(isOS)
    {
      for(size_t i = 0; i < selLeptons.size(); ++i)
      {
        if(i == (size_t)leptonIndex)
          continue;
        if(abs(selLeptons[i].id) != 11 && selLeptons[i].dZ > 0.2)  // If muon
          continue;
        isMultilepton = true;
        break;
      }

      // Set up channels
      if(!isMultilepton)
      {
        if(abs(selLeptons[leptonIndex].id) == 11)
        {
          chTags.push_back("etau");
          isetau = true;
        }
        else
        {
          chTags.push_back("mutau");
          ismutau = true;
        }

        #if defined(DEBUG_EVENT)
        if(debugEvent || true)
        {
          bool isMuTau = false;
          bool isAll = false;
          for(auto tag = chTags.begin(); tag != chTags.end(); ++tag)
          {
            if(*tag == "all")
              isAll = true;
            if(*tag == "mutau")
              isMuTau = true;
          }

          if(!isAll)
          {
            myCout << " Found a ";
            if(isMuTau)
              myCout << "mu-tau";
            else
              myCout << "e-tau";
            myCout << " event";
            if(isAll)
              myCout << " and it is tagged as all" << std::endl;
            else
              myCout << " but it is not tagged as all" << std::endl;
          }
        }
        #endif
      }
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Computing angular variables" << std::endl;
    #endif
    // Get angular variables
    if(isOS && !isMultilepton)
    {
      /**    LAB FRAME    **/
      TLorentzVector lep(selLeptons[leptonIndex].Px(), selLeptons[leptonIndex].Py(), selLeptons[leptonIndex].Pz(), selLeptons[leptonIndex].E());
      TLorentzVector tau(selTaus[tauIndex].Px(), selTaus[tauIndex].Py(), selTaus[tauIndex].Pz(), selTaus[tauIndex].E());
      TLorentzVector Tmet(met.Px(), met.Py(), met.Pz(), met.E());

      deltaAlphaLepTau = lep.Angle(tau.Vect());
      deltaRLepTau = deltaR(selTaus[tauIndex], selLeptons[leptonIndex]);
      deltaPhiLepTauMET = Tmet.DeltaPhi(lep + tau);
      deltaPhiLepTau = lep.DeltaPhi(tau);

      minDeltaPhiMETJetPt40 = 5;
      for(auto &jet : selJets)
      {
        TLorentzVector tJet(jet.Px(), jet.Py(), jet.Pz(), jet.E());

        if(tJet.Pt() < 40)
          break;

        double temp = Tmet.DeltaPhi(tJet);
        if(temp < minDeltaPhiMETJetPt40)
          minDeltaPhiMETJetPt40 = temp;
      }

      double posSign = Tmet.CosTheta();
      cosThetaTau = tau.CosTheta();
      cosThetaLep = lep.CosTheta();
      if(posSign < 0)
      {
        cosThetaTau = -cosThetaTau;
        cosThetaLep = -cosThetaLep;
      }

      /**   CS FRAME   **/
      TLorentzVector tauSystem = tau + lep;
      TLorentzVector tauCS = tau, lepCS = lep, metCS = Tmet;
      TLorentzVector beam1(0, 0, 0, 0);
      TLorentzVector beam2(0, 0, 0, 0);
      double energy = sqrtS * 500.; // sqrtS / 2 * 1000
      double mom = sqrt(energy*energy + 0.938*0.938);
      if(posSign > 0)
      {
        beam1.SetPxPyPzE(0, 0,  mom, energy);
        beam2.SetPxPyPzE(0, 0, -mom, energy);
      }
      else
      {
        beam1.SetPxPyPzE(0, 0, -mom, energy);
        beam2.SetPxPyPzE(0, 0,  mom, energy);
      }

      TVector3 boost = -tauSystem.BoostVector();
      //tauSystem.Boost(boost); // By construction, this will be 0
      tauCS.Boost(boost);
      lepCS.Boost(boost);
      metCS.Boost(boost);
      beam1.Boost(boost);
      beam2.Boost(boost);

      TLorentzVector SQA = beam1 - beam2; // Spin quantization axis
      TRotation rotation;

      TVector3 newZAxis = SQA.Vect().Unit();
      TVector3 targetZaxis(0, 0, 1);
      TVector3 rotAxis = targetZaxis.Cross(newZAxis);
      double rotAngle = targetZaxis.Angle(newZAxis);
      rotation.Rotate(-rotAngle, rotAxis);

      SQA.Transform(rotation);
      tauCS.Transform(rotation);
      lepCS.Transform(rotation);
      metCS.Transform(rotation);
      beam1.Transform(rotation);
      beam2.Transform(rotation);

      cosThetaCS = lepCS.CosTheta();
      deltaPhiLepMETCS = metCS.DeltaPhi(lepCS);
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Computing pair mass" << std::endl;
    #endif
    // Tau-Lepton pair mass calculation
    isSVfit = doSVfit;
    if(isOS && !isMultilepton)
    {
      auto selLepton = selLeptons[leptonIndex];
      auto selTau    = selTaus[tauIndex];
      invMass = (selLepton+selTau).M(); // Invariant mass

      if(doSVfit)
      {
        TMatrixD covMET(2, 2);
        std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
        svFitStandalone::Vector measuredMET(met.px(), met.py(), 0);

        covMET[0][0] = met.sigx2;
        covMET[0][1] = met.sigxy;
        covMET[1][0] = met.sigxy;
        covMET[1][1] = met.sigy2;

        if(abs(selLepton.id) == 11)
          measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, svFitStandalone::LorentzVector(selLepton.px(), selLepton.py(), selLepton.pz(), selLepton.E())));
        else
          measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, svFitStandalone::LorentzVector(selLepton.px(), selLepton.py(), selLepton.pz(), selLepton.E())));
        measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, svFitStandalone::LorentzVector(selTau.px(), selTau.py(), selTau.pz(), selTau.E())));

        SVfitStandaloneAlgorithm SVfit_algo(measuredTauLeptons, measuredMET, covMET, 0);
        //SVfit_algo.maxObjFunctionCalls(10000) // To change the max number of iterations before minimization is terminated, default 5000
        SVfit_algo.addLogM(false); // To not use the LogM penalty, it is used by default
        //SVfit_algo.metPower(0.5); // Additional power to enhance MET likelihood, default is 1.
        //SVfit_algo.fit();
        //SVfit_algo.integrate();
        SVfit_algo.integrateVEGAS();
        //SVfit_algo.integrateMarkovChain();
        if(SVfit_algo.isValidSolution())
          mass = SVfit_algo.mass();
        else
          isSVfit = false;
      }
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Computing MT and deconstructed MT" << std::endl;
    #endif
    // MT and deconstructed MT calculation
    if(isOS && !isMultilepton && (!doSVfit || isSVfit))
    {
      auto& selLepton = selLeptons[leptonIndex];
      double cosDeltaPhi = cos(deltaPhi(selLepton.phi(), met.phi()));
      double fac = 2 * met.pt() * selLepton.pt();

      mt = sqrt(fac * (1 - cosDeltaPhi));
      Q80 = 1 - (80.0*80.0) / fac;
      Q100 = 1 - (100.0*100.0) / fac;
      cosPhi = cosDeltaPhi;

      auto& selTau = selTaus[tauIndex];
      cosDeltaPhi = cos(deltaPhi(selTau.phi(), met.phi()));
      fac = 2 * met.pt() * selTau.pt();

      mtTau = sqrt(fac * (1 - cosDeltaPhi));
      Q80Tau = 1 - (80.0*80.0) / fac;
      Q100Tau = 1 - (100.0*100.0) / fac;
      cosPhiTau = cosDeltaPhi;

      sumMt = mt + mtTau;
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Computing MT2" << std::endl;
    #endif
    // MT2 calculation
    if(isOS && !isMultilepton && (!doSVfit || isSVfit))
    {
      auto selLepton = selLeptons[leptonIndex];
      auto selTau    = selTaus[tauIndex];
      double pa[3];
      double pb[3];
      double pmiss[3];
      double mn;

      pa[0] = selLepton.M();
      pa[1] = selLepton.px();
      pa[2] = selLepton.py();
      pb[0] = selTau.M();
      pb[1] = selTau.px();
      pb[2] = selTau.py();
      pmiss[0] = 0;
      pmiss[1] = met.px();
      pmiss[2] = met.py();
      mn = 0;

      mt2_bisect::mt2 mt2_event;
      mt2_event.set_momenta(pa,pb,pmiss);
      mt2_event.set_mn(mn);
      mt2 = mt2_event.get_mt2();

/*      mn = 50;
      mt2_event.set_momenta(pa,pb,pmiss);
      mt2_event.set_mn(mn);
      mt2_50 = mt2_event.get_mt2();

      mn = 150;
      mt2_event.set_momenta(pa,pb,pmiss);
      mt2_event.set_mn(mn);
      mt2_150 = mt2_event.get_mt2(); // */
    }

    bool stauPlot = false;
    if(stauMass == stauMtoPlot && neutralinoMass == neutralinoMtoPlot)
      stauPlot = true;

    bool isIPM = false;
/*    if(invMass > 15 && (invMass < 45 || invMass > 75) && mt2 > 40 && minDeltaPhiMETJetPt40 > 1)
    {
      isIPM = true;
      chTags.push_back("IPM");
      if(abs(selLeptons[leptonIndex].id) == 11)
        chTags.push_back("IPM-etau");
      else
        chTags.push_back("IPM-mutau");
    }// */

    bool keep = false;
    if(doDDBkg)
      keep = true;
    if(!isMC)
      keep = true;
    if(isPromptTau)
      keep = true;
    if(isStauStau)
      keep = true;

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Filling histograms" << std::endl;
    #endif
    bool plotThisEvent = !isStauStau || stauPlot;
    if(plotThisEvent)
    {
      //mon.fillHisto("nvtxAll", chTags, nvtx, weight);
      if(triggeredOn && keep)
      {
        mon.fillHisto("eventflow", chTags, 0, weight);
        if(met.pt() > 30)
//        if(true)
        {
          mon.fillHisto("eventflow", chTags, 1, weight);
          if(selLeptons.size() > 0)
          {
            mon.fillHisto("eventflow", chTags, 2, weight);
            mon.fillHisto("nbjets", chTags, selBJets.size(), weight);
            if(selTaus.size() > 0)
            {
              mon.fillHisto("eventflow", chTags, 3, weight);
              if(selBJets.size() == 0)
              {
                mon.fillHisto("eventflow", chTags, 4, weight);
                if(isOS)
                {
                  mon.fillHisto("eventflow", chTags, 5, weight);
                  if(!isMultilepton)
                  {
                    mon.fillHisto("eventflow", chTags, 6, weight);
                    if(!doSVfit || isSVfit)
                    {
                      if(keep)
                      {
                      mon.fillHisto("eventflow", chTags, 7, weight);

                      mon.fillHisto("nvtx", chTags, nvtx, weight);
                      mon.fillHisto("nvtxraw", chTags, nvtx, weight/puWeight);
                      mon.fillHisto("nup", "", genEv.nup, 1);

                      mon.fillHisto("rho", chTags, rho, weight);
                      mon.fillHisto("rho25", chTags, rho25, weight);

                      mon.fillHisto("MET", chTags, met.pt(), weight);

                      mon.fillHisto("MT", chTags, mt, weight);
                      mon.fillHisto("Q80", chTags, Q80, weight);
                      mon.fillHisto("Q100", chTags, Q100, weight);
                      mon.fillHisto("cosPhi", chTags, cosPhi, weight);
                      mon.fillHisto("Q80VsCosPhi", chTags, cosPhi, Q80, weight);
                      mon.fillHisto("Q100VsCosPhi", chTags, cosPhi, Q100, weight);

                      mon.fillHisto("MTTau", chTags, mtTau, weight);
                      mon.fillHisto("Q80Tau", chTags, Q80Tau, weight);
                      mon.fillHisto("Q100Tau", chTags, Q100Tau, weight);
                      mon.fillHisto("cosPhiTau", chTags, cosPhiTau, weight);
                      mon.fillHisto("Q80VsCosPhiTau", chTags, cosPhiTau, Q80Tau, weight);
                      mon.fillHisto("Q100VsCosPhiTau", chTags, cosPhiTau, Q100Tau, weight);

                      mon.fillHisto("SumMT", chTags, sumMt, weight);

                      mon.fillHisto("MT2", chTags, mt2, weight);
                      if(doSVfit)
                        mon.fillHisto("SVFitMass", chTags, mass, weight);
                      mon.fillHisto("InvMass", chTags, invMass, weight);

                      mon.fillHisto("deltaAlphaLepTau", chTags, deltaAlphaLepTau, weight);
                      mon.fillHisto("deltaRLepTau", chTags, deltaRLepTau, weight);
                      mon.fillHisto("deltaPhiLepTauMET", chTags, deltaPhiLepTauMET, weight);
                      mon.fillHisto("deltaPhiLepTau", chTags, deltaPhiLepTau, weight);
                      mon.fillHisto("cosThetaTau", chTags, cosThetaTau, weight);
                      mon.fillHisto("cosThetaLep", chTags, cosThetaLep, weight);
                      mon.fillHisto("cosThetaCS", chTags, cosThetaCS, weight);
                      mon.fillHisto("deltaPhiLepMETCS", chTags, deltaPhiLepMETCS, weight);
                      mon.fillHisto("minDeltaPhiMETJetPt40", chTags, minDeltaPhiMETJetPt40, weight);

                      mon.fillHisto("metVsPtl", chTags, selLeptons[leptonIndex].pt(), met.pt(), weight);
                      mon.fillHisto("metVsPtTau", chTags, selTaus[tauIndex].pt(), met.pt(), weight);
                      mon.fillHisto("metPtVsmetEt", chTags, met.pt(), met.Et(), weight);

                      mon.fillHisto("nlep", chTags, selLeptons.size(), weight);
                      double eta = selLeptons[leptonIndex].eta();
                      if(abs(selLeptons[leptonIndex].id) == 11) eta = selLeptons[leptonIndex].electronInfoRef->sceta;
                      mon.fillHisto("etaSelectedLep", chTags, eta, weight);
                      mon.fillHisto("ptSelectedLep", chTags, selLeptons[leptonIndex].pt(), weight);
                      mon.fillHisto("chargeSelectedLep", chTags, (selLeptons[leptonIndex].id > 0)?(-1):(1), weight);

                      mon.fillHisto("ntaus", chTags, selTaus.size(), weight);
                      mon.fillHisto("ptSelectedTau", chTags, selTaus[tauIndex].pt(), weight);
                      mon.fillHisto("ptSelectedTauExtended", chTags, selTaus[tauIndex].pt(), weight);
                      mon.fillHisto("etaSelectedTau", chTags, selTaus[tauIndex].eta(), weight);
                      mon.fillHisto("chargeSelectedTau", chTags, (selTaus[tauIndex].id > 0)?(-1):(1), weight);
                      mon.fillHisto("emfracSelectedTau", chTags, selTaus[tauIndex].emfraction, weight);
                      mon.fillHisto("dzSelectedTau", chTags, selTaus[tauIndex].dZ, weight);

                      mon.fillHisto("njets", chTags, selJets.size(), weight);
                      if(selJets.size() != 0)
                      {
                        mon.fillHisto("jetleadpt", chTags, selJets[0].pt(), weight);
                        mon.fillHisto("jetleadeta", chTags, selJets[0].eta(), weight);
                      }

                      for(auto i = selJets.begin(); i != selJets.end(); ++i)
                      {
                        mon.fillHisto("jetcsv", chTags, i->origcsv, weight);
                      }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if(triggeredOn && met.pt() > 30 && selLeptons.size() > 0 && selBJets.size() == 0 && selTaus.size() > 0 && isOS && !isMultilepton && (!doSVfit || isSVfit))
//    if(triggeredOn && selLeptons.size() > 0 && selBJets.size() == 0 && selTaus.size() > 0 && isOS && !isMultilepton && (!doSVfit || isSVfit))
      selected = true;

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Filling TTree" << std::endl;
    NAN_WARN(rho)
    NAN_WARN(rho25)
    NAN_WARN(weight)
    NAN_WARN(weight_plus)
    NAN_WARN(weight_minus)
    NAN_WARN(puWeight)
    NAN_WARN(triggerSF)
    NAN_WARN(leptonIdIsoSF)
    NAN_WARN(tauSF)
    NAN_WARN(mass)
    NAN_WARN(invMass)
    NAN_WARN(mt)
    NAN_WARN(Q80)
    NAN_WARN(Q100)
    NAN_WARN(cosPhi)
    NAN_WARN(mt2)
    NAN_WARN(stauMass)
    NAN_WARN(neutralinoMass)
    NAN_WARN(deltaAlphaLepTau)
    NAN_WARN(deltaRLepTau)
    NAN_WARN(deltaPhiLepTauMET)
    NAN_WARN(deltaPhiLepTau)
    NAN_WARN(cosThetaTau)
    NAN_WARN(cosThetaLep)
    NAN_WARN(deltaPhiLepMETCS)
    NAN_WARN(cosThetaCS)
//    NAN_WARN(minDeltaPhiMetJet40)
    NAN_WARN(tauLeadPt)
    NAN_WARN(lepLeadPt)
    NAN_WARN(maxPtSum)
    #endif
//    if(saveSummaryTree)
    bool outputEvent = true;
//    if(isMC && !isStauStau)
//      if(!isPromptTau)
//        outputEvent = false;
    if(saveSummaryTree && outputEvent)// */
    {
      TDirectory* cwd = gDirectory;
      summaryOutFile->cd();
      summaryTree->Fill();
      cwd->cd();
    }
    
    
    if(outputEventList)
    {
      eventListFile << std::setw(EVENTLISTWIDTH) << ev.eventAuxiliary().run() << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << ev.eventAuxiliary().luminosityBlock() << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << ev.eventAuxiliary().event() << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << (selected?("True"):("False")) << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << weight << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << puWeight << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << xsecWeight << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << triggerSF << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << leptonIdIsoSF << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << tauSF << "|";
      if(doDDBkg)
        eventListFile << std::setw(EVENTLISTWIDTH) << "NA" << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << met.pt() << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << jets.size() << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << selLeptons.size() << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << selTaus.size() << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << selJets.size() << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << selBJets.size() << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << mt << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << mtTau << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << sumMt << "|";
      eventListFile << std::setw(EVENTLISTWIDTH) << mt2 << "|";
      eventListFile << std::endl;
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Finished event " << iev << std::endl;
    #endif
//    break;
    if(limit != 0)
      if(iev >= limit - 1)
        break;
  }

  // Output temporary buffer and restore cout and cerr behaviour
  std::cout.rdbuf(coutbuf);
  std::cerr.rdbuf(cerrbuf);
  std::cout << std::endl;
//  std::cout << buffer.str();

  std::cout << "totalEntries: " << totalEntries << "; vs nInitEvent: " << nInitEvent << ";" << std::endl;



  /***************************************************************************/
  /*                        Saving Histograms to File                        */
  /***************************************************************************/
  std::cout << "Saving results in " << outUrl << std::endl;
  TFile* outfile = new TFile(outUrl.c_str(), "RECREATE");
  mon.Write();
  outfile->Close();
  delete outfile;

  if(saveSummaryTree)
  {
    TDirectory* cwd = gDirectory;
    summaryOutFile->cd();
    summaryTree->Write();
    summaryOutFile->Close();
    delete summaryOutFile;
    cwd->cd();
  }
  
  if(outputEventList)
  {
    eventListFile.close();
  }

  return 0;
}

double stauCrossSec(double stauM, double neutM)
{
  //Points taken from  http://arxiv.org/abs/1204.2379
  //Mstau == 100 => 0.1
  //Mstau == 125 => 0.05
  //Mstau == 145 => 0.03
  //Mstau == 195 => 0.01
  //Mstau == 240 => 0.005
  //Mstau == 275 => 0.003
  //Mstau == 300 => 0.002
  //Mstau == 360 => 0.001
  //Mstau == 425 => 0.0005
  double a = 0.2979;
  double b = 17.626;
  double c = 67.632;
  double d = 3.463;
  return a / (1 + std::pow((stauM - b) / c, d));
}

bool electronMVAID(double mva, llvvLepton& lepton, ID_Type id)
{
  if(id == MediumID)
    id = TightID;

  double eta = lepton.electronInfoRef->sceta;
  bool pass = false;

  switch(id)
  {
  case LooseID:
    if(lepton.pt() < 20)
    {
      if(abs(eta) < 0.8)
      {
        if(mva > 0.925)
          pass = true;
      }
      else
      {
        if(abs(eta) < 1.479)
        {
          if(mva > 0.915)
            pass = true;
        }
        else
        {
          if(mva > 0.965)
            pass = true;
        }
      }
    }
    else
    {
      if(abs(eta) < 0.8)
      {
        if(mva > 0.905)
          pass = true;
      }
      else
      {
        if(abs(eta) < 1.479)
        {
          if(mva > 0.955)
            pass = true;
        }
        else
        {
          if(mva > 0.975)
            pass = true;
        }
      }
    }
    break;
  case TightID:
  default:
    if(lepton.pt() >= 20)
    {
      if(abs(eta) < 0.8)
      {
        if(mva > 0.925)
          pass = true;
      }
      else
      {
        if(abs(eta) < 1.479)
        {
          if(mva > 0.975)
            pass = true;
        }
        else
        {
          if(mva > 0.985)
            pass = true;
        }
      }
    }
    break;
  }

  return pass;
}

double tauSF(llvvTau& tau, llvvGenParticleCollection& genPartColl, TAU_E_ID eId)
{
  double scaleFactor = 1;

  // No correction necessary for tau ID
  // No correction necessary for normalization of Jet->Tau fake (if doing shape analysis, should investigate pt dependence of this)
  // No correction necessary for mu->Tau fake if using tight muon discriminator
  // Hadronic tau energy scale, no correction necessary
  // Tau charge misidentification rate, no correction necessary
  // This leaves only e->Tau fake, which must be corrected according to the anti-e discriminator used

  bool isElectronFakingTau = false;

  for(auto genPart = genPartColl.begin(); genPart != genPartColl.end(); ++genPart)
  {
    if(abs(genPart->id) == 11 && genPart->status == 3) // If the gen particle is a stable electron
    {
      if(deltaR(tau, *genPart) < 0.3)
      {
        isElectronFakingTau = true;
        break;
      }
    }
  }

  if(isElectronFakingTau)
  {
    double barrelSF = 1;
    double endcapSF = 1;

    switch(eId)
    {
    case antiELoose: // Both are 1, so no change
      break;
    case antiEMedium:
      barrelSF = 0.95;
      endcapSF = 0.75;
      break;
    case antiETight:
      barrelSF = 0.90;
      endcapSF = 0.70;
      break;
    case antiEMva:
      barrelSF = 0.85;
      endcapSF = 0.65;
      break;
    case antiEMva3Loose:
      barrelSF = 1.4; // +- 0.3
      endcapSF = 0.8; // +- 0.3
      break;
    case antiEMva3Medium:
      barrelSF = 1.6; // +- 0.3
      endcapSF = 0.8; // +- 0.3
      break;
    case antiEMva3Tight:
      barrelSF = 2.0; // +- 0.4
      endcapSF = 1.2; // +- 0.4
      break;
    case antiEMva3VTight:
      barrelSF = 2.4; // +- 0.5
      endcapSF = 1.2; // +- 0.5
      break;
    case antiEMva5Medium: // 1.6 +/- 0.3 for the barrel (abs(tauEta) < 1.5) and 1.1 +/- 0.3 for the endcap.
    default:
      barrelSF = 1.6;
      endcapSF = 1.1;
      break;
    }

    if(tau.eta() < 1.5)
    {
      scaleFactor = barrelSF;
    }
    else
    {
      scaleFactor = endcapSF;
    }
  }

  return scaleFactor;
}

double leptonIdAndIsoScaleFactor(llvvLepton& lepton)
{
  double scaleFactor = 1;

  if(abs(lepton.id) == 11) // If an electron
  {
    double isoSF = 0;
    double idSF  = 0;

    double pt = lepton.pt();
    double eta = lepton.electronInfoRef->sceta;
    if(abs(eta) < 1.479)  // Electron in barrel
    {
      if(pt < 30)
      {
        idSF  = 0.8999; // +- 0.0018
        isoSF = 0.9417; // +- 0.0019
      }
      else
      {
        idSF  = 0.9486; // +- 0.0003
        isoSF = 0.9804; // +- 0.0003
      }
    }
    else // Electron in endcap
    {
      if(pt < 30)
      {
        idSF  = 0.7945; // +- 0.0055
        isoSF = 0.9471; // +- 0.0037
      }
      else
      {
        idSF  = 0.8866; // +- 0.0001
        isoSF = 0.9900; // +- 0.0002
      }
    }

    scaleFactor = isoSF * idSF;
  }
  else // If a muon
  {
    double isoSF = 0;
    double idSF  = 0;

    double eta = lepton.eta();
    double pt  = lepton.pt();
    if(abs(eta) < 0.8) // Barrel muons
    {
      if(pt < 30)
      {
        idSF  = 0.9818; // +- 0.0005
        isoSF = 0.9494; // +- 0.0015
      }
      else
      {
        idSF  = 0.9852; // +- 0.0001
        isoSF = 0.9883; // +- 0.0003
      }
    }
    else
    {
      if(abs(eta) < 1.2) // Transition muons
      {
        if(pt < 30)
        {
          idSF  = 0.9829; // +- 0.0009
          isoSF = 0.9835; // +- 0.0020
        }
        else
        {
          idSF  = 0.9852; // +- 0.0002
          isoSF = 0.9937; // +- 0.0004
        }
      }
      else // Endcap muons
      {
        if(pt < 30)
        {
          idSF  = 0.9869; // +- 0.0007
          isoSF = 0.9923; // +- 0.0013
        }
        else
        {
          idSF  = 0.9884; // +- 0.0001
          isoSF = 0.9996; // +- 0.0005
        }
      }
    }

    scaleFactor = isoSF * idSF;
  }

  return scaleFactor;
}

double leptonTauTriggerScaleFactor(llvvLepton& lepton, llvvTau& tau)
{
  double scaleFactor = 1;
  double m0[2], sigma[2], alpha[2], n[2], norm[2]; // Index 0 - Data; Index 1 - MC
  double pt, eta;

  if(abs(lepton.id) == 11) // eTau channel
  {
    // Electron leg
    eta = lepton.electronInfoRef->sceta;
    pt  = lepton.pt();
    if(abs(eta) < 1.479) // In barrel
    {
      m0[0]    = 22.9704;
      m0[1]    = 21.7243;
      sigma[0] = 1.0258;
      sigma[1] = 0.619015;
      alpha[0] = 1.26889;
      alpha[1] = 0.739301;
      n[0]     = 1.31024;
      n[1]     = 1.34903;
      norm[0]  = 1.06409;
      norm[1]  = 1.02594;
    }
    else // In endcap
    {
      m0[0] = 21.9816;
      m0[1] = 22.1217;
      sigma[0] = 1.40993;
      sigma[1] = 1.34054;
      alpha[0] = 0.978597;
      alpha[1] = 1.8885;
      n[0] = 2.33144;
      n[1] = 1.01855;
      norm[0] = 0.937552;
      norm[1] = 4.7241;
    }

    double electronSF = 1;
    if(pt >= 20) // Do not apply for electrons with pt below threshold
    {
      double electronDataEff = efficiency(pt, m0[0], sigma[0], alpha[0], n[0], norm[0]);
      double electronMCEff   = efficiency(pt, m0[1], sigma[1], alpha[1], n[1], norm[1]);
      electronSF = electronDataEff/electronMCEff;
    }

    // Tau leg
    eta = tau.eta();
    pt  = tau.pt();
    if(abs(eta) < 1.5) // In barrel
    {
      m0[0]    = 18.538229;
      m0[1]    = 18.605055;
      sigma[0] = 0.651562;
      sigma[1] = 0.264062;
      alpha[0] = 0.324869;
      alpha[1] = 0.139561;
      n[0]     = 13.099048;
      n[1]     = 4.792849;
      norm[0]  = 0.902365;
      norm[1]  = 0.915035;
    }
    else // In endcap
    {
      m0[0]    = 18.756548;
      m0[1]    = 18.557810;
      sigma[0] = 0.230732;
      sigma[1] = 0.280908;
      alpha[0] = 0.142859;
      alpha[1] = 0.119282;
      n[0]     = 3.358497;
      n[1]     = 17.749043;
      norm[0]  = 0.851919;
      norm[1]  = 0.865756;
    }

    double tauSF = 1;
    if(pt >= 20)
    {
      double tauDataEff = efficiency(pt, m0[0], sigma[0], alpha[0], n[0], norm[0]);
      double tauMCEff   = efficiency(pt, m0[1], sigma[1], alpha[1], n[1], norm[1]);
      tauSF = tauDataEff/tauMCEff;
    }

    scaleFactor = electronSF * tauSF;
  }
  else // muTau channel
  {
    // Muon leg
    eta = lepton.eta();
    pt  = lepton.pt();
    if(eta < -1.2)
    {
      m0[0]    = 15.9977;
      m0[1]    = 16.0051;
      sigma[0] = 7.64004e-05;
      sigma[1] = 2.45144e-05;
      alpha[0] = 6.4951e-08;
      alpha[1] = 4.3335e-09;
      n[0]     = 1.57403;
      n[1]     = 1.66134;
      norm[0]  = 0.865325;
      norm[1]  = 0.87045;
    }
    else if(eta < -0.8)
    {
      m0[0]    = 17.3974;
      m0[1]    = 17.3135;
      sigma[0] = 0.804001;
      sigma[1] = 0.747636;
      alpha[0] = 1.47145;
      alpha[1] = 1.21803;
      n[0]     = 1.24295;
      n[1]     = 1.40611;
      norm[0]  = 0.928198;
      norm[1]  = 0.934983;
    }
    else if(eta < 0)
    {
      m0[0]    = 16.4307;
      m0[1]    = 15.9556;
      sigma[0] = 0.226312;
      sigma[1] = 0.0236127;
      alpha[0] = 0.265553;
      alpha[1] = 0.00589832;
      n[0]     = 1.55756;
      n[1]     = 1.75409;
      norm[0]  = 0.974462;
      norm[1]  = 0.981338;
    }
    else if(eta < 0.8)
    {
      m0[0]    = 17.313;
      m0[1]    = 15.9289;
      sigma[0] = 0.662731;
      sigma[1] = 0.0271317;
      alpha[0] = 1.3412;
      alpha[1] = 0.00448573;
      n[0]     = 1.05778;
      n[1]     = 1.92101;
      norm[0]  = 1.26624;
      norm[1]  = 0.978625;
    }
    else if(eta < 1.2)
    {
      m0[0]    = 16.9966;
      m0[1]    = 16.5678;
      sigma[0] = 0.550532;
      sigma[1] = 0.328333;
      alpha[0] = 0.807863;
      alpha[1] = 0.354533;
      n[0]     = 1.55402;
      n[1]     = 1.67085;
      norm[0]  = 0.885134;
      norm[1]  = 0.916992;
    }
    else
    {
      m0[0]    = 15.9962;
      m0[1]    = 15.997;
      sigma[0] = 0.000106195;
      sigma[1] = 7.90069e-05;
      alpha[0] = 4.95058e-08;
      alpha[1] = 4.40036e-08;
      n[0]     = 1.9991;
      n[1]     = 1.66272;
      norm[0]  = 0.851294;
      norm[1]  = 0.884502;
    }

    double muonSF = 1;
    if(pt >= 17)
    {
      double muonDataEff = efficiency(pt, m0[0], sigma[0], alpha[0], n[0], norm[0]);
      double muonMCEff   = efficiency(pt, m0[1], sigma[1], alpha[1], n[1], norm[1]);
      muonSF = muonDataEff/muonMCEff;
    }

    // Tau leg
    eta = tau.eta();
    pt  = tau.pt();
    if(abs(eta) < 1.5) // In barrel
    {
      m0[0]    = 18.604910;
      m0[1]    = 18.532997;
      sigma[0] = 0.276042;
      sigma[1] = 1.027880;
      alpha[0] = 0.137039;
      alpha[1] = 2.262950;
      n[0]     = 2.698437;
      n[1]     = 1.003322;
      norm[0]  = 0.940721;
      norm[1]  = 5.297292;
    }
    else // In endcap
    {
      m0[0]    = 18.701715;
      m0[1]    = 18.212782;
      sigma[0] = 0.216523;
      sigma[1] = 0.338119;
      alpha[0] = 0.148111;
      alpha[1] = 0.122828;
      n[0]     = 2.245081;
      n[1]     = 12.577926;
      norm[0]  = 0.895320;
      norm[1]  = 0.893975;
    }

    double tauSF = 1;
    if(pt >= 20)
    {
      double tauDataEff = efficiency(pt, m0[0], sigma[0], alpha[0], n[0], norm[0]);
      double tauMCEff   = efficiency(pt, m0[1], sigma[1], alpha[1], n[1], norm[1]);
      tauSF = tauDataEff/tauMCEff;
    }

    scaleFactor = muonSF * tauSF;
  }

  return scaleFactor;
}

// Following function from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2012#ETau_MuTau_trigger_turn_on_Joshu
// it parametrizes a trigger efficiency turn on curve. m is the pT of the object
double efficiency(double m, double m0, double sigma, double alpha, double n, double norm)
{
  const double sqrtPiOver2 = 1.2533141373;
  const double sqrt2 = 1.4142135624;
  double sig = fabs((double) sigma);
  double t = (m - m0)/sig;
  if(alpha < 0)
    t = -t;
  double absAlpha = fabs(alpha/sig);
  double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
  double b = absAlpha - n/absAlpha;
  double ApproxErf;
  double arg = absAlpha / sqrt2;
  if (arg > 5.) ApproxErf = 1;
  else if (arg < -5.) ApproxErf = -1;
  else ApproxErf = TMath::Erf(arg);
  double leftArea = (1 + ApproxErf) * sqrtPiOver2;
  double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
  double area = leftArea + rightArea;
  if( t <= absAlpha )
  {
    arg = t / sqrt2;
    if(arg > 5.) ApproxErf = 1;
    else if (arg < -5.) ApproxErf = -1;
    else ApproxErf = TMath::Erf(arg);
    return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
  }
  else
  {
    return norm * (leftArea + a * (1/TMath::Power(t-b,n-1) - 1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area;
  }
}
