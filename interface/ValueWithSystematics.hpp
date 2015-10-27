// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2015-10-27</date>
// <summary>Implementation file for the ValueWithSystematics class</summary>
//
// <description>
//  Implementation file with the definitions of the ValueWithSystematics class and derivates.
//  This class is made to function exactly like the builtin types (and other types when the operators are defined),
// except the types now have associated named systematic uncertainties, which are treated idependently. The systematic
// uncertainties are correctly handled when performing computations and many other tasks, only requiring to be
// handled by hand in specific circumstances.
//  Several methods have been implemented to allow to handle the uncertainties by hand when needed. Also, some other
// methods have been implemented to allow functionality that is normally available to specific types or builtin functions.
// </description>

#ifndef VALUE_WITH_SYSTEMATICS_HPP
#define VALUE_WITH_SYSTEMATICS_HPP

//#include "UserCode/llvv_fwk/interface/ValueWithSystematics.h"

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
      throw ValueException("Unable to add metadata \""+key+":"+value+"\" after locking.");

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
      throw ValueException("Unable to add systematic \""+name+"\" after locking.");
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
      throw ValueException("Unable to add systematic \""+name+"\" after locking.");
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


// Specializations --------------------------------------------------------------------------------
//****************** std::vector *******************************************************************
template<class T>
ValueWithSystematics<int> ValueWithSystematics<std::vector<T>>::size() const
{
  ValueWithSystematics<int> retVal(value.size());

  for(auto& kv: systematics)
    retVal(kv.first) = kv.second.size();

  return retVal;
}

template<class T>
ValueWithSystematics<T> ValueWithSystematics<std::vector<T>>::first() const
{
  ValueWithSystematics<T> retVal;

  for(auto& kv: systematics)
  {
    retVal(kv.first);
    if(kv.second.size() > 0)
      retVal(kv.first) = kv.second.at(0);
  }

  if(value.size() > 0)
    retVal.value = value.at(0);

  return retVal;
}

//****************** TLorentzVector ****************************************************************
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

//****************** LorentzVectorF ****************************************************************
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

//****************** double ************************************************************************
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

#endif

