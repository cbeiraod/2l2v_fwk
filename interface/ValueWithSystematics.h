// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2015-10-27</date>
// <summary>Header file for the ValueWithSystematics class</summary>
//
// <description>
//  Header file with the declarations of the ValueWithSystematics class and derivates.
//  This class is made to function exactly like the builtin types (and other types when the operators are defined),
// except the types now have associated named systematic uncertainties, which are treated idependently. The systematic
// uncertainties are correctly handled when performing computations and many other tasks, only requiring to be
// handled by hand in specific circumstances.
//  Several methods have been implemented to allow to handle the uncertainties by hand when needed. Also, some other
// methods have been implemented to allow functionality that is normally available to specific types or builtin functions.
// </description>

#ifndef VALUE_WITH_SYSTEMATICS_INCLUDED
#define VALUE_WITH_SYSTEMATICS_INCLUDED

#undef WITH_UNLOCK

#include <exception>

class ValueException: public exception
{
public:
  ValueException(std::string mess): message(mess)
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

// Internal class (do not use)
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


// External classes (to be used)
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

// forward declaration of double specialization
template<>
class ValueWithSystematics<double>;

// vector specialization
template<class T>
class ValueWithSystematics<std::vector<T>>: public ValueWithSystematicsInternal<std::vector<T>>
{
public:
  ValueWithSystematics(): ValueWithSystematicsInternal<std::vector<T>>(std::vector<T>(0)) {};
  ValueWithSystematics(std::vector<T> val): ValueWithSystematicsInternal<std::vector<T>>(val) {};
  ValueWithSystematics(const ValueWithSystematics<std::vector<T>>& val): ValueWithSystematicsInternal<std::vector<T>>(val) {}; // Copy constructor
  ValueWithSystematics(const ValueWithSystematicsInternal<std::vector<T>>& val): ValueWithSystematicsInternal<std::vector<T>>(val) {}; // Copy constructor
  // This second copy contructor should be removed (idem for all the derived user classes) [actually it should probably be moved to private]

  ValueWithSystematics<int> size() const;

private:
protected:
  using ValueWithSystematicsInternal<std::vector<T>>::systematics;
  using ValueWithSystematicsInternal<std::vector<T>>::value;
};

// specialization for TLorentzVector derivates
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

// specialization for LorentzVectorF and derivates
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

// double specialization
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

#endif