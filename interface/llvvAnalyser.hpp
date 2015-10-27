// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2015-10-27</date>
// <summary>Implementation file for the Analyser class</summary>
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

#ifndef LLVV_ANALYSER_H
#define LLVV_ANALYSER_H

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

#endif


