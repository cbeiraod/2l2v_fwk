// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2015-10-27</date>
// <summary>Implementation file for the Analyser class</summary>
//
// <description>
//  Implementation file with the definitions of the llvvAnalyser class and related classes.
//  The base analyser from which other analyses can be derived from is defined within.
// </description>

#include "UserCode/llvv_fwk/interface/llvvAnalyser.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

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

ValueWithSystematics<double>& EventInfo::GetDouble(std::string name)
{
  if(eventDoubles.count(name) == 0)
    throw AnalyserException("Tried to access non-existing value: "+name);
  return eventDoubles.at(name);
}

ValueWithSystematics<int>&    EventInfo::GetInt   (std::string name)
{
  if(eventInts.count(name) == 0)
    throw AnalyserException("Tried to access non-existing value: "+name);
  return eventInts.at(name);
}

ValueWithSystematics<bool>&   EventInfo::GetBool  (std::string name)
{
  if(eventBools.count(name) == 0)
    throw AnalyserException("Tried to access non-existing value: "+name);
  return eventBools.at(name);
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
  mctruthmode = 0;

  if(cfgOptions.exists("applyScaleFactors"))
    applyScaleFactors = cfgOptions.getParameter<bool>("applyScaleFactors");
  if(cfgOptions.exists("debug"))
    debug             = cfgOptions.getParameter<bool>("debug");
  if(cfgOptions.exists("doDDBkg"))
    doDDBkg           = cfgOptions.getParameter<bool>("doDDBkg");
  if(cfgOptions.exists("outputEventList"))
    outputEventList   = cfgOptions.getParameter<bool>("outputEventList");
  if(cfgOptions.exists("mctruthmode"))
    mctruthmode = cfgOptions.getParameter<int>("mctruthmode");

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
          //crossSection_("topMass_UP")   += 0.52;
          //crossSection_("topMass_DOWN") -= 0.52;
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
          //crossSection_("topMass_UP")   += 0.30;
          //crossSection_("topMass_DOWN") -= 0.30;
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
  if(cfgOptions.exists("mctruthmode"))
  {
    char tmp = '0' + mctruthmode;
    outFile = outDir + "/" + outFileUrl + "_filt" + tmp + ".root";
  }
  else
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
    if(iev < static_cast<size_t>(skipEvents))
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
      //jesCor->setNPV(nvtx); ?

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
    //metHandle.getByLabel(ev, "llvvObjectProducersUsed", "pfType1CorrectedMet");
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
    if(outputEventList && (!eventlistSelected || eventContent.GetBool("selected").Value()))
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
  
  eventContent.AddInt("nvtx", -1);
  
  eventContent.AddDouble("rho", -1);
  eventContent.AddDouble("rho25", -1);
  
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

  //eventContent.Lock(); // It should only be locked after the first iteration through the event loop
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
  if(isMC && runSystematics)
  {
    ValueWithSystematics<double> lumiUnc(1.0);
    lumiUnc("lumi_UP")   *= 1.026;
    lumiUnc("lumi_DOWN") *= 0.974;

    eventContent.GetDouble("weight") *= lumiUnc;
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



