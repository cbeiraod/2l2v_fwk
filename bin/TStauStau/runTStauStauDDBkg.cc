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



#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
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
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TNtuple.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRotation.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/shared_ptr.hpp>
#include <Math/VectorUtil.h>
#include <bitset>
#include <cctype>
#include <cmath>


#ifndef DEBUG_EVENT
//#define DEBUG_EVENT true
#endif

#define NAN_WARN(X) if(std::isnan(X)) std::cout << "  Warning: " << #X << " is nan" << std::endl;

enum ID_Type {LooseID, MediumID, TightID};
enum TAU_E_ID {antiELoose, antiEMedium, antiETight, antiEMva, antiEMva3Loose, antiEMva3Medium, antiEMva3Tight, antiEMva3VTight, antiEMva5Medium};

bool electronMVAID(double mva, llvvLepton& lepton, ID_Type id);

/*****************************************************************************/
/* Return Codes:                                                             */
/*   0 - Everything OK                                                       */
/*   1 - Missing parameters_cfg.py configuration file                        */
/*****************************************************************************/
int main(int argc, char* argv[])
{
  if(argc < 2)
    std::cout << "Usage: " << argv[0] << " parameters_cfg.py" << std::endl, exit(1);

  size_t limit = 0;

  #if defined(DEBUG_EVENT)
  bool debugEvent = false;
  int skipEvents = 0;
  #endif

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

      #if defined(DEBUG_EVENT)
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
      #endif
    }
  }

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  // Read parameters from the configuration file
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[fileIndex])->getParameter<edm::ParameterSet>("runProcess");

  bool isMC = runProcess.getParameter<bool>("isMC");
  double xsec = runProcess.getParameter<double>("xsec");
  std::vector<std::string> urls = runProcess.getParameter<std::vector<std::string> >("input");
  std::string baseDir = runProcess.getParameter<std::string>("dirName");
  std::string outdir = runProcess.getParameter<std::string>("outdir");
  std::string jecDir = runProcess.getParameter<std::string>("jecDir");
  bool saveSummaryTree = runProcess.getParameter<bool>("saveSummaryTree");
  bool doPrompt = false;
  if(runProcess.exists("doPrompt"))
    doPrompt = runProcess.getParameter<bool>("doPrompt");
  bool debug = false;
  if(runProcess.exists("debug"))
    debug = runProcess.getParameter<bool>("debug");

  if(debug)
    std::cout << "Finished loading config file" << std::endl;

  // Hardcoded Values
  double sqrtS          =  8;      // Center of mass energy
  double minElPt        = 30;      // Selected electron pT and eta
  double maxElEta       =  2.1;
  double ECALGap_MinEta =  1.4442; // ECAL gap parameters
  double ECALGap_MaxEta =  1.5660;
  double minMuPt        = 27;      // Selected muon pT and eta
  double maxMuEta       =  2.1;
  double minTauPt       = 20;      // Selected tau pT and eta (I was using 25)
  double maxTauEta      =  2.3;
  double maxJetEta      =  4.7;    // Selected jet eta

  // Setting up -------------------------------------------------------------------------
  if(debug)
    std::cout << "Setting up" << std::endl;
  gSystem->Exec(("mkdir -p " + outdir).c_str());
  std::string url = urls[0];
  std::string outFileUrl(gSystem->BaseName(url.c_str()));
  while(outFileUrl.find(".root", 0) != std::string::npos)
    outFileUrl.replace(outFileUrl.find(".root", 0), 5, "");
  std::string outUrl = outdir;
  outUrl += "/";
  outUrl += outFileUrl + ".root";

  TString turl(url);
  bool isV0JetsMC(isMC && (turl.Contains("DYJetsToLL_50toInf") || turl.Contains("WJets")));
  bool isStauStau(isMC && turl.Contains("TStauStau"));
  bool isSingleMuPD(!isMC && (turl.Contains("Data8TeV_SingleMu2012")));

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



  /***************************************************************************/
  /*                         Initializing Histograms                         */
  /***************************************************************************/
  if(debug)
    std::cout << "Initializing histograms" << std::endl;
  SmartSelectionMonitor mon;
  TH1D *eventflow = (TH1D*)mon.addHistogram(new TH1D("eventflow", ";;Events", 8, 0, 8));
  eventflow->GetXaxis()->SetBinLabel(1, "HLT");
  eventflow->GetXaxis()->SetBinLabel(2, "> 1l");
  eventflow->GetXaxis()->SetBinLabel(3, "met > 30");
  eventflow->GetXaxis()->SetBinLabel(4, "Loose #tau");
  eventflow->GetXaxis()->SetBinLabel(5, "");
  eventflow->GetXaxis()->SetBinLabel(6, "");
  eventflow->GetXaxis()->SetBinLabel(7, "");
  eventflow->GetXaxis()->SetBinLabel(8, "");

  TH1D *genParticleStatus = static_cast<TH1D*>(mon.addHistogram(new TH1D("genStatus", ";genParticleStatus;Events", 8, 0, 8)));
  genParticleStatus->GetXaxis()->SetBinLabel(1, "pp");
  genParticleStatus->GetXaxis()->SetBinLabel(2, "pf");
  genParticleStatus->GetXaxis()->SetBinLabel(3, "fp");
  genParticleStatus->GetXaxis()->SetBinLabel(4, "ff");
  genParticleStatus->GetXaxis()->SetBinLabel(5, "");
  genParticleStatus->GetXaxis()->SetBinLabel(6, "data");
  genParticleStatus->GetXaxis()->SetBinLabel(7, "");
  genParticleStatus->GetXaxis()->SetBinLabel(8, "err");

  TH1D *genParticleStatusTight = static_cast<TH1D*>(mon.addHistogram(new TH1D("genStatusTight", ";genParticleStatus;Events", 8, 0, 8)));
  genParticleStatusTight->GetXaxis()->SetBinLabel(1, "pp");
  genParticleStatusTight->GetXaxis()->SetBinLabel(2, "pf");
  genParticleStatusTight->GetXaxis()->SetBinLabel(3, "fp");
  genParticleStatusTight->GetXaxis()->SetBinLabel(4, "ff");
  genParticleStatusTight->GetXaxis()->SetBinLabel(5, "");
  genParticleStatusTight->GetXaxis()->SetBinLabel(6, "data");
  genParticleStatusTight->GetXaxis()->SetBinLabel(7, "");
  genParticleStatusTight->GetXaxis()->SetBinLabel(8, "err");

  TH1D *genTauStatus = static_cast<TH1D*>(mon.addHistogram(new TH1D("genTauStatus", ";genTauStatus;Taus", 6, 0, 6)));
  genTauStatus->GetXaxis()->SetBinLabel(1, "p");
  genTauStatus->GetXaxis()->SetBinLabel(2, "f");
  genTauStatus->GetXaxis()->SetBinLabel(3, "");
  genTauStatus->GetXaxis()->SetBinLabel(4, "data");
  genTauStatus->GetXaxis()->SetBinLabel(5, "");
  genTauStatus->GetXaxis()->SetBinLabel(6, "err");

  TH1D *genTauStatusTight = static_cast<TH1D*>(mon.addHistogram(new TH1D("genTauStatusTight", ";genTauStatus;Taus", 6, 0, 6)));
  genTauStatusTight->GetXaxis()->SetBinLabel(1, "p");
  genTauStatusTight->GetXaxis()->SetBinLabel(2, "f");
  genTauStatusTight->GetXaxis()->SetBinLabel(3, "");
  genTauStatusTight->GetXaxis()->SetBinLabel(4, "data");
  genTauStatusTight->GetXaxis()->SetBinLabel(5, "");
  genTauStatusTight->GetXaxis()->SetBinLabel(6, "err");

  mon.addHistogram(new TH1D("nup", ";NUP;Events", 10, 0, 10));

  // Pile Up
  //mon.addHistogram(new TH1D("nvtxAll", ";Vertices;Events", 50, -0.5, 49.5));
  mon.addHistogram(new TH1D("nvtx", ";Vertices;Events", 50, -0.5, 49.5));
  mon.addHistogram(new TH1D("nvtxraw", ";Vertices;Events", 50, -0.5, 49.5));
  mon.addHistogram(new TH1D("rho", ";#rho;Events", 25, 0, 25));
  mon.addHistogram(new TH1D("rho25", ";#rho(#eta<2.5);Events", 25, 0, 25));


  // Leptons
  TH1D *leptonCutFlow = (TH1D*)mon.addHistogram(new TH1D("leptonCutFlow", ";;Leptons", 4, 0, 4));
  leptonCutFlow->GetXaxis()->SetBinLabel(1, "All");
  leptonCutFlow->GetXaxis()->SetBinLabel(2, "ID");
  leptonCutFlow->GetXaxis()->SetBinLabel(3, "Kin");
  leptonCutFlow->GetXaxis()->SetBinLabel(4, "Iso");


  // Taus
  mon.addHistogram(new TH1D("ptSelectedTau", ";p_{T}^{#tau};Events", 20, 0, 100));
  mon.addHistogram(new TH1D("ptSelectedTauExtended", ";p_{T}^{#tau};Events", 50, 0, 250));
  mon.addHistogram(new TH1D("ptSelectedTauTight", ";p_{T}^{#tau};Events", 20, 0, 100));
  mon.addHistogram(new TH1D("ptSelectedTauExtendedTight", ";p_{T}^{#tau};Events", 50, 0, 250));
  mon.addHistogram(new TH1D("etaSelectedTau", ";#eta^{#tau};Events", 26, -2.6, 2.6));
  mon.addHistogram(new TH1D("etaSelectedTauTight", ";#eta^{#tau};Events", 26, -2.6, 2.6));
  mon.addHistogram(new TH2D("ptetaSelectedTau", ";p_{T}^{#tau};#eta^{#tau}", 20, 0, 100, 26, -2.6, 2.6))->SetOption("colz");
  mon.addHistogram(new TH2D("ptetaSelectedTauTight", ";p_{T}^{#tau};#eta^{#tau}", 20, 0, 100, 26, -2.6, 2.6))->SetOption("colz");
  mon.addHistogram(new TH2D("ptetaSelectedTauExtended", ";p_{T}^{#tau};#eta^{#tau}", 100, 0, 500, 26, -2.6, 2.6))->SetOption("colz");
  mon.addHistogram(new TH2D("ptetaSelectedTauExtendedTight", ";p_{T}^{#tau};#eta^{#tau}", 100, 0, 500, 26, -2.6, 2.6))->SetOption("colz");
  mon.addHistogram(new TH2D("ptabsetaSelectedTau", ";p_{T}^{#tau};|#eta^{#tau}|", 20, 0, 100, 13, 0, 2.6))->SetOption("colz");
  mon.addHistogram(new TH2D("ptabsetaSelectedTauTight", ";p_{T}^{#tau};|#eta^{#tau}|", 20, 0, 100, 13, 0, 2.6))->SetOption("colz");
  mon.addHistogram(new TH2D("ptabsetaSelectedTauExtended", ";p_{T}^{#tau};|#eta^{#tau}|", 100, 0, 500, 13, 0, 2.6))->SetOption("colz");
  mon.addHistogram(new TH2D("ptabsetaSelectedTauExtendedTight", ";p_{T}^{#tau};|#eta^{#tau}|", 100, 0, 500, 13, 0, 2.6))->SetOption("colz");
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

  // MET
  mon.addHistogram(new TH1D("MET", ";MET [GeV];Events", 25, 0, 250));
  mon.addHistogram(new TH1D("METTight", ";MET [GeV];Events", 25, 0, 250));
  mon.addHistogram(new TH1D("TauMET", ";MET [GeV];Events", 25, 0, 250));
  mon.addHistogram(new TH1D("TauMETTight", ";MET [GeV];Events", 25, 0, 250));



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

  if(debug)
    std::cout << "Getting PU weights" << std::endl;
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
  if(debug)
    std::cout << "Done with PU weights" << std::endl;

  gROOT->cd(); //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

  DuplicatesChecker duplicatesChecker;
//  int nDuplicates(0);
  int step = int(totalEntries/50);

  // Redirect stdout and stderr to a temporary buffer, then output buffer after event loop
  std::ostream myCout(std::cout.rdbuf());
  std::stringstream buffer;
  std::streambuf *coutbuf = std::cout.rdbuf();
  std::streambuf *cerrbuf = std::cerr.rdbuf();
  std::cout.rdbuf(buffer.rdbuf());
  std::cerr.rdbuf(buffer.rdbuf());

  // Variables used in loop
  if(debug)
    myCout << "  Declaring all variables used in loop" << std::endl;
  int nvtx = 0;
  bool selected = false;
  bool istight  = false;
  int genStatus = 7;
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
  llvvLeptonCollection selLeptons;
  llvvTauCollection selTaus;
  int tauIndex = -1, leptonIndex = -1;
  bool isOS = false;

  // Prepare summary tree
  if(saveSummaryTree)
  {
    if(debug)
      std::cout << "  Defining all branches in output root file" << std::endl;

    TDirectory* cwd = gDirectory;
    summaryOutFile->cd();

    // Dataset specific variables
    summaryTree->Branch("isMC", &isMC);
    summaryTree->Branch("crossSection", &crossSection);
    // Add here nInitEvent
//    summaryTree->Branch("xSecWeight", &xsecWeight);

    // Event specific variables
    summaryTree->Branch("selected", &selected);
    summaryTree->Branch("istight",  &istight);
    summaryTree->Branch("genStatus",&genStatus);
    summaryTree->Branch("nvtx", &nvtx);
    summaryTree->Branch("triggeredOn", &triggeredOn);
    summaryTree->Branch("rho", &rho);
    summaryTree->Branch("rho25", &rho25);
    summaryTree->Branch("met", &met);
    summaryTree->Branch("weight", &weight);
    summaryTree->Branch("weight_plus", &weight_plus);
    summaryTree->Branch("weight_minus", &weight_minus);
    summaryTree->Branch("puWeight", &puWeight);
    summaryTree->Branch("isOS", &isOS);

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
    istight  = false;
    genStatus = 7;
    nvtx = 0;
    weight = 1.;
    weight_plus = 1.;
    weight_minus = 1.;
    puWeight = 1.;
    chTags.clear();
    selLeptons.clear();
    selTaus.clear();
    tauIndex = -1, leptonIndex = -1;
    isOS = false;

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
    if(isV0JetsMC)
    {
      if(genEv.nup > 5) // Drop V+{1,2,3,4}Jets from VJets samples to avoid double counting (but keep V+0Jets) [V = W,Z]
        continue;
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

    bool singleETrigger = triggerBits[13]; // HLT_Ele27_WP80_v*
    bool singleMuTrigger = triggerBits[15]; // HLT_IsoMu24_v*

    triggeredOn = singleETrigger || singleMuTrigger;
    if(isSingleMuPD) // Remove repeated events from different Primary Dataset
    {
      if(singleETrigger)
        triggeredOn = false;
    }

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
      puWeight     = LumiWeights->weight(genEv.ngenITpu) * PUNorm[0];
      weight       = xsecWeight*puWeight;
      weight_plus  = PuShifters[utils::cmssw::PUUP ]->Eval(genEv.ngenITpu) * (PUNorm[2]/PUNorm[0]);
      weight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(genEv.ngenITpu) * (PUNorm[1]/PUNorm[0]);
    }


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
      int lepId = leptons[i].id;

      if(abs(lepId) == 13 && muCor)
      {
        TLorentzVector p4(leptons[i].px(), leptons[i].py(), leptons[i].pz(), leptons[i].energy());
        muCor->applyPtCorrection(p4, (lepId>0)?1:-1);
        if(isMC)
          muCor->applyPtSmearing(p4, (lepId>0)?1:-1, false);
        leptons[i].SetPxPyPzE(p4.Px(), p4.Py(), p4.Pz(), p4.Energy());
      }

      lepId = abs(lepId);

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
        if(deltaR(tau, selLeptons[lep]) < 0.1)
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
//      if(!doTightTauID)
//      {
//        if(!tau.passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits)) passID = false;
//      }
//      else
//      {
//        if(!tau.passId(llvvTAUID::byTightCombinedIsolationDeltaBetaCorr3Hits)) passID = false;
//      }
//      if(!tau.passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits)) passID = false;
      if(!tau.passId(llvvTAUID::againstMuonTight3)) passID = false;
      if(!tau.passId(llvvTAUID::againstElectronMediumMVA5)) passID = false;

      if(passID && passKin && passIso && passQual && tau.isPF)
        selTaus.push_back(tau);
      if(!triggeredOn)
        continue;

      // Fill control histograms
      mon.fillHisto("tauCutFlow", chTags, 0, weight);
      if(tau.isPF)
      {
        mon.fillHisto("tauCutFlow", chTags, 1, weight);
        mon.fillHisto("tauID", chTags, 0, weight);
        if(tau.passId(llvvTAUID::byTightCombinedIsolationDeltaBetaCorr3Hits))
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

    if(selLeptons.size() != 0)
      std::sort(selLeptons.begin(), selLeptons.end(), sort_llvvObjectByPt);

    if(selTaus.size() != 0)
      std::sort(selTaus.begin(), selTaus.end(), sort_llvvObjectByPt);



    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Requiring an opposite sign pair" << std::endl;
    #endif
    // Opposite Sign requirements
    double maxPtSum = 0;
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
          }
        }
        if(PtSum < maxPtSum) // Skip a few iterations if it is not expected that we will find a better candidate
          break;
      }
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Filling histograms" << std::endl;
    #endif

    if(triggeredOn && selLeptons.size() > 0 && selTaus.size() > 0 && ((!doPrompt && met.pt() > 50) || doPrompt))
    {
      selected = true;
      chTags.push_back("selected");
    }

    if(triggeredOn)
    {
      chTags.push_back("HLT");
      mon.fillHisto("eventflow", chTags, 0, weight);
      if(selLeptons.size() > 0)
      {
        chTags.push_back("1lepton");
        mon.fillHisto("eventflow", chTags, 1, weight);
        if(doPrompt || (!doPrompt && met.pt() > 50))
        {
          chTags.push_back("met");
          mon.fillHisto("eventflow", chTags, 2, weight);
          if(selTaus.size() > 0)
          {
            chTags.push_back("1tau");
            mon.fillHisto("eventflow", chTags, 3, weight);
            bool hasTight = false;

            if(abs(selLeptons[0].id) == 11)
              chTags.push_back("leadingE");
            else
              chTags.push_back("leadingMu");
            for(auto &tau : selTaus)
            {
              bool isTight = tau.passId(llvvTAUID::byTightCombinedIsolationDeltaBetaCorr3Hits);

              bool isPrompt = false;
              int status = 3;
              if(isMC)
              {
                status = 1;
                for(auto &genPart : gen)
                {
                  if(genPart.status != 3)
                    continue;

                  if(genPart.id == tau.id)
                  {
                    if(deltaR(tau, genPart) < 0.3)
                    {
                      isPrompt = true;
                      status = 0;
                    }
                  }
                }
              }

              mon.fillHisto("genTauStatus", chTags, status, weight);
              mon.fillHisto("ptSelectedTau", chTags, tau.pt(), weight);
              mon.fillHisto("ptSelectedTauExtended", chTags, tau.pt(), weight);
              mon.fillHisto("etaSelectedTau", chTags, tau.eta(), weight);
              mon.fillHisto("TauMET", chTags, met.pt(), weight);
              mon.fillHisto("ptetaSelectedTau", chTags, tau.pt(), tau.eta(), weight);
              mon.fillHisto("ptetaSelectedTauExtended", chTags, tau.pt(), tau.eta(), weight);
              mon.fillHisto("ptabsetaSelectedTau", chTags, tau.pt(), abs(tau.eta()), weight);
              mon.fillHisto("ptabsetaSelectedTauExtended", chTags, tau.pt(), abs(tau.eta()), weight);

              if(isTight)
              {
                hasTight = true;
                mon.fillHisto("genTauStatusTight", chTags, status, weight);
                mon.fillHisto("ptSelectedTauTight", chTags, tau.pt(), weight);
                mon.fillHisto("ptSelectedTauExtendedTight", chTags, tau.pt(), weight);
                mon.fillHisto("etaSelectedTauTight", chTags, tau.eta(), weight);
                mon.fillHisto("TauMETTight", chTags, met.pt(), weight);
                mon.fillHisto("ptetaSelectedTauTight", chTags, tau.pt(), tau.eta(), weight);
                mon.fillHisto("ptetaSelectedTauExtendedTight", chTags, tau.pt(), tau.eta(), weight);
                mon.fillHisto("ptabsetaSelectedTauTight", chTags, tau.pt(), abs(tau.eta()), weight);
                mon.fillHisto("ptabsetaSelectedTauExtendedTight", chTags, tau.pt(), abs(tau.eta()), weight);
              }
            }

            mon.fillHisto("MET", chTags, met.pt(), weight);
            if(hasTight)
            {
              mon.fillHisto("METTight", chTags, met.pt(), weight);
              istight = true;
            }

          }
        }
      }
    }

    mon.fillHisto("genStatus", chTags, genStatus, weight);
    if(istight)
    {
        mon.fillHisto("genStatusTight", chTags, genStatus, weight);
    }

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
    NAN_WARN(minDeltaPhiMetJet40)
    NAN_WARN(tauLeadPt)
    NAN_WARN(lepLeadPt)
    NAN_WARN(maxPtSum)
    #endif
    if(saveSummaryTree)
    {
      TDirectory* cwd = gDirectory;
      summaryOutFile->cd();
      summaryTree->Fill();
      cwd->cd();
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
  std::cout << buffer.str();

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
//    std::string summaryOutUrl = outUrl;
//    summaryOutUrl.replace(summaryOutUrl.find(".root", 0), 5, "_summary.root");
//    std::cout << "Saving summary results in " << summaryOutUrl << std::endl;
//    summaryOutFile = new TFile(summaryOutUrl.c_str(), "RECREATE");
    // Write Tuple/Tree
//    summaryTree->SetDirectory(summaryOutFile);
    TDirectory* cwd = gDirectory;
    summaryOutFile->cd();
    summaryTree->Write();
    summaryOutFile->Close();
    delete summaryOutFile;
    cwd->cd();
  }

  return 0;
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
