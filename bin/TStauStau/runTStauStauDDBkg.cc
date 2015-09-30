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
double tauSF(llvvTau& tau, llvvGenParticleCollection& genPartColl, TAU_E_ID eId);
double leptonIdAndIsoScaleFactor(llvvLepton& lepton);
double leptonTauTriggerScaleFactor(llvvLepton& lepton, llvvTau& tau);
double efficiency(double m, double m0, double sigma, double alpha, double n, double norm);

double getFromTH2(TH2* hist, double x, double y);

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
  bool debug = false;
  if(runProcess.exists("debug"))
    debug = runProcess.getParameter<bool>("debug");

  if(debug)
    std::cout << "Finished loading config file" << std::endl;

  // Hardcoded Values
  double sqrtS          =  8;      // Center of mass energy
//  double minElPt        = 30;      // Selected electron pT and eta
  double minElPt        = 25;      // Selected electron pT and eta
  double maxElEta       =  2.1;
  double ECALGap_MinEta =  1.4442; // ECAL gap parameters
  double ECALGap_MaxEta =  1.5660;
//  double minMuPt        = 27;      // Selected muon pT and eta
  double minMuPt        = 20;      // Selected muon pT and eta
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

  TH2D *etauFR = NULL;
  TH2D *mutauFR = NULL;
  TH2D *etauPR = NULL;
  TH2D *mutauPR = NULL;
  TDirectory* cwd = gDirectory;

  std::string FRFileName = gSystem->ExpandPathName("$CMSSW_BASE/src/UserCode/llvv_fwk/data/TStauStau/rates.root");
  std::string PRFileName = gSystem->ExpandPathName("$CMSSW_BASE/src/UserCode/llvv_fwk/data/TStauStau/rates.root");
  std::cout << "Trying to open FR file: " << FRFileName << std::endl;
  TFile FRFile(FRFileName.c_str(), "READ");
  std::cout << "Trying to open PR file: " << PRFileName << std::endl;
  TFile PRFile(PRFileName.c_str(), "READ");
  cwd->cd();

//  etauFR  = static_cast<TH2D*>(FRFile.Get("data-Z/data-Z_leadingE_varptetaSelectedTau_FR")->Clone("etauFR"));
//  mutauFR = static_cast<TH2D*>(FRFile.Get("data-Z/data-Z_leadingMu_varptetaSelectedTau_FR")->Clone("mutauFR"));
  etauFR  = static_cast<TH2D*>(FRFile.Get("data-Zprompt/data-Zprompt_InvMET_OS_etaSelectedTau_FR")->Clone("etauFR"));
  mutauFR = static_cast<TH2D*>(FRFile.Get("data-Zprompt/data-Zprompt_InvMET_OS_etaSelectedTau_FR")->Clone("mutauFR"));

  etauPR  = static_cast<TH2D*>(PRFile.Get("Z #rightarrow ll/Zrightarrowll_InvMET_OS_etaSelectedTau_FR")->Clone("etauPR"));
  mutauPR = static_cast<TH2D*>(PRFile.Get("Z #rightarrow ll/Zrightarrowll_InvMET_OS_etaSelectedTau_FR")->Clone("mutauPR"));

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

  TH1D *genTauStatus = static_cast<TH1D*>(mon.addHistogram(new TH1D("genTauStatus", ";genTauStatus;Taus", 4, 0, 4)));
  genTauStatus->GetXaxis()->SetBinLabel(1, "prompt");
  genTauStatus->GetXaxis()->SetBinLabel(2, "fake");
  genTauStatus->GetXaxis()->SetBinLabel(3, "");
  genTauStatus->GetXaxis()->SetBinLabel(4, "data");
//  genTauStatus->GetXaxis()->SetBinLabel(5, "");
//  genTauStatus->GetXaxis()->SetBinLabel(6, "err");

  TH1D *genTauStatusTight = static_cast<TH1D*>(mon.addHistogram(new TH1D("genTauStatusTight", ";genTauStatus;Taus", 4, 0, 4)));
  genTauStatusTight->GetXaxis()->SetBinLabel(1, "prompt");
  genTauStatusTight->GetXaxis()->SetBinLabel(2, "fake");
  genTauStatusTight->GetXaxis()->SetBinLabel(3, "");
  genTauStatusTight->GetXaxis()->SetBinLabel(4, "data");
//  genTauStatusTight->GetXaxis()->SetBinLabel(5, "");
//  genTauStatusTight->GetXaxis()->SetBinLabel(6, "err");

  TH1D *tauTypes = static_cast<TH1D*>(mon.addHistogram(new TH1D("tauTypes", ";tauType;Taus", 6, 0, 6)));
  tauTypes->GetXaxis()->SetBinLabel(1, "data");
  tauTypes->GetXaxis()->SetBinLabel(2, "prompt");
  tauTypes->GetXaxis()->SetBinLabel(3, "wrong sign");
  tauTypes->GetXaxis()->SetBinLabel(4, "quark");
  tauTypes->GetXaxis()->SetBinLabel(5, "gluon");
  tauTypes->GetXaxis()->SetBinLabel(6, "unknown");

  mon.addHistogram(new TH1D("nup", ";NUP;Events", 10, 0, 10));

  mon.addHistogram(new TH1D("chargedParticleSig", ";#part;Taus", 4, 0, 4));
  mon.addHistogram(new TH1D("chargedParticleIso", ";#part;Taus", 10, 0, 10));

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

  Double_t etabinning[] = {-2.3, -2.1, -1.9, -1.7, -1.4, -1.1, -0.7, 0, 0.7, 1.1, 1.4, 1.7, 1.9, 2.1, 2.3};
  int etabins = sizeof(etabinning)/sizeof(Double_t) - 1;
  Double_t absetabinning[] = {0, 0.7, 1.1, 1.4, 1.7, 1.9, 2.1, 2.3};
  int absetabins = sizeof(absetabinning)/sizeof(Double_t) - 1;
  Double_t ptbinning[] = {20, 22.5, 26, 30, 35, 40, 45, 55, 70, 85, 100};
  int ptbins = sizeof(ptbinning)/sizeof(Double_t) - 1;
  Double_t ptextendedbinning[] = {20, 22.5, 26, 30, 35, 40, 45, 55, 70, 85, 100, 125, 150, 250, 500};
  int ptextendedbins = sizeof(ptextendedbinning)/sizeof(Double_t) - 1;

  // Taus
  mon.addHistogram(new TH1D("ptSelectedTau", ";p_{T}^{#tau};Taus", 16, 20, 100));
  mon.addHistogram(new TH1D("ptSelectedTauExtended", ";p_{T}^{#tau};Taus", 46, 20, 250));
  mon.addHistogram(new TH1D("ptSelectedTauTight", ";p_{T}^{#tau};Taus", 16, 20, 100));
  mon.addHistogram(new TH1D("ptSelectedTauExtendedTight", ";p_{T}^{#tau};Taus", 46, 20, 250));
  mon.addHistogram(new TH1D("etaSelectedTau", ";#eta^{#tau};Taus", 23, -2.3, 2.3));
  mon.addHistogram(new TH1D("etaSelectedTauTight", ";#eta^{#tau};Taus", 23, -2.3, 2.3));
  mon.addHistogram(new TH1D("cosPhiSelectedTau", ";cos#Phi;Taus", 20, -1, 1));
  mon.addHistogram(new TH1D("cosPhiSelectedTauTight", ";cos#Phi;Taus", 20, -1, 1));
  mon.addHistogram(new TH1D("deltaPhi", ";#Delta#Phi_{l-#tau};Taus", 20, -3.14, 3.14));
  mon.addHistogram(new TH1D("absDeltaPhi", ";|#Delta#Phi_{l-#tau}|;Taus", 20, 0, 3.14));
  mon.addHistogram(new TH1D("deltaR", ";#DeltaR_{l-#tau};Taus", 20, 0, 5));
  mon.addHistogram(new TH1D("deltaPhiTight", ";#Delta#Phi_{l-#tau};Taus", 20, -3.14, 3.14));
  mon.addHistogram(new TH1D("absDeltaPhiTight", ";|#Delta#Phi_{l-#tau}|;Taus", 20, 0, 3.14));
  mon.addHistogram(new TH1D("deltaRTight", ";#DeltaR_{l-#tau};Taus", 20, 0, 5));


  mon.addHistogram(new TH2D("ptetaSelectedTau", ";p_{T}^{#tau};#eta^{#tau}", 16, 20, 100, 23, -2.3, 2.3))->SetOption("colz");
  mon.addHistogram(new TH2D("ptetaSelectedTauTight", ";p_{T}^{#tau};#eta^{#tau}", 16, 20, 100, 23, -2.3, 2.3))->SetOption("colz");
  mon.addHistogram(new TH2D("ptetaSelectedTauExtended", ";p_{T}^{#tau};#eta^{#tau}", 96, 20, 500, 23, -2.3, 2.3))->SetOption("colz");
  mon.addHistogram(new TH2D("ptetaSelectedTauExtendedTight", ";p_{T}^{#tau};#eta^{#tau}", 96, 20, 500, 23, -2.3, 2.3))->SetOption("colz");
  mon.addHistogram(new TH2D("ptabsetaSelectedTau", ";p_{T}^{#tau};|#eta^{#tau}|", 16, 20, 100, 11, 0, 2.3))->SetOption("colz");
  mon.addHistogram(new TH2D("ptabsetaSelectedTauTight", ";p_{T}^{#tau};|#eta^{#tau}|", 16, 20, 100, 11, 0, 2.3))->SetOption("colz");
  mon.addHistogram(new TH2D("ptabsetaSelectedTauExtended", ";p_{T}^{#tau};|#eta^{#tau}|", 96, 20, 500, 11, 0, 2.3))->SetOption("colz");
  mon.addHistogram(new TH2D("ptabsetaSelectedTauExtendedTight", ";p_{T}^{#tau};|#eta^{#tau}|", 96, 20, 500, 11, 0, 2.3))->SetOption("colz");

  mon.addHistogram(new TH1D("varptSelectedTau", ";p_{T}^{#tau};Taus", ptbins, ptbinning));
  mon.addHistogram(new TH1D("varptSelectedTauExtended", ";p_{T}^{#tau};Taus", ptextendedbins, ptextendedbinning));
  mon.addHistogram(new TH1D("varptSelectedTauTight", ";p_{T}^{#tau};Taus", ptbins, ptbinning));
  mon.addHistogram(new TH1D("varptSelectedTauExtendedTight", ";p_{T}^{#tau};Taus", ptextendedbins, ptextendedbinning));
  mon.addHistogram(new TH1D("varetaSelectedTau", ";#eta^{#tau};Taus", etabins, etabinning));
  mon.addHistogram(new TH1D("varetaSelectedTauTight", ";#eta^{#tau};Taus", etabins, etabinning));
  mon.addHistogram(new TH2D("varptetaSelectedTau", ";p_{T}^{#tau};#eta^{#tau}", ptbins, ptbinning, etabins, etabinning))->SetOption("colz");
  mon.addHistogram(new TH2D("varptetaSelectedTauTight", ";p_{T}^{#tau};#eta^{#tau}", ptbins, ptbinning, etabins, etabinning))->SetOption("colz");
  mon.addHistogram(new TH2D("varptetaSelectedTauExtended", ";p_{T}^{#tau};#eta^{#tau}", ptextendedbins, ptextendedbinning, etabins, etabinning))->SetOption("colz");
  mon.addHistogram(new TH2D("varptetaSelectedTauExtendedTight", ";p_{T}^{#tau};#eta^{#tau}", ptextendedbins, ptextendedbinning, etabins, etabinning))->SetOption("colz");
  mon.addHistogram(new TH2D("varptabsetaSelectedTau", ";p_{T}^{#tau};|#eta^{#tau}|", ptbins, ptbinning, absetabins, absetabinning))->SetOption("colz");
  mon.addHistogram(new TH2D("varptabsetaSelectedTauTight", ";p_{T}^{#tau};|#eta^{#tau}|", ptbins, ptbinning, absetabins, absetabinning))->SetOption("colz");
  mon.addHistogram(new TH2D("varptabsetaSelectedTauExtended", ";p_{T}^{#tau};|#eta^{#tau}|", ptextendedbins, ptextendedbinning, absetabins, absetabinning))->SetOption("colz");
  mon.addHistogram(new TH2D("varptabsetaSelectedTauExtendedTight", ";p_{T}^{#tau};|#eta^{#tau}|", ptextendedbins, ptextendedbinning, absetabins, absetabinning))->SetOption("colz");
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

  //Leptons
  mon.addHistogram(new TH1D("ptSelectedLep", ";p_{T}^{#tau};Taus", 16, 20, 100));
  mon.addHistogram(new TH1D("ptSelectedLepExtended", ";p_{T}^{#tau};Taus", 46, 20, 250));
  mon.addHistogram(new TH1D("ptSelectedLepTight", ";p_{T}^{#tau};Taus", 16, 20, 100));
  mon.addHistogram(new TH1D("ptSelectedLepExtendedTight", ";p_{T}^{#tau};Taus", 46, 20, 250));
  mon.addHistogram(new TH1D("etaSelectedLep", ";#eta^{#tau};Taus", 23, -2.3, 2.3));
  mon.addHistogram(new TH1D("etaSelectedLepTight", ";#eta^{#tau};Taus", 23, -2.3, 2.3));
  mon.addHistogram(new TH1D("cosPhiSelectedLep", ";cos#Phi;Taus", 20, -1, 1));
  mon.addHistogram(new TH1D("cosPhiSelectedLepTight", ";cos#Phi;Taus", 20, -1, 1));

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
  if(debug)
    std::cout << "  Made chain" << std::endl;
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
  std::vector<TString> ddTags;
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
  double ddWeight = 1.;
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
    summaryTree->Branch("ddWeight", &ddWeight);
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
    ddWeight = 1.;
    chTags.clear();
    ddTags.clear();
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

    bool TauPlusE2012A = triggerBits[18]; // HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*
    bool TauPlusMu2012A = triggerBits[22]; // HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v*
    bool TauPlusE2012B = triggerBits[17]; // HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v*
    bool TauPlusMu2012B = triggerBits[21]; // HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*
    bool TauPlusETrigger = TauPlusE2012A || TauPlusE2012B;
    bool TauPlusMuTrigger = TauPlusMu2012A || TauPlusMu2012B;

    triggeredOn = TauPlusETrigger || TauPlusMuTrigger;
//    triggeredOn = singleETrigger || singleMuTrigger;
//    if(isSingleMuPD) // Remove repeated events from different Primary Dataset
//    {
//      if(singleETrigger)
//        triggeredOn = false;
//    }


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

    if(isMC)
    {
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
//          triggerSF *= leptonTauTriggerScaleFactor(*trigE, *trigTau);
          weight *= leptonTauTriggerScaleFactor(*trigE, *trigTau);
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
//          triggerSF *= leptonTauTriggerScaleFactor(*trigMu, *trigTau);
          weight *= leptonTauTriggerScaleFactor(*trigMu, *trigTau);
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
      if(passKin && passID && passIso)
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
      if(selLeptons.size() > 0)
//        if(deltaR(tau, selLeptons[0]) < 0.2)
        if(deltaR(tau, selLeptons[0]) < 0.5)
          passIso = false;

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
      if(!tau.passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits)) passID = false;
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



/*    #if defined(DEBUG_EVENT)
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
    }// */

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Filling histograms" << std::endl;
    #endif

    if(triggeredOn && selLeptons.size() > 0 && selTaus.size() > 0)
    {
      selected = true;
      chTags.push_back("selected");
      ddTags.push_back("selectedClosure");
    }

    if(triggeredOn)
    {
//      chTags.push_back("HLT");
      mon.fillHisto("eventflow", chTags, 0, weight);
      if(selLeptons.size() > 0)
//      if(selLeptons.size() == 1)
      {
//        chTags.push_back("1lepton");
        mon.fillHisto("eventflow", chTags, 1, weight);
//        if((doPrompt && met.pt() < 50) || (!doPrompt && met.pt() > 50))
        if(true)
        {
//          chTags.push_back("met");
          mon.fillHisto("eventflow", chTags, 2, weight);
          if(selTaus.size() > 0)
          {
//            chTags.push_back("1tau");
            mon.fillHisto("eventflow", chTags, 3, weight);
            bool hasTight = false;
            bool hasMet = false;
            if(met.pt() > 30)
              hasMet = true;

            if(abs(selLeptons[0].id) == 11)
            {
              chTags.push_back("leadingE");
              ddTags.push_back("leadingEClosure");
            }
            else
            {
              chTags.push_back("leadingMu");
              ddTags.push_back("leadingMuClosure");
            }

//            double leptonIdIsoSF = leptonIdAndIsoScaleFactor(selLeptons[0]);
            if(isMC)
              weight *= leptonIdAndIsoScaleFactor(selLeptons[0]);

            for(auto &tau : selTaus)
            {
              double tauScaleFactor = 1.0;
              if(isMC)
                tauScaleFactor = tauSF(tau, gen, antiEMva5Medium);

              std::vector<TString> tauTags = chTags;

              bool isTight = tau.passId(llvvTAUID::byTightCombinedIsolationDeltaBetaCorr3Hits);
              bool isPrompt = false;
              int status = 3;
              int tauType = 0;
              if(isMC)
              {
                status = 1;
                tauType = 5;
                double closest = 7.0;
                int closestID = 0;

                // Gen matching with hard interaction (status 3 particles)
                for(auto &genPart : gen)
                {
                  if(genPart.status != 3)
                    continue;

                  double dist = deltaR(tau, genPart);
                  if(dist < 0.3)
                  {
                    if(dist < closest)
                    {
                      closest = dist;
                      closestID = genPart.id;
                    }
                  }
                }

                if(closestID == tau.id)
                {
                  isPrompt = true;
                  status = 0;
                  tauType = 1;
                }
                else
                {
                  switch(abs(closestID))
                  {
                    case 15: // Tau but wrong sign
                      isPrompt = true;
                      status = 0;
                      tauType = 2;
                      break;
                    case 9: // gluon
                    case 21: // gluon
                      tauType = 3;
                      break;
                    case 1: // quark down
                    case 2: // quark up
                    case 3: // quark strange
                    case 4: // quark charm
                    case 5: // quark bottom
                    case 6: // quark top
                      tauType = 4;
                      break;
                    default:
                      std::cout << "Got an unexpected tau type. PDGID = " << closestID << std::endl;
                      break;
                  }
                }

                // Gen matching didn't seem to give very good results, trying jet matching and then using the jet genflavour
                if(!isPrompt)
                  tauType = 5;
                closest = 7;
                closestID = 0;
                for(auto &jet : jets)
                {
                  double dist = deltaR(tau.jet, jet);
                  if(dist < closest)
                  {
                    closest = dist;
                    closestID = jet.genflav;
                  }
                }

                if(closest < 0.3 && !isPrompt)
                {
                  std::cout << "Found jet, with ID:" << closestID << "; with a deltaR of " << closest << std::endl;

                  switch(abs(closestID))
                  {
//                    case 15: // Tau
//                      isPrompt = true;
//                      status = 0;
//                      tauType = 2;
//                      break;
                    case 9: // gluon
                    case 21: // gluon
                      tauType = 3;
                      break;
                    case 1: // quark down
                    case 2: // quark up
                    case 3: // quark strange
                    case 4: // quark charm
                    case 5: // quark bottom
                    case 6: // quark top
                      tauType = 4;
                      break;
                    default:
                      std::cout << "Got an unexpected jet type. PDGID = " << closestID << std::endl;
                      break;
                  }
                }
              }

//              if(tauType == 3)
//                tauTags.push_back("Gluon");
//              if(tauType == 4)
//                tauTags.push_back("Quark");
              if(!hasMet)
              {
                tauTags.push_back("InvMET");
                if(isPrompt)
                  tauTags.push_back("InvMET_Prompt");
              }
              if(selLeptons[0].id * tau.id < 0)
              {
                tauTags.push_back("OS");
                if(hasMet)
                {
//                  tauTags.push_back("MET_OS");
//                  if(isPrompt)
//                    tauTags.push_back("MET_OS_Prompt");
                }
                else
                {
                  tauTags.push_back("InvMET_OS");
                  if(isPrompt)
                    tauTags.push_back("InvMET_OS_Prompt");
                  else
                    tauTags.push_back("InvMET_OS_Fake");
                }
/*                if(tauType == 3)
                  tauTags.push_back("OS_Gluon");
                if(tauType == 4)
                  tauTags.push_back("OS_Quark");
                if(selLeptons[0].id > 0)
                  tauTags.push_back("lepPlus");
                else
                  tauTags.push_back("lepMinus");

                if(abs(selLeptons[0].id) == 11)
                {
                  tauTags.push_back("OS_leadingE");
                  if(hasMet)
                    tauTags.push_back("MET_OS_leadingE");
                }
                else
                {
                  tauTags.push_back("OS_leadingMu");
                  if(hasMet)
                    tauTags.push_back("MET_OS_leadingMu");
                }// */
              }
              else
              {
//                tauTags.push_back("SS");
//                if(hasMet)
//                {
//                  tauTags.push_back("MET_SS");
//                  if(isPrompt)
//                    tauTags.push_back("MET_SS_Prompt");
//                }
//                else
//                {
//                  tauTags.push_back("InvMET_SS");
//                  if(isPrompt)
//                    tauTags.push_back("InvMET_SS_Prompt");
//                }
//                if(tau.id < 0)
//                  tauTags.push_back("mm");
//                else
//                  tauTags.push_back("pp");
//                if(tauType == 3)
//                  tauTags.push_back("SS_Gluon");
//                if(tauType == 4)
//                  tauTags.push_back("SS_Quark");
/*                if(abs(selLeptons[0].id) == 11)
                {
                  tauTags.push_back("SS_leadingE");
                  if(hasMet)
                    tauTags.push_back("MET_SS_leadingE");
                }
                else
                {
                  tauTags.push_back("SS_leadingMu");
                  if(hasMet)
                    tauTags.push_back("MET_SS_leadingMu");
                }// */
              }

              if(isPrompt)
              {
                tauTags.push_back("Prompt");
                if(selLeptons[0].id * tau.id < 0)
                  tauTags.push_back("OS_Prompt");
//                else
//                  tauTags.push_back("SS_Prompt");
              }
              else
              {
//                tauTags.push_back("Fake");
//                if(selLeptons[0].id * tau.id < 0)
//                  tauTags.push_back("OS_Fake");
//                else
//                  tauTags.push_back("SS_Fake");
              }

              if(tau.numChargedParticlesSigCone == 1)
              {
                tauTags.push_back("1Prong");
                if(selLeptons[0].id * tau.id < 0)
                  tauTags.push_back("OS_1Prong");
                else
                  tauTags.push_back("SS_1Prong");
              }
              else
              {
                tauTags.push_back("3Prong");
                if(selLeptons[0].id * tau.id < 0)
                  tauTags.push_back("OS_3Prong");
                else
                  tauTags.push_back("SS_3Prong");
              }

/*              if(selLeptons[0].id * tau.id < 0) // If opposite sign pair
                continue;
              if(selLeptons[0].id * tau.id < 0 && !doPrompt) // If opposite sign pair
                continue;
              if(selLeptons[0].id * tau.id > 0 && !doPrompt) // If same sign pair
                continue;

              if(doPrompt && !isPrompt)
              {
                continue;
              }// */

//              mon.fillHisto("chargedParticleSig",          tauTags, tau.numChargedParticlesSigCone, weight*tauScaleFactor);
//              mon.fillHisto("chargedParticleIso",          tauTags, tau.numChargedParticlesIsoCone, weight*tauScaleFactor);
//              mon.fillHisto("tauTypes",                    tauTags, tauType, weight*tauScaleFactor);
//              mon.fillHisto("genTauStatus",                tauTags, status, weight*tauScaleFactor);
              mon.fillHisto("ptSelectedTau",               tauTags, tau.pt(), weight*tauScaleFactor);
              mon.fillHisto("ptSelectedTauExtended",       tauTags, tau.pt(), weight*tauScaleFactor);
              mon.fillHisto("etaSelectedTau",              tauTags, tau.eta(), weight*tauScaleFactor);
//              mon.fillHisto("cosPhiSelectedTau",           tauTags, cos(deltaPhi(tau.phi(), met.phi())), weight*tauScaleFactor);
//              mon.fillHisto("TauMET",                      tauTags, met.pt(), weight*tauScaleFactor);
//              mon.fillHisto("ptetaSelectedTau",            tauTags, tau.pt(), tau.eta(), weight*tauScaleFactor);
//              mon.fillHisto("ptetaSelectedTauExtended",    tauTags, tau.pt(), tau.eta(), weight*tauScaleFactor);
//              mon.fillHisto("ptabsetaSelectedTau",         tauTags, tau.pt(), abs(tau.eta()), weight*tauScaleFactor);
//              mon.fillHisto("ptabsetaSelectedTauExtended", tauTags, tau.pt(), abs(tau.eta()), weight*tauScaleFactor);
//              mon.fillHisto("deltaPhi",                    tauTags, deltaPhi(selLeptons[0].phi(), tau.phi()), weight*tauScaleFactor);
//              mon.fillHisto("absDeltaPhi",                 tauTags, abs(deltaPhi(selLeptons[0].phi(), tau.phi())), weight*tauScaleFactor);
//              mon.fillHisto("deltaR",                      tauTags, deltaR(tau, selLeptons[0]), weight*tauScaleFactor);

//              mon.fillHisto("varptSelectedTau",               tauTags, tau.pt(), weight*tauScaleFactor);
//              mon.fillHisto("varptSelectedTauExtended",       tauTags, tau.pt(), weight*tauScaleFactor);
//              mon.fillHisto("varetaSelectedTau",              tauTags, tau.eta(), weight*tauScaleFactor);
//              mon.fillHisto("varptetaSelectedTau",            tauTags, tau.pt(), tau.eta(), weight*tauScaleFactor);
//              mon.fillHisto("varptetaSelectedTauExtended",    tauTags, tau.pt(), tau.eta(), weight*tauScaleFactor);
//              mon.fillHisto("varptabsetaSelectedTau",         tauTags, tau.pt(), abs(tau.eta()), weight*tauScaleFactor);
//              mon.fillHisto("varptabsetaSelectedTauExtended", tauTags, tau.pt(), abs(tau.eta()), weight*tauScaleFactor);

//              mon.fillHisto("ptSelectedLep",         tauTags, selLeptons[0].pt(), weight*tauScaleFactor);
//              mon.fillHisto("ptSelectedLepExtended", tauTags, selLeptons[0].pt(), weight*tauScaleFactor);
//              mon.fillHisto("etaSelectedLep",        tauTags, selLeptons[0].eta(), weight*tauScaleFactor);
//              mon.fillHisto("cosPhiSelectedLep",     tauTags, cos(deltaPhi(selLeptons[0].phi(), met.phi())), weight*tauScaleFactor);

              double f = 0;
              double p = 0;
/*              if(abs(selLeptons[0].id) == 11)
              {
                f = getFromTH2(etauFR, tau.pt(), tau.eta());
                p = getFromTH2(etauPR, tau.pt(), tau.eta());
              }
              else
              {
                f = getFromTH2(mutauFR, tau.pt(), tau.eta());
                p = getFromTH2(mutauPR, tau.pt(), tau.eta());
              }
              ddWeight = f*p/(p-f); // */
              ddWeight = 1;

              if(isTight)
              {
                ddWeight = f*(1-p)/(p-f);
                hasTight = true;
//                mon.fillHisto("genTauStatusTight",                tauTags, status, weight*tauScaleFactor);
                mon.fillHisto("ptSelectedTauTight",               tauTags, tau.pt(), weight*tauScaleFactor);
                mon.fillHisto("ptSelectedTauExtendedTight",       tauTags, tau.pt(), weight*tauScaleFactor);
                mon.fillHisto("etaSelectedTauTight",              tauTags, tau.eta(), weight*tauScaleFactor);
//                mon.fillHisto("cosPhiSelectedTauTight",           tauTags, cos(deltaPhi(tau.phi(), met.phi())), weight*tauScaleFactor);
//                mon.fillHisto("TauMETTight",                      tauTags, met.pt(), weight*tauScaleFactor);
//                mon.fillHisto("ptetaSelectedTauTight",            tauTags, tau.pt(), tau.eta(), weight*tauScaleFactor);
//                mon.fillHisto("ptetaSelectedTauExtendedTight",    tauTags, tau.pt(), tau.eta(), weight*tauScaleFactor);
//                mon.fillHisto("ptabsetaSelectedTauTight",         tauTags, tau.pt(), abs(tau.eta()), weight*tauScaleFactor);
//                mon.fillHisto("ptabsetaSelectedTauExtendedTight", tauTags, tau.pt(), abs(tau.eta()), weight*tauScaleFactor);
//                mon.fillHisto("deltaPhiTight",                    tauTags, deltaPhi(selLeptons[0].phi(), tau.phi()), weight*tauScaleFactor);
//                mon.fillHisto("absDeltaPhiTight",                 tauTags, abs(deltaPhi(selLeptons[0].phi(), tau.phi())), weight*tauScaleFactor);
//                mon.fillHisto("deltaRTight",                      tauTags, deltaR(tau, selLeptons[0]), weight*tauScaleFactor);

//                mon.fillHisto("varptSelectedTauTight",               tauTags, tau.pt(), weight*tauScaleFactor);
//                mon.fillHisto("varptSelectedTauExtendedTight",       tauTags, tau.pt(), weight*tauScaleFactor);
//                mon.fillHisto("varetaSelectedTauTight",              tauTags, tau.eta(), weight*tauScaleFactor);
//                mon.fillHisto("varptetaSelectedTauTight",            tauTags, tau.pt(), tau.eta(), weight*tauScaleFactor);
//                mon.fillHisto("varptetaSelectedTauExtendedTight",    tauTags, tau.pt(), tau.eta(), weight*tauScaleFactor);
//                mon.fillHisto("varptabsetaSelectedTauTight",         tauTags, tau.pt(), abs(tau.eta()), weight*tauScaleFactor);
//                mon.fillHisto("varptabsetaSelectedTauExtendedTight", tauTags, tau.pt(), abs(tau.eta()), weight*tauScaleFactor);

//                mon.fillHisto("ptSelectedLepTight",         tauTags, selLeptons[0].pt(), weight*tauScaleFactor);
//                mon.fillHisto("ptSelectedLepExtendedTight", tauTags, selLeptons[0].pt(), weight*tauScaleFactor);
//                mon.fillHisto("etaSelectedLepTight",        tauTags, selLeptons[0].eta(), weight*tauScaleFactor);
//                mon.fillHisto("cosPhiSelectedLepTight",     tauTags, cos(deltaPhi(selLeptons[0].phi(), met.phi())), weight*tauScaleFactor);
              }

//              mon.fillHisto("genTauStatus", ddTags, status, weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("ptSelectedTau", ddTags, tau.pt(), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("ptSelectedTauExtended", ddTags, tau.pt(), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("etaSelectedTau", ddTags, tau.eta(), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("cosPhiSelectedTau", ddTags, cos(deltaPhi(tau.phi(), met.phi())), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("TauMET", ddTags, met.pt(), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("ptetaSelectedTau", ddTags, tau.pt(), tau.eta(), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("ptetaSelectedTauExtended", ddTags, tau.pt(), tau.eta(), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("ptabsetaSelectedTau", ddTags, tau.pt(), abs(tau.eta()), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("ptabsetaSelectedTauExtended", ddTags, tau.pt(), abs(tau.eta()), weight*ddWeight*tauScaleFactor);

//              mon.fillHisto("varptSelectedTau", ddTags, tau.pt(), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("varptSelectedTauExtended", ddTags, tau.pt(), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("varetaSelectedTau", ddTags, tau.eta(), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("varptetaSelectedTau", ddTags, tau.pt(), tau.eta(), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("varptetaSelectedTauExtended", ddTags, tau.pt(), tau.eta(), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("varptabsetaSelectedTau", ddTags, tau.pt(), abs(tau.eta()), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("varptabsetaSelectedTauExtended", ddTags, tau.pt(), abs(tau.eta()), weight*ddWeight*tauScaleFactor);

//              mon.fillHisto("ptSelectedLep", ddTags, selLeptons[0].pt(), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("ptSelectedLepExtended", ddTags, selLeptons[0].pt(), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("etaSelectedLep", ddTags, selLeptons[0].eta(), weight*ddWeight*tauScaleFactor);
//              mon.fillHisto("cosPhiSelectedLep", ddTags, cos(deltaPhi(selLeptons[0].phi(), met.phi())), weight*ddWeight*tauScaleFactor);

              //Keep only leading tau
              //break;
            }

//            mon.fillHisto("MET", chTags, met.pt(), weight);
            if(hasTight)
            {
//              mon.fillHisto("METTight", chTags, met.pt(), weight);
              istight = true;
            }

          }
        }
      }
    }

//    mon.fillHisto("genStatus", chTags, genStatus, weight);
    if(istight)
    {
//        mon.fillHisto("genStatusTight", chTags, genStatus, weight);
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

double getFromTH2(TH2* hist, double x, double y)
{
  if(x < hist->GetXaxis()->GetXmin())
    x = hist->GetXaxis()->GetXmin();
  if(x > hist->GetXaxis()->GetXmax())
    x = hist->GetXaxis()->GetXmax();
  if(y < hist->GetYaxis()->GetXmin())
    y = hist->GetYaxis()->GetXmin();
  if(y > hist->GetYaxis()->GetXmax())
    y = hist->GetYaxis()->GetXmax();

  Int_t bin = hist->FindBin(x, y);

  return hist->GetBinContent(bin);
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
