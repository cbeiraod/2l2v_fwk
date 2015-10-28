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
#include <algorithm>

// Include MT2 library:
// http://particle.physics.ucdavis.edu/hefti/projects/doku.php?id=wimpmass    ** Code from here
// http://www.hep.phy.cam.ac.uk/~lester/mt2/    ** Other libraries
#include "UserCode/llvv_fwk/interface/mt2_bisect.h"


//#define WITH_UNLOCK
#include "UserCode/llvv_fwk/interface/ValueWithSystematics.h"
#include "UserCode/llvv_fwk/interface/llvvAnalyser.h"


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

class PepAnalyser : public Analyser
{
public:
  PepAnalyser(std::string cfgFile);
  
  enum class IDType {LooseID, MediumID, TightID};
  enum class TAU_E_ID {antiELoose, antiEMedium, antiETight, antiEMva, antiEMva3Loose, antiEMva3Medium, antiEMva3Tight, antiEMva3VTight, antiEMva5Medium};

private:

protected:
  bool exclusiveRun;

  double minElPt;
  double maxElEta;
  double ECALGap_MinEta;
  double ECALGap_MaxEta;
  double minMuPt;
  double maxMuEta;
  double minJetPt;
  double maxJetEta;
  double maxElDz;
  double maxElD0;
  double maxMuDz;
  double maxMuD0;
  double elIso;
  double muIso;
  double genMatchRCone;

  virtual void UserLoadCfgOptions();
  virtual void UserSetup();
  virtual void UserProcessEvent(size_t iev);
  virtual void UserInitHistograms();
  virtual void UserEventContentSetup();
  virtual void UserFillHistograms();

  bool electronMVAID(double mva, llvvLepton& lepton, IDType id);

  int isZTauTau();
};

bool PepAnalyser::electronMVAID(double mva, llvvLepton& lepton, IDType id)
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

int PepAnalyser::isZTauTau()
{
  /*for(auto& genPart : gen)
  {
    if(abs(genPart.id) == 23) // If a Z boson
    {
      int nTau = 0;

      for(int i = 0; i < genPart.numberOfDaughters(); ++i)
      {
        if(abs(genPart.daughter(i)->pdgId()) == 15)
          nTau++;
      }

      if(nTau > 0)
        return nTau;
    }
  }// */

  return 0;
}

PepAnalyser::PepAnalyser(std::string cfgFile): Analyser(cfgFile)
{
}

void PepAnalyser::UserLoadCfgOptions()
{
  exclusiveRun = cfgOptions.getParameter<bool>("exclusiveRun");

  // Consider setting here the cut values etc, will have to be added to the cfg file
  minElPt        = 30;      // Selected electron pT and eta
  maxElEta       =  2.1;
  ECALGap_MinEta =  1.4442; // ECAL gap parameters
  ECALGap_MaxEta =  1.5660;
  minMuPt        = 20;      // Selected muon pT and eta
  maxMuEta       =  2.1;
  minJetPt       = 30;
  maxJetEta      =  4.7;    // Selected jet eta

  maxElDz        =  0.1;
  maxElD0        =  0.045;

  maxMuDz        =  0.5;
  maxMuD0        =  0.2;

  elIso          =  0.1;
  muIso          =  0.1;

  genMatchRCone  =  0.3;

  if(debug)
    std::cout << "Finished PepAnalyser::LoadCfgOptions()" << std::endl;

  return;
}

void PepAnalyser::UserSetup()
{
  //TDirectory* cwd = gDirectory;

  if(debug)
    std::cout << "Finished PepAnalyser::UserSetup()" << std::endl;

  return;
}

void PepAnalyser::UserProcessEvent(size_t iev)
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

  bool isDYTauTau = false;
  int nTau = isZTauTau();
  if(nTau > 0)
    isDYTauTau = true;

  if((mctruthmode == 0 && isDYTauTau) || (mctruthmode == 1 && !isDYTauTau))
    dropEvent = true;

  bool singleETrigger  = triggerBits[13]; // HLT_Ele27_WP80_v*
  bool singleMuTrigger = triggerBits[15]; // HLT_IsoMu24_v*
  //bool TauPlusE2012A  = triggerBits[18]; // HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*
  //bool TauPlusMu2012A = triggerBits[22]; // HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v*
  //bool TauPlusE2012B  = triggerBits[17]; // HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v*
  //bool TauPlusMu2012B = triggerBits[21]; // HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*

  //bool TauPlusETrigger = TauPlusE2012A || TauPlusE2012B;
  //bool TauPlusMuTrigger = TauPlusMu2012A || TauPlusMu2012B;
  auto& triggeredOn = eventContent.GetBool("triggeredOn");
  triggeredOn = singleETrigger || singleMuTrigger;

  // Get Leptons
  if(debugEvent)
  {
    analyserCout << " Finished computing PU weight" << std::endl;
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
      muCor->applyPtCorrection(p4, (lep.id>0)?1:-1);
      if(isMC)
        muCor->applyPtSmearing(p4, (lep.id>0)?1:-1, false);
      lep.SetPxPyPzE(p4.Px(), p4.Py(), p4.Pz(), p4.Energy());
    }

    // Lepton Kinematics
    double eta = (lepId == 11)?(lep.electronInfoRef->sceta):(lep.eta());
    ValueWithSystematics<bool> passKin(true);
    if(lepId == 11) // If Electron
    {
      if(std::abs(eta) > maxElEta)
        passKin = false;

      if(std::abs(eta) > ECALGap_MinEta && std::abs(eta) < ECALGap_MaxEta) // Remove electrons that fall in ECAL Gap
      {
        passKin = false;
      }
      
      if(runSystematics && isMC)
      {
        passKin("LES_UP");
        passKin("LES_DOWN");
      }
      if(lep.pt() < minElPt)
        passKin.Value() = false;
      if(runSystematics && isMC)
      {
        if(lep.pt()*(1+sf) < minElPt)
          passKin("LES_UP") = false;
        if(lep.pt()*(1-sf) < minElPt)
          passKin("LES_DOWN") = false;
      }
    }
    else // If Muon
    {
      if(std::abs(eta) > maxMuEta)
        passKin = false;

      if(runSystematics && isMC)
      {
        passKin("LES_UP");
        passKin("LES_DOWN");
      }
      if(lep.pt() < minMuPt)
        passKin.Value() = false;
      if(runSystematics && isMC)
      {
        if(lep.pt()*(1+sf) < minMuPt)
          passKin("LES_UP") = false;
        if(lep.pt()*(1-sf) < minMuPt)
          passKin("LES_DOWN") = false;
      }
    }

    // Lepton ID
    bool passID = true;
    Int_t idbits = lep.idbits;
    if(lepId == 11)
    {
      // bool isTight = electronMVAID(lep.electronInfoRef->mvanontrigv0, lep, IDType::TightID);
      // bool isLoose = electronMVAID(lep.electronInfoRef->mvanontrigv0, lep, IDType::LooseID);
      // bool isLoose = ((idbits >> 4) & 0x1);
      // bool isTight = ((idbits >> 6) & 0x1);
      passID = electronMVAID(lep.electronInfoRef->mvanontrigv0, lep, IDType::LooseID);
      if(lep.d0 > maxElD0)
        passID = false;
      if(lep.dZ > maxElDz)
        passID = false;

      if(lep.electronInfoRef->isConv)
      {
        passID = false;
      }
      if(lep.trkLostInnerHits > 0)
      {
        passID = false;
      }
    }
    else
    {
      // bool isLoose = ((idbits >> 8) & 0x1);
      // bool isTight = ((idbits >> 10) & 0x1);
      passID = ((idbits >> 10) & 0x1);
      if(lep.d0 > maxMuD0)
        passID = false;
      if(lep.dZ > maxMuDz)
        passID = false;
    }

    // Lepton Isolation
    bool passIso = true;
    double relIso = utils::cmssw::relIso(lep, eventContent.GetDouble("rho").Value());
    if(lepId == 11)
    {
      if(relIso > elIso)
        passIso = false;
    }
    else
    {
      if(relIso > muIso)
        passIso = false;
    }

    // Keep desired leptons
    if(static_cast<bool>(passKin) && passID && passIso)
    {
      if(runSystematics && isMC)
      {
        if(passKin.GetSystematicOrValue("LES_UP"))
          selLeptons.Systematic("LES_UP").push_back(lep*(1+sf));
        if(passKin.GetSystematicOrValue("LES_DOWN"))
          selLeptons.Systematic("LES_DOWN").push_back(lep*(1-sf));
      }
      if(passKin.Value())
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
    //bool hasBtagCorr = false;
    if(jet.csv > 0.679)
    {
      isBJet = true;
      //hasBtagCorr = true;
    }

    ValueWithSystematics<bool> passIso(true);
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
    loadSystematics(tmpLoop, selJets);
  }
  for(auto& val: tmpLoop)
  {
    auto& lleptons = selLeptons.GetSystematicOrValue(val);
    if(lleptons.size() != 0)
      std::sort(lleptons.begin(), lleptons.end(), sort_llvvObjectByPt);

    auto& ljets = selJets.GetSystematicOrValue(val);
    if(ljets.size() != 0)
      std::sort(ljets.begin(), ljets.end(), sort_llvvObjectByPt);

    //if(val != "Value")
    //  nBJet.Systematic(val) = selBJets.GetSystematicOrValue(val).size();
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

    analyserCout << " selLeptons systematics:\n";
    for(auto& syst: selLeptons.Systematics())
      analyserCout << "    " << syst << ": " << selLeptons(syst).size() << "\n";
    analyserCout << std::endl;
  }

  auto& nBJet = eventContent.GetInt("nBJet");
  nBJet = selBJets.size();
  auto& nJet = eventContent.GetInt("nJet");
  nJet = selJets.size();
  auto& nLep = eventContent.GetInt("nLep");
  nLep = selLeptons.size();
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

  const ValueWithSystematics<double> unit(1);

  auto met = MET;
  auto lep = selLeptons.first().ToTLorentzVector();

  auto& cosDeltaPhiLep = eventContent.GetDouble("cosPhiLep");
  auto& MTLep = eventContent.GetDouble("MTLep");

  cosDeltaPhiLep = (lep.DeltaPhi(met)).Cos();
  ValueWithSystematics<double> fac = met.Pt() * lep.Pt() * 2;
  MTLep = fac * (unit - cosDeltaPhiLep);
  MTLep = MTLep.Sqrt();

  if(debugEvent)
    analyserCout << " Is the event selected?" << std::endl;
  auto& selected = eventContent.GetBool("selected");
  selected = triggeredOn && (nBJet >= 1) && (nJet >= 4) && (nLep >= 1) && (eventContent.GetDouble("MET") > 80) && (MTLep > 100);
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

void PepAnalyser::UserInitHistograms()
{
  histMonitor.addHistogram(new TH1D("nlep", ";# lep;Events", 10, 0, 10));
  histMonitor.addHistogram(new TH1D("nel", ";# e;Events", 10, 0, 10));
  histMonitor.addHistogram(new TH1D("nmu", ";# #mu;Events", 10, 0, 10));
  histMonitor.addHistogram(new TH1D("njets", ";# jets;Events", 6, 0, 6));
  histMonitor.addHistogram(new TH1D("nbjets", ";# jets_{b};Events", 6, 0, 6));

  // Eventflow
  TH1D *eventflow = (TH1D*)histMonitor.addHistogram(new TH1D("eventflow", ";;Events", 6, 0, 6));
  eventflow->GetXaxis()->SetBinLabel(1, "HLT");
  eventflow->GetXaxis()->SetBinLabel(2, "MET > 80");
  eventflow->GetXaxis()->SetBinLabel(3, "1 lepton");
  eventflow->GetXaxis()->SetBinLabel(4, "4 jets");
  eventflow->GetXaxis()->SetBinLabel(5, "1 b-jet");
  eventflow->GetXaxis()->SetBinLabel(6, "MT_{l} > 100");

  // Leptons
  histMonitor.addHistogram(new TH1D("ptSelectedLep", ";p_{T}^{l};Events", 50, 0, 100));
  histMonitor.addHistogram(new TH1D("etaSelectedLep", ";#eta^{l};Events", 25, -2.6, 2.6));
  histMonitor.addHistogram(new TH1D("chargeSelectedLep", ";q^{l};Events", 5, -2, 2));
  TH1D *leptonCutFlow = (TH1D*)histMonitor.addHistogram(new TH1D("leptonCutFlow", ";;Leptons", 4, 0, 4));
  leptonCutFlow->GetXaxis()->SetBinLabel(1, "All");
  leptonCutFlow->GetXaxis()->SetBinLabel(2, "ID");
  leptonCutFlow->GetXaxis()->SetBinLabel(3, "Kin");
  leptonCutFlow->GetXaxis()->SetBinLabel(4, "Iso");

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
  histMonitor.addHistogram(new TH1D("MTLep", ";MT_{l} [GeV];Events", 25, 0, 200));

  if(debug)
    std::cout << "Finished PepAnalyser::UserInitHistograms()" << std::endl;

  return;
}

void PepAnalyser::UserFillHistograms()
{
  auto& weight = eventContent.GetDouble("weight").Value();
  //auto& puWeight = eventContent.GetDouble("PUweight").Value();
  auto& selected = eventContent.GetBool("selected").Value();
  auto& dropEvent = eventContent.GetBool("dropEvent").Value();

  // Eventflow
  if(!dropEvent && eventContent.GetBool("triggeredOn").Value())
  {
    histMonitor.fillHisto("eventflow", chTags, 0, weight);
    if(eventContent.GetDouble("MET").Value() > 80)
    {
      histMonitor.fillHisto("eventflow", chTags, 1, weight);
      if(eventContent.GetInt("nLep").Value() >= 1)
      {
        histMonitor.fillHisto("eventflow", chTags, 2, weight);
        if(eventContent.GetInt("nJet").Value() >= 4)
        {
          histMonitor.fillHisto("eventflow", chTags, 3, weight);
          if(eventContent.GetInt("nBJet").Value() >= 1)
          {
            histMonitor.fillHisto("eventflow", chTags, 4, weight);
            if(eventContent.GetDouble("MTLep").Value() > 100)
            {
              histMonitor.fillHisto("eventflow", chTags, 5, weight);
            }
          }
        }
      }
    }
  }


  if(selected)
  {
    histMonitor.fillHisto("nbjets", chTags, eventContent.GetInt("nBJet").Value(), weight);
    histMonitor.fillHisto("nlep", chTags, eventContent.GetInt("nLep").Value(), weight);
    histMonitor.fillHisto("nel", chTags, eventContent.GetInt("nEl").Value(), weight);
    histMonitor.fillHisto("nmu", chTags, eventContent.GetInt("nMu").Value(), weight);
    histMonitor.fillHisto("ntau", chTags, eventContent.GetInt("nTau").Value(), weight);
    histMonitor.fillHisto("njets", chTags, eventContent.GetInt("nJet").Value(), weight);

    histMonitor.fillHisto("MET", chTags, eventContent.GetDouble("MET").Value(), weight);
    
    // Selected Lepton
    histMonitor.fillHisto("ptSelectedLep", chTags, eventContent.GetDouble("LeptonPt").Value(), weight);
    //histMonitor.fillHisto("etaSelectedLep", chTags, eventContent.GetInt("nBJet").Value(), weight);

    // MT
    histMonitor.fillHisto("MTLep", chTags, eventContent.GetDouble("MTLep").Value(), weight);
  }
}

void PepAnalyser::UserEventContentSetup()
{
  eventContent.AddBool("triggeredOn", false);
  eventContent.AddBool("dropEvent", false);

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

  eventContent.AddDouble("cosPhiLep", -999);
  auto& MTLep = eventContent.AddDouble("MTLep", -999);
  MTLep.AddMetadata("eventlist", "true");

  if(debug)
    std::cout << "Finished PepAnalyser::UserEventContentSetup()" << std::endl;

  return;
}

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
  //bool doOld = false;
  bool eventlistSelected = false;

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
      
      if(arg.find("--eventlistSelected") != std::string::npos)
      {
        eventlistSelected = true;
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

  PepAnalyser testing(argv[fileIndex]);
  if(limit != 0)
    testing.SetEventLimit(limit);
  if(keepAllEvents)
    testing.KeepAllEvents();
  if(debugEvent)
    testing.SetDebugEvent(debugEvent);
  if(skipEvents != 0)
    testing.SetSkipEvents(skipEvents);
  if(eventlistSelected)
    testing.SetEventlistSelected(eventlistSelected);
  testing.LoopOverEvents();

  return 0;
}
