#ifndef llvvObjects_H
#define llvvObjects_H

#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "Math/LorentzVector.h"
#include <vector>
#include <boost/cstdint.hpp>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVectorF;




class llvvGenEvent 
{
   public:
   // constructor
   llvvGenEvent(){};
   ~llvvGenEvent(){};

   //member variables
   public:
  int ngenITpu, ngenOOTpu, ngenOOTpum1, ngenTruepu;
  float pthat, genWeight, qscale, x1,x2;
  int id1, id2, nup;
};
typedef  std::vector<llvvGenEvent> llvvGenEventCollection;
typedef  edm::Ref<llvvGenEventCollection> llvvGenEventRef;
typedef  edm::RefProd<llvvGenEventCollection> llvvGenEventRefProd;
typedef  edm::RefVector<llvvGenEventCollection> llvvGenEventRefVector;


class llvvGenParticle : public LorentzVectorF
{
   public:
   // constructor
   llvvGenParticle(){};
   ~llvvGenParticle(){};

   //member variables
   public:
   float lxy;
   int id, status;
};
typedef  std::vector<llvvGenParticle> llvvGenParticleCollection;
typedef  edm::Ref<llvvGenParticleCollection> llvvGenParticleRef;
typedef  edm::RefProd<llvvGenParticleCollection> llvvGenParticleRefProd;
typedef  edm::RefVector<llvvGenParticleCollection> llvvGenParticleRefVector;


class llvvMuonInfo {
   public:
   // constructor
   llvvMuonInfo(){};
/*   llvvMuonInfo(llvvMuonInfo val): nMatches(val.nMatches), nMatchedStations(val.nMatchedStations), validMuonHits(val.validMuonHits),
                                   innerTrackChi2(val.innerTrackChi2), trkLayersWithMeasurement(val.trkLayersWithMeasurement),
                                   pixelLayersWithMeasurement(val.pixelLayersWithMeasurement) {};// */
   ~llvvMuonInfo(){};

   //member variables
   public:
   float nMatches, nMatchedStations, validMuonHits, innerTrackChi2, trkLayersWithMeasurement, pixelLayersWithMeasurement;

};
typedef  std::vector<llvvMuonInfo> llvvMuonInfoCollection;
typedef  edm::Ref<llvvMuonInfoCollection> llvvMuonInfoRef;
typedef  edm::RefProd<llvvMuonInfoCollection> llvvMuonInfoRefProd;
typedef  edm::RefVector<llvvMuonInfoCollection> llvvMuonInfoRefVector;


class llvvElectronInfo {
   public:
   // constructor
   llvvElectronInfo(){};
/*   llvvElectronInfo(llvvElectronInfo val): isConv(val.isConv), hoe(val.hoe), h2te(val.h2te), dphiin(val.dphiin), detain(val.detain),
                                           sihih(val.sihih), sipip(val.sipip), sihip(val.sihip), eopin(val.eopin), eopout(val.eopout),
                                           r9(val.r9), fbrem(val.fbrem), sce(val.sce), sceta(val.sceta), scphi(val.scphi), ooemoop(val.ooemoop),
                                           mvatrigv0(val.mvatrigv0), mvanontrigv0(val.mvanontrigv0) {};// */
   ~llvvElectronInfo(){};

   //member variables
   public:
   bool  isConv;
   float hoe,h2te,dphiin,detain,sihih,sipip,sihip, eopin, eopout,r9,fbrem;
   float sce,sceta,scphi, ooemoop;
   float mvatrigv0, mvanontrigv0;
};
typedef  std::vector<llvvElectronInfo> llvvElectronInfoCollection;
typedef  edm::Ref<llvvElectronInfoCollection> llvvElectronInfoRef;
typedef  edm::RefProd<llvvElectronInfoCollection> llvvElectronInfoRefProd;
typedef  edm::RefVector<llvvElectronInfoCollection> llvvElectronInfoRefVector;


class llvvLepton : public LorentzVectorF
{
   public:
   // constructor
   llvvLepton(){};
/*   llvvLepton(llvvLepton val): LorentzVectorF(val), id(val.id), idbits(val.idbits), genid(val.genid),
                               Tbits(val.Tbits), isPF(val.isPF), gen(val.gen), ecalIso03(val.ecalIso03),
                               hcalIso03(val.hcalIso03), trkIso03(val.trkIso03), gIso03(val.gIso03),
                               chIso03(val.chIso03), puchIso03(val.puchIso03), nhIso03(val.nhIso03),
                               ecalIso04(val.ecalIso04), hcalIso04(val.hcalIso04), trkIso04(val.trkIso04),
                               gIso04(val.gIso04), chIso04(val.chIso04), puchIso04(val.puchIso04),
                               nhIso04(val.nhIso04), d0(val.d0), dZ(val.dZ), ip3d(val.ip3d), ip3dsig(val.ip3dsig),
                               trkchi2(val.trkchi2), trkValidPixelHits(val.trkValidPixelHits),
                               trkValidTrackerHits(val.trkValidTrackerHits), trkLostInnerHits(val.trkLostInnerHits),
                               trkPtErr(val.trkPtErr), trk(val.trk), muonInfoRef(val.muonInfoRef),
                               electronInfoRef(val.electronInfoRef){};// */
   ~llvvLepton(){};

   llvvLepton& operator*=(const double& val) {LorentzVectorF::operator*=(val); return *this;};
   llvvLepton operator*(double val) {return llvvLepton(*this) *= val;};

   //member variables
   public:
   int id,          idbits,    genid, Tbits;//, pid;
   int isPF;
   LorentzVectorF gen;
   float ecalIso03, hcalIso03, trkIso03;
   float gIso03,    chIso03,   puchIso03, nhIso03; 
   float ecalIso04, hcalIso04, trkIso04;
   float gIso04,    chIso04,   puchIso04, nhIso04; 
   float d0,        dZ,        ip3d,      ip3dsig;
   float trkchi2, trkValidPixelHits, trkValidTrackerHits, trkLostInnerHits, trkPtErr;
   LorentzVectorF trk;

   //Specific Lepton Information
   llvvMuonInfoRef     muonInfoRef;
   llvvElectronInfoRef electronInfoRef;

   //functions

};
typedef  std::vector<llvvLepton> llvvLeptonCollection;
typedef  edm::Ref<llvvLeptonCollection> llvvLeptonRef;
typedef  edm::RefProd<llvvLeptonCollection> llvvLeptonRefProd;
typedef  edm::RefVector<llvvLeptonCollection> llvvLeptonRefVector;

class llvvTau  : public LorentzVectorF {
   public:
   // constructor
   llvvTau(){};
/*   llvvTau(llvvTau val): id(val.id), genid(val.genid), Tbits(val.Tbits), isPF(val.isPF), idbits(val.idbits), gen(val.gen),
                         d0(val.d0), dZ(val.dZ), ip3d(val.ip3d), ip3dsig(val.ip3dsig), dxy(val.dxy), dxySig(val.dxySig),
                         flightLength(val.flightLength), flightLengthSig(val.flightLengthSig), hasSV(val.hasSV),
                         trkchi2(val.trkchi2), trkValidPixelHits(val.trkValidPixelHits), trkValidTrackerHits(val.trkValidTrackerHits),
                         trkLostInnerHits(val.trkLostInnerHits), trkPtErr(val.trkPtErr), tracks(val.tracks), pi0s(val.pi0s),
                         vz(val.vz), z_expo(val.z_expo), emfraction(val.emfraction), hcalEnergy(val.hcalEnergy),
                         ecalEnergy(val.ecalEnergy), jet(val.jet), numChargedParticlesSigCone(val.numChargedParticlesSigCone),
                         numNeutralHadronsSigCone(val.numNeutralHadronsSigCone), numPhotonsSigCone(val.numPhotonsSigCone),
                         numPiZeroSigCone(val.numPiZeroSigCone), numParticlesSigCone(val.numParticlesSigCone),
                         numChargedParticlesIsoCone(val.numChargedParticlesIsoCone), numNeutralHadronsIsoCone(val.numNeutralHadronsIsoCone),
                         numPhotonsIsoCone(val.numPhotonsIsoCone), numParticlesIsoCone(val.numParticlesIsoCone),
                         ptSumChargedParticlesIsoCone(val.ptSumChargedParticlesIsoCone), ptSumPhotonsIsoCone(val.ptSumPhotonsIsoCone),
                         mva_e_pi(val.mva_e_pi), mva_pi_mu(val.mva_pi_mu), mva_e_mu(val.mva_e_mu) {};// */
   ~llvvTau(){};

   llvvTau& operator*=(const double& val) {LorentzVectorF::operator*=(val); return *this;};
   llvvTau operator*(double val) {return llvvTau(*this) *= val;};

   //member variables
   public:
   int id,  genid, Tbits;
   bool isPF;
   uint64_t idbits;
   LorentzVectorF gen;
   float d0,        dZ,        ip3d,      ip3dsig;
   float dxy, dxySig, flightLength, flightLengthSig;
   bool hasSV;

   float trkchi2, trkValidPixelHits, trkValidTrackerHits, trkLostInnerHits, trkPtErr;
   std::vector<LorentzVectorF> tracks;
   std::vector<LorentzVectorF> pi0s;

   float vz, z_expo;
   float emfraction, hcalEnergy, ecalEnergy;
   LorentzVectorF jet;
   int   numChargedParticlesSigCone, numNeutralHadronsSigCone, numPhotonsSigCone, numPiZeroSigCone, numParticlesSigCone;
   int   numChargedParticlesIsoCone, numNeutralHadronsIsoCone, numPhotonsIsoCone,                   numParticlesIsoCone;
   float ptSumChargedParticlesIsoCone, ptSumPhotonsIsoCone;
   float mva_e_pi, mva_pi_mu, mva_e_mu;
 
   //functions
   bool passId(uint64_t IdBit){return ((idbits&((uint64_t)1<<IdBit))>0);}
};
typedef  std::vector<llvvTau> llvvTauCollection;
typedef  edm::Ref<llvvTauCollection> llvvTauRef;
typedef  edm::RefProd<llvvTauCollection> llvvTauRefProd;
typedef  edm::RefVector<llvvTauCollection> llvvTauRefVector;



class llvvPhoton : public LorentzVectorF
{
   public:
   // constructor
   llvvPhoton(){};
   ~llvvPhoton(){};

   //member variables
   public:
   int idbits, pid, Tbits;
   float ecalIso03, hcalIso03, trkIso03, gIso03,    chIso03,   puchIso03, nhIso03;
   float ecalIso04, hcalIso04, trkIso04, gIso04,    chIso04,   puchIso04, nhIso04;

   bool  isConv;
   float hoe,h2te,sihih,sipip,sihip,r9;
   float sce,sceta,scphi;
};
typedef  std::vector<llvvPhoton> llvvPhotonCollection;
typedef  edm::Ref<llvvPhotonCollection> llvvPhotonRef;
typedef  edm::RefProd<llvvPhotonCollection> llvvPhotonRefProd;
typedef  edm::RefVector<llvvPhotonCollection> llvvPhotonRefVector;

class llvvJet : public LorentzVectorF
{
   public:
   // constructor
   llvvJet(){};
/*   llvvJet(llvvJet val): idbits(val.idbits), pfstart(val.pfstart), pfend(val.pfend), Tbits(val.Tbits), torawsf(val.torawsf),
                         neutHadFrac(val.neutHadFrac), neutEmFrac(val.neutEmFrac), chHadFrac(val.chHadFrac),
                         muFrac(val.muFrac), area(val.area), tchp(val.tchp), jp(val.jp), origcsv(val.origcsv), csv(val.csv),
                         jpcsv(val.jpcsv), slcsv(val.slcsv), supercsv(val.supercsv), ssvhe(val.ssvhe), ivf(val.ivf),
                         svxPx(val.svxPx), svxPy(val.svxPy), svxPz(val.svxPz), svxM(val.svxM), svxNtrk(val.svxNtrk),
                         svxLxy(val.svxLxy), svxLxyErr(val.svxLxyErr), ivfPx(val.ivfPx), ivfPy(val.ivfPy), ivfPz(val.ivfPz),
                         ivfM(val.ivfM), ivfNtrk(val.ivfNtrk), ivfLxy(val.ivfLxy), ivfLxyErr(val.ivfLxyErr), puMVA(val.puMVA),
                         qgMVA(val.qgMVA), beta(val.beta), betaStar(val.betaStar), dRMean(val.dRMean), dR2Mean(val.dR2Mean),
                         ptRMS(val.ptRMS), ptD(val.ptD), etaW(val.etaW), phiW(val.phiW), genflav(val.genflav)
                         genid(val.genid), gen(val.gen), genj(val.genj\) {};// */
   ~llvvJet(){};

   llvvJet& operator*=(const double& val) {LorentzVectorF::operator*=(val); return *this;};
   llvvJet operator*(double val) {return llvvJet(*this) *= val;};

   //member variables
   public:
   int idbits, pfstart, pfend, Tbits;
   float torawsf;
   float neutHadFrac, neutEmFrac, chHadFrac, muFrac, area;
   float tchp, jp, origcsv, csv, jpcsv, slcsv, supercsv, ssvhe, ivf;
   float svxPx, svxPy, svxPz, svxM, svxNtrk, svxLxy, svxLxyErr;
   float ivfPx, ivfPy, ivfPz, ivfM, ivfNtrk, ivfLxy, ivfLxyErr;
   float puMVA, qgMVA;
   float beta,betaStar, dRMean, dR2Mean, ptRMS,ptD, etaW, phiW;
   int   genflav, genid;
   LorentzVectorF gen, genj;
};
typedef  std::vector<llvvJet> llvvJetCollection;
typedef  edm::Ref<llvvJetCollection> llvvJetRef;
typedef  edm::RefProd<llvvJetCollection> llvvJetRefProd;
typedef  edm::RefVector<llvvJetCollection> llvvJetRefVector;

class llvvMet : public LorentzVectorF
{
   public:
   // constructor
   llvvMet(){};
   ~llvvMet(){};

   llvvMet& operator*=(const double& val) {LorentzVectorF::operator*=(val); return *this;};
   llvvMet operator*(double val) {return llvvMet(*this) *= val;};

   //member variables
   public:
   float sig,sigx2,sigxy,sigy2;
};
typedef  std::vector<llvvMet> llvvMetCollection;
typedef  edm::Ref<llvvMetCollection> llvvMetRef;
typedef  edm::RefProd<llvvMetCollection> llvvMetRefProd;
typedef  edm::RefVector<llvvMetCollection> llvvMetRefVector;


class llvvPFParticle : public LorentzVectorF
{
   public:
   // constructor
   llvvPFParticle(){};
   ~llvvPFParticle(){};

   //member variables
   public:
   int id, charge;
};
typedef  std::vector<llvvPFParticle> llvvPFParticleCollection;
typedef  edm::Ref<llvvPFParticleCollection> llvvPFParticleRef;
typedef  edm::RefProd<llvvPFParticleCollection> llvvPFParticleRefProd;
typedef  edm::RefVector<llvvPFParticleCollection> llvvPFParticleRefVector;




///////////////////////////////////////////////////////////////////////////////////////////////////////////
// // //  NOT USED TO SAVE THE OBJECTS IN THE EDM FORMAT BUT CONVENIENT WHEN ANALYZING THE OBJECTS // // //
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//Just an ugly way to save additional info to the object
struct llvvJetInfo{};
class llvvJetExt : public llvvJet{
   public:
   llvvJetExt(): llvvJet(){
   };
   llvvJetExt(llvvJet jet_): llvvJet(jet_), jer(1), jerup(1), jerdown(1), jesup(1), jesdown(1) {};
/*   llvvJetExt(llvvJet jet_){
      SetPxPyPzE(jet_.px(), jet_.py(), jet_.pz(), jet_.energy());
      idbits=jet_.idbits; pfstart=jet_.pfstart; pfend=jet_.pfend;
      Tbits=jet_.Tbits;
      torawsf=jet_.torawsf;
      neutHadFrac=jet_.neutHadFrac; neutEmFrac=jet_.neutEmFrac; chHadFrac=jet_.chHadFrac; muFrac=jet_.muFrac; area=jet_.area;
      tchp=jet_.tchp; jp=jet_.jp; origcsv=jet_.origcsv; csv=jet_.csv; jpcsv=jet_.jpcsv; slcsv=jet_.slcsv; supercsv=jet_.supercsv; ssvhe=jet_.ssvhe; ivf=jet_.ivf;
      svxPx=jet_.svxPx; svxPy=jet_.svxPy; svxPz=jet_.svxPz; svxM=jet_.svxM; svxNtrk=jet_.svxNtrk; svxLxy=jet_.svxLxy; svxLxyErr=jet_.svxLxyErr;
      ivfPx=jet_.ivfPx; ivfPy=jet_.ivfPy; ivfPz=jet_.ivfPz; ivfM=jet_.ivfM; ivfNtrk=jet_.ivfNtrk; ivfLxy=jet_.ivfLxy; ivfLxyErr=jet_.ivfLxyErr;
      puMVA=jet_.puMVA; qgMVA=jet_.qgMVA;
      beta=jet_.beta; betaStar=jet_.betaStar; dRMean=jet_.dRMean; dR2Mean=jet_.dR2Mean; ptRMS=jet_.ptRMS; ptD=jet_.ptD; etaW=jet_.etaW; phiW=jet_.phiW;
      genflav=jet_.genflav; genid=jet_.genid;
      gen=jet_.gen; genj=jet_.genj;
   };// */
   ~llvvJetExt(){};

   llvvJetExt& operator*=(const double& val) {LorentzVectorF::operator*=(val); return *this;};
   llvvJetExt operator*(double val) {return llvvJetExt(*this) *= val;};

   public: 
   double jer; double jerup; double jerdown; double jesup; double jesdown;
};
typedef  std::vector<llvvJetExt> llvvJetExtCollection;

//ALL Sorting functions
inline bool sort_llvvObjectByPt(const LorentzVectorF &a, const LorentzVectorF &b)  { return a.pt()>b.pt(); }
inline bool sort_llvvJetByCSV(const llvvJet &a, const llvvJet &b) { return a.supercsv>b.supercsv; }

//ONLY ADD STUFF AT THE END... CAN HOST UP TO 64 VARIABLES
enum llvvTAUID {
   decayModeFindingNewDMs,
   decayModeFindingOldDMs,
   decayModeFinding,
//   byLooseIsolation,
//   byVLooseCombinedIsolationDeltaBetaCorr,
//   byLooseCombinedIsolationDeltaBetaCorr,
   byMediumCombinedIsolationDeltaBetaCorr,
   byTightCombinedIsolationDeltaBetaCorr,
   byCombinedIsolationDeltaBetaCorrRaw,
   byLooseCombinedIsolationDeltaBetaCorr3Hits,
   byMediumCombinedIsolationDeltaBetaCorr3Hits,
   byTightCombinedIsolationDeltaBetaCorr3Hits,
   byCombinedIsolationDeltaBetaCorrRaw3Hits,
   chargedIsoPtSum,
   neutralIsoPtSum,
   puCorrPtSum,
   byIsolationMVA3oldDMwoLTraw,
   byVLooseIsolationMVA3oldDMwoLT,
   byLooseIsolationMVA3oldDMwoLT,
   byMediumIsolationMVA3oldDMwoLT,
   byTightIsolationMVA3oldDMwoLT,
   byVTightIsolationMVA3oldDMwoLT,
   byVVTightIsolationMVA3oldDMwoLT,
   byIsolationMVA3oldDMwLTraw,
   byVLooseIsolationMVA3oldDMwLT,
   byLooseIsolationMVA3oldDMwLT,
   byMediumIsolationMVA3oldDMwLT,
   byTightIsolationMVA3oldDMwLT,
   byVTightIsolationMVA3oldDMwLT,
   byVVTightIsolationMVA3oldDMwLT,
   byIsolationMVA3newDMwoLTraw,
   byVLooseIsolationMVA3newDMwoLT,
   byLooseIsolationMVA3newDMwoLT,
   byMediumIsolationMVA3newDMwoLT,
   byTightIsolationMVA3newDMwoLT,
   byVTightIsolationMVA3newDMwoLT,
   byVVTightIsolationMVA3newDMwoLT,
   byIsolationMVA3newDMwLTraw,
   byVLooseIsolationMVA3newDMwLT,
   byLooseIsolationMVA3newDMwLT,
   byMediumIsolationMVA3newDMwLT,
   byTightIsolationMVA3newDMwLT,
   byVTightIsolationMVA3newDMwLT,
   byVVTightIsolationMVA3newDMwLT,
   againstElectronLoose,
   againstElectronMedium,
   againstElectronTight,
   againstElectronMVA5raw,
   againstElectronMVA5category,
   againstElectronVLooseMVA5,
   againstElectronLooseMVA5,
   againstElectronMediumMVA5,
   againstElectronTightMVA5,
   againstElectronVTightMVA5,
   againstElectronDeadECAL,
   againstMuonLoose,
   againstMuonMedium,
   againstMuonTight,
   againstMuonLoose2,
   againstMuonMedium2,
   againstMuonTight2,
   againstMuonLoose3,
   againstMuonTight3,
   againstMuonMVAraw,
   againstMuonLooseMVA,
   againstMuonMediumMVA,
   againstMuonTightMVA
};





class llvvTauLepton  : public llvvLepton {
   public:
   // constructor
   ~llvvTauLepton(){};
   llvvTauLepton(llvvLepton lep_){lep = lep_;  SetPxPyPzE(lep_.px(), lep_.py(), lep_.pz(), lep_.energy()); id=lep_.id; };
   llvvTauLepton(llvvTau    tau_){tau = tau_;  SetPxPyPzE(tau_.px(), tau_.py(), tau_.pz(), tau_.energy()); id=tau_.id; };
   llvvLepton lep;
   llvvTau tau;
};
typedef  std::vector<llvvTauLepton> llvvTauLeptonCollection;
typedef  edm::Ref<llvvTauLeptonCollection> llvvTauLeptonRef;
typedef  edm::RefProd<llvvTauLeptonCollection> llvvTauLeptonRefProd;
typedef  edm::RefVector<llvvTauLeptonCollection> llvvTauLeptonRefVector;





#endif
