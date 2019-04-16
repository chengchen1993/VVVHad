// system include files
#include <iostream>
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"  
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <DataFormats/JetReco/interface/Jet.h>
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include <TFormula.h>
#define Pi 3.141593
using namespace std;
//
// class declaration
//
class EDBRTreeMaker : public edm::EDAnalyzer {
public:
  explicit EDBRTreeMaker(const edm::ParameterSet&);
  ~EDBRTreeMaker();
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
  virtual void endRun(const edm::Run&, const edm::EventSetup&) override;

  virtual bool looseJetID( const pat::Jet& j);
  virtual const reco::Candidate* findLastW(const reco::Candidate *particle,int IDpdg);
  virtual const reco::Candidate* findFirstW(const reco::Candidate *particle);

  virtual bool tightJetID( const pat::Jet& j);
  virtual float dEtaInSeed( const pat::Electron* ele) ;
  virtual void initJetCorrFactors( void );
  virtual void addTypeICorr( edm::Event const & event );
  virtual double getJEC( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
  virtual double getJECOffset( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
  math::XYZTLorentzVector getNeutrinoP4(double& MetPt, double& MetPhi, TLorentzVector& lep, int lepType);
  std::vector<std::string>                    jecAK8PayloadNames_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK8_            ;
  std::vector<std::string>                    jecAK8PayloadNamesGroomed_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK8Groomed_            ;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK8GroomedSD_            ;
  std::vector<std::string>                    jecAK8puppiPayloadNames_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK8puppi_            ;
  std::vector<std::string>                    jecAK8puppiPayloadNamesGroomed_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK8puppiGroomed_            ;
  std::vector<std::string>                    jecAK4PayloadNames_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK4_            ;
  std::vector<std::string> offsetCorrLabel_;
  boost::shared_ptr<FactorizedJetCorrector> jecOffset_;
  edm::Handle< double >  rho_;
  edm::InputTag  METsRawLabel_;
  edm::Handle<pat::METCollection>  METs_;
  edm::Handle<pat::JetCollection> jets_;
  edm::Handle<reco::VertexCollection> vertices_;
  edm::EDGetTokenT<pat::MuonCollection> muons_;
  edm::Handle<pat::METCollection>  reclusteredMETs_;
  edm::Handle<edm::View<reco::PFMET> >     pfMET_ ;
  edm::EDGetTokenT<pat::JetCollection> prunedjetInputToken_;
  edm::EDGetTokenT<pat::JetCollection> softdropjetInputToken_;
  edm::EDGetTokenT<pat::JetCollection> fatjetInputToken_;
  edm::EDGetTokenT<pat::JetCollection> puppijetInputToken_;

  edm::EDGetTokenT<pat::METCollection>  metInputToken_;
  edm::EDGetTokenT<pat::METCollection>  reclusteredmetInputToken_;
  std::vector<edm::EDGetTokenT<pat::METCollection>> mettokens;
  edm::Handle<pat::JetCollection> prunedjets_;
  edm::Handle<pat::JetCollection> softdropjets_;
  edm::Handle<pat::JetCollection> puppijets_;

  std::vector<edm::EDGetTokenT<pat::JetCollection>> jetTokens;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<pat::METCollection> reclusteredmetToken_;
  edm::EDGetTokenT<pat::JetCollection> jetToken_;
  edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
  edm::EDGetTokenT<pat::JetCollection> prunedjetToken_;
  edm::EDGetTokenT<pat::JetCollection> softdropjetToken_;
  edm::EDGetTokenT<pat::JetCollection> puppijetToken_;
  edm::Handle<pat::JetCollection> fatjets_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  std::vector<std::string> jetCorrLabel_;
  std::vector<std::string> jecAK4Labels;
  std::vector<std::string> jecAK8Labels;
  bool doCorrOnTheFly_;
  edm::EDGetTokenT<edm::TriggerResults> 		     noiseFilterToken_;
  edm::Handle< edm::TriggerResults> 			     noiseFilterBits_;
  std::string HBHENoiseFilter_Selector_;
  std::string HBHENoiseIsoFilter_Selector_;
  std::string GlobalHaloNoiseFilter_Selector_;
  std::string ECALDeadCellNoiseFilter_Selector_;
  std::string GoodVtxNoiseFilter_Selector_;
  std::string EEBadScNoiseFilter_Selector_;
  edm::EDGetTokenT<bool> badMuonNoiseFilter_Selector_;
  edm::EDGetTokenT<bool> badChargedHadronNoiseFilter_Selector_;
  // ----------member data ---------------------------
  TTree* outTree_;
  TTree* outTreew_;
  double MW_;
  int nmetmatch, nmetno;
  int nevent, run, ls;
  int nVtx;
  int numCands;
  int nLooseEle, nLooseMu;
  int nVetoEle, nVetoMu;
  int qnumber, qnumber2, qnumber3;
  double met, metPhi;
  double jetAK8puppi_ptJEC, jetAK8puppi_eta, jetAK8puppi_phi;
  double jetAK8puppi_tau1,  jetAK8puppi_tau2, jetAK8puppi_tau3, jetAK8puppi_tau4;
  double jetAK8puppi_tau21, jetAK8puppi_tau42;
  double jetAK8puppi_sd, jetAK8puppi_sdJEC, jetAK8puppi_sdcorr;
  double jetAK8puppi_ptJEC_2, jetAK8puppi_eta_2, jetAK8puppi_phi_2;
  double jetAK8puppi_tau1_2,  jetAK8puppi_tau2_2, jetAK8puppi_tau3_2, jetAK8puppi_tau4_2;    
  double jetAK8puppi_tau21_2, jetAK8puppi_tau42_2;
  double jetAK8puppi_sd_2, jetAK8puppi_sdJEC_2, jetAK8puppi_sdcorr_2;
  double jetAK8puppi_ptJEC_3, jetAK8puppi_eta_3, jetAK8puppi_phi_3;
  double jetAK8puppi_tau1_3,  jetAK8puppi_tau2_3, jetAK8puppi_tau3_3, jetAK8puppi_tau4_3;    
  double jetAK8puppi_tau21_3, jetAK8puppi_tau42_3;
  double jetAK8puppi_sd_3, jetAK8puppi_sdJEC_3, jetAK8puppi_sdcorr_3;
  double jetAK8puppi_ptJEC_4, jetAK8puppi_eta_4, jetAK8puppi_phi_4;
  double jetAK8puppi_tau1_4,  jetAK8puppi_tau2_4, jetAK8puppi_tau3_4, jetAK8puppi_tau4_4;
  double jetAK8puppi_tau21_4, jetAK8puppi_tau42_4;
  double jetAK8puppi_sd_4, jetAK8puppi_sdJEC_4, jetAK8puppi_sdcorr_4;
  double jetAK8puppi_ptJEC_5, jetAK8puppi_eta_5, jetAK8puppi_phi_5;
  double jetAK8puppi_tau1_5,  jetAK8puppi_tau2_5, jetAK8puppi_tau3_5, jetAK8puppi_tau4_5;
  double jetAK8puppi_tau21_5, jetAK8puppi_tau42_5;
  double jetAK8puppi_sd_5, jetAK8puppi_sdJEC_5, jetAK8puppi_sdcorr_5;
  double jetAK8puppi_ptJEC_6, jetAK8puppi_eta_6, jetAK8puppi_phi_6;
  double jetAK8puppi_tau1_6,  jetAK8puppi_tau2_6, jetAK8puppi_tau3_6, jetAK8puppi_tau4_6;
  double jetAK8puppi_tau21_6, jetAK8puppi_tau42_6;
  double jetAK8puppi_sd_6, jetAK8puppi_sdJEC_6, jetAK8puppi_sdcorr_6;
  double jetAK8puppi_ptJEC_7, jetAK8puppi_eta_7, jetAK8puppi_phi_7;
  double jetAK8puppi_tau1_7,  jetAK8puppi_tau2_7, jetAK8puppi_tau3_7, jetAK8puppi_tau4_7;
  double jetAK8puppi_tau21_7, jetAK8puppi_tau42_7;
  double jetAK8puppi_sd_7, jetAK8puppi_sdJEC_7, jetAK8puppi_sdcorr_7;
  double jetAK8puppi_ptJEC_8, jetAK8puppi_eta_8, jetAK8puppi_phi_8;
  double jetAK8puppi_tau1_8,  jetAK8puppi_tau2_8, jetAK8puppi_tau3_8, jetAK8puppi_tau4_8;
  double jetAK8puppi_tau21_8, jetAK8puppi_tau42_8;
  double jetAK8puppi_sd_8, jetAK8puppi_sdJEC_8, jetAK8puppi_sdcorr_8;
  double triggerWeight, lumiWeight, pileupWeight;
  double ptq11, etaq11, phiq11, massq11;
  double ptq21, etaq21, phiq21, massq21;
  double ptq31, etaq31, phiq31, massq31;
  double ptq12, etaq12, phiq12, massq12;
  double ptq22, etaq22, phiq22, massq22;
  double ptq32, etaq32, phiq32, massq32;
  double delPhijetmet, delPhijetmet_2, delPhijetmet_3, delPhijetmet_4, delPhijetmet_5, delPhijetmet_6, delPhijetmet_7, delPhijetmet_8;
  double pt_graviton,pt_graviton1;
  double candMasspuppiJEC;
  double massww[3];
  double pweight[211];
    TLorentzVector ak8sj11,ak8sj21,ak8sj31,ak8sj41,ak8sj51;
    TLorentzVector ak8sj12,ak8sj22,ak8sj32,ak8sj42,ak8sj52;
    TLorentzVector ak8sj13,ak8sj23,ak8sj33,ak8sj43,ak8sj53;
    TLorentzVector ak8sj14,ak8sj24,ak8sj34,ak8sj44,ak8sj54;
    TLorentzVector ak8sj15,ak8sj25,ak8sj35,ak8sj45,ak8sj55;

  double theWeight;
  double  nump=0;
  double  numm=0;
  double  npT, npIT;
  int     nBX;
  bool isHEEP;
  double iso, isoCut, et, trackIso;
  //Gen Level
  double gen_gra_m, gen_gra_pt, gen_gra_eta;
  double gen_rad_m, gen_rad_pt, gen_rad_eta;
  double gen_ele_pt, gen_ele_eta, gen_ele_phi, gen_ele_e;
  double gen_mu_pt, gen_mu_eta, gen_mu_phi, gen_mu_e;
  double gen_ele_pt_2, gen_ele_eta_2, gen_ele_phi_2, gen_ele_e_2;
  double gen_mu_pt_2, gen_mu_eta_2, gen_mu_phi_2, gen_mu_e_2;
  double gen_ele_pt_3, gen_ele_eta_3, gen_ele_phi_3, gen_ele_e_3;
  double gen_mu_pt_3, gen_mu_eta_3, gen_mu_phi_3, gen_mu_e_3;
  double gen_tau_pt, gen_tau_eta, gen_tau_phi, gen_tau_e;
  double gen_tau_pt_2, gen_tau_eta_2, gen_tau_phi_2, gen_tau_e_2;
  double gen_tau_pt_3, gen_tau_eta_3, gen_tau_phi_3, gen_tau_e_3;
  double gentop_pt, gentop_eta, gentop_phi, gentop_mass;
  double genantitop_pt, genantitop_eta, genantitop_phi, genantitop_mass;
  double ptGenVlep, etaGenVlep, phiGenVlep, massGenVlep;
  double ptGenVlep_2, etaGenVlep_2, phiGenVlep_2, massGenVlep_2;
  double ptGenVlep_3, etaGenVlep_3, phiGenVlep_3, massGenVlep_3;
  double ptGenVhad, etaGenVhad, phiGenVhad, massGenVhad;
  double ptGenV_2, etaGenV_2, phiGenV_2, massGenV_2;
  double ptGenV_3, etaGenV_3, phiGenV_3, massGenV_3;
  int status_1,status_2, status_3;

  double useless;

//  JEC
  double corr_AK8, corr_AK81[8];
  double corr_AK8puppi[8],corr_AK8puppiSD[8];
  double jetAK8_pt,jetAK8_mass,jetAK8_jec,jetAK8_e,jetAK8_eta,jetAK8_phi;
  double jetAK8_pt1[8], jetAK8_mass1[8], jetAK8_SF_mass1[8], jetAK8_SF_mass2[8], jetAK8_jec1[8],jetAK8_eta1[8];
  double jetAK8puppi_pt1[8], jetAK8puppi_mass1[8], jetAK8puppi_eta1[8], jetAK8puppi_jec1[8], jetAK8puppiSD_jec1[8];

  double corr;
  double METraw_et, METraw_phi, METraw_sumEt;
  double MET_et, MET_phi, MET_sumEt, MET_corrPx, MET_corrPy;
  // AK4 Jets
  int ak4jet_hf[8],ak4jet_pf[8],ak4jet_hf_2[8],ak4jet_pf_2[8];
  double ak4jet_pt[8],ak4jet_pt_uncorr[8],ak4jet_eta[8],ak4jet_phi[8],ak4jet_e[8], ak4jet_dr[8]; 
  double ak4jet_csv[8],ak4jet_icsv[8], ak4jet_IDLoose[8], ak4jet_IDTight[8];

    


  void setDummyValues();

  /// Parameters to steer the treeDumper
  int originalNEvents_;
  double crossSectionPb_;
  double targetLumiInvPb_;
  bool isGen_;
  bool isJEC_;
  bool RunOnMC_;
  std::vector<JetCorrectorParameters> vPar;
  std::map<std::string,double>  TypeICorrMap_;
  edm::InputTag mets_;

  //High Level Trigger
  HLTConfigProvider hltConfig;
  edm::EDGetTokenT<edm::TriggerResults> hltToken_;
  std::vector<std::string> muPaths1_, muPaths2_, muPaths3_, muPaths4_, muPaths5_, muPaths6_, muPaths7_, muPaths8_, muPaths9_;
  std::vector<std::string> muPaths1, muPaths2, muPaths3, muPaths4, muPaths5, muPaths6, muPaths7, muPaths8, muPaths9;
  std::vector<std::string> el1_, el2_, el3_, mu1_, mu2_, mu3_, mu4_;
  std::vector<std::string> el1, el2, el3, mu1, mu2, mu3, mu4;
  int  HLT_Mu1, HLT_Mu2, HLT_Mu3, HLT_Mu4, HLT_Mu5, HLT_Mu6, HLT_Mu7, HLT_Mu8, HLT_Mu9, HLT_Mu10, HLT_Mu11, HLT_Mu12, HLT_Mu13, HLT_Mu14, HLT_Mu15, HLT_Mu16;
  bool IDLoose, IDTight, IDLoose_2, IDTight_2, IDLoose_3, IDTight_3, IDLoose_4, IDTight_4;
          
// filter
  bool passFilter_HBHE_                   ;
  bool passFilter_HBHEIso_                ;
  bool passFilter_GlobalHalo_             ;
  bool passFilter_ECALDeadCell_           ;
  bool passFilter_GoodVtx_                ;
  bool passFilter_EEBadSc_                ;
  bool passFilter_badMuon_                ;
  bool passFilter_badChargedHadron_       ;

  edm::EDGetTokenT<edm::View<reco::Candidate>> leptonicVSrc_;
  edm::EDGetTokenT<edm::View<pat::Jet>> hadronicVSrc_;
  edm::EDGetTokenT<edm::View<pat::Jet>> ak4jetsSrc_;
  edm::EDGetTokenT<edm::View<pat::Electron> > looseelectronToken_ ;
  edm::EDGetTokenT<edm::View<pat::Muon>> loosemuonToken_;
  edm::EDGetTokenT<edm::View<pat::Electron> > vetoelectronToken_ ;
  edm::EDGetTokenT<edm::View<pat::Muon>> vetomuonToken_;
  edm::EDGetTokenT<edm::View<pat::Muon>> t1muSrc_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> gravitonSrc_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> metSrc_;
  edm::EDGetTokenT<GenEventInfoProduct> GenToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genSrc_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PUToken_;
  edm::EDGetTokenT<edm::View<pat::Jet>> jetsAK8Label_;
  edm::EDGetTokenT<LHEEventProduct> LheToken_;
};

//
// constructors and destructor
//
EDBRTreeMaker::EDBRTreeMaker(const edm::ParameterSet& iConfig):
  hltToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("hltToken"))),
  muPaths1_(iConfig.getParameter<std::vector<std::string>>("muPaths1")),
  muPaths2_(iConfig.getParameter<std::vector<std::string>>("muPaths2")),
  muPaths3_(iConfig.getParameter<std::vector<std::string>>("muPaths3")),
  muPaths4_(iConfig.getParameter<std::vector<std::string>>("muPaths4")),
  muPaths5_(iConfig.getParameter<std::vector<std::string>>("muPaths5")),
  muPaths6_(iConfig.getParameter<std::vector<std::string>>("muPaths6")),
  muPaths7_(iConfig.getParameter<std::vector<std::string>>("muPaths7")),
  muPaths8_(iConfig.getParameter<std::vector<std::string>>("muPaths8")),
  muPaths9_(iConfig.getParameter<std::vector<std::string>>("muPaths9")),
  el1_(iConfig.getParameter<std::vector<std::string>>("el1")),
  el2_(iConfig.getParameter<std::vector<std::string>>("el2")),
  el3_(iConfig.getParameter<std::vector<std::string>>("el3")),
  mu1_(iConfig.getParameter<std::vector<std::string>>("mu1")),
  mu2_(iConfig.getParameter<std::vector<std::string>>("mu2")),
  mu3_(iConfig.getParameter<std::vector<std::string>>("mu3")),
  mu4_(iConfig.getParameter<std::vector<std::string>>("mu4"))
{
  originalNEvents_ = iConfig.getParameter<int>("originalNEvents");
  crossSectionPb_  = iConfig.getParameter<double>("crossSectionPb");
  targetLumiInvPb_ = iConfig.getParameter<double>("targetLumiInvPb");
  isGen_           = iConfig.getParameter<bool>("isGen");
  isJEC_           = iConfig.getParameter<bool>("isJEC");
  RunOnMC_           = iConfig.getParameter<bool>("RunOnMC");
  leptonicVSrc_=consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>( "leptonicVSrc") ) ;
  jetsAK8Label_      = consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>( "ak8JetSrc") ) ;
  looseelectronToken_    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("looseElectronSrc"))) ;
  loosemuonToken_    = (consumes<edm::View<pat::Muon> > (iConfig.getParameter<edm::InputTag>("looseMuonSrc")));
  vetoelectronToken_    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("vetoElectronSrc"))) ;
  vetomuonToken_    = (consumes<edm::View<pat::Muon> > (iConfig.getParameter<edm::InputTag>("vetoMuonSrc")));
  LheToken_        = consumes<LHEEventProduct> (iConfig.getParameter<edm::InputTag>( "lhe") ) ;
  ak4jetsSrc_      = consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>( "ak4jetsSrc") ) ;
  hadronicVSrc_ = consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("hadronicVSrc") ) ;
  jetToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  puppijetToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("puppijets"));
  fatjetToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"));
  prunedjetToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("prunedjets"));
  softdropjetToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("softdropjets"));
  rhoToken_  = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
  vtxToken_  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  GenToken_=consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>( "generator") ) ;
  genSrc_      = consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>( "genSrc") ) ;
  PUToken_=consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileup") ) ;
  metSrc_      = consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>( "metSrc") ) ;
  gravitonSrc_      = consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>( "gravitonSrc") ) ;
  metToken_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  t1muSrc_      = consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>( "t1muSrc") ) ;
  noiseFilterToken_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("noiseFilter"));
  HBHENoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_HBHENoiseFilter");
  HBHENoiseIsoFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_HBHENoiseIsoFilter");
  GlobalHaloNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_GlobalTightHaloFilter");
  ECALDeadCellNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter");
  GoodVtxNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_goodVertices");
  EEBadScNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_eeBadScFilter");
  badMuonNoiseFilter_Selector_ = consumes<bool>(iConfig.getParameter<edm::InputTag> ("noiseFilterSelection_badMuon"));
  badChargedHadronNoiseFilter_Selector_  = consumes<bool>(iConfig.getParameter<edm::InputTag> ("noiseFilterSelection_badChargedHadron"));

  std::string jecpath = iConfig.getParameter<std::string>("jecpath");
  std::string tmpString;
  std::vector<std::string> tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK8chsPayloadNames");
  for( unsigned int v = 0; v < tmpVec.size(); ++v ){
     tmpString = jecpath + tmpVec[v];
     jecAK8Labels.push_back(tmpString);
  }
  std::vector<std::string> jecAK8LabelsGroomed;
  tmpVec.clear(); tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK8chsPayloadNamesGroomed");
  for( unsigned int v = 0; v < tmpVec.size(); ++v ){
     tmpString = jecpath + tmpVec[v];
     jecAK8LabelsGroomed.push_back(tmpString);
  }

  std::vector<std::string> jecAK8Labelspuppi;
  tmpVec.clear(); tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK8puppiPayloadNames");
  for( unsigned int v = 0; v < tmpVec.size(); ++v ){
     tmpString = jecpath + tmpVec[v];
     jecAK8Labelspuppi.push_back(tmpString);
  }

  std::vector<std::string> jecAK8LabelspuppiGroomed;
  tmpVec.clear(); tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK8puppiPayloadNamesGroomed");
  for( unsigned int v = 0; v < tmpVec.size(); ++v ){
     tmpString = jecpath + tmpVec[v];
     jecAK8LabelspuppiGroomed.push_back(tmpString);
  }



  std::vector<std::string> jecAK4Labels;
  tmpVec.clear(); tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK4chsPayloadNames");
  for( unsigned int v = 0; v < tmpVec.size(); ++v ){
     tmpString = jecpath + tmpVec[v];
     jecAK4Labels.push_back(tmpString);
  }

  /*=======================================================================================*/
  MW_=80.385;
  nmetmatch = 0;
  nmetno = 0;
  mettokens.push_back( metToken_ );
  mettokens.push_back( reclusteredmetToken_ );
  jetTokens.push_back( jetToken_ );
  jetTokens.push_back( fatjetToken_         );
  jetTokens.push_back( prunedjetToken_      );
  jetTokens.push_back( softdropjetToken_    );
  jetTokens.push_back( puppijetToken_      );

// add 3 up
 

  metInputToken_ = mettokens[0]; 
  reclusteredmetInputToken_ = mettokens[1];

  jetCorrLabel_ = jecAK4Labels;
  offsetCorrLabel_.push_back(jetCorrLabel_[0]);
 
  doCorrOnTheFly_ = false;
  if( jecAK4Labels.size() != 0 && jecAK8Labels.size() != 0 ){

     jecAK4PayloadNames_ = jecAK4Labels;
     //jecAK4PayloadNames_.pop_back();

     jecAK8PayloadNames_ = jecAK8Labels;
     //jecAK8PayloadNames_.pop_back();

     jecAK8PayloadNamesGroomed_ = jecAK8LabelsGroomed;
     //jecAK8PayloadNamesGroomed_.pop_back();

     jecAK8puppiPayloadNames_ = jecAK8Labelspuppi;
     jecAK8puppiPayloadNamesGroomed_ = jecAK8LabelspuppiGroomed;


  fatjetInputToken_ = jetTokens[1];
  prunedjetInputToken_ = jetTokens[2];
  softdropjetInputToken_ = jetTokens[3];
  puppijetInputToken_ = jetTokens[4];
// add 3 up
     initJetCorrFactors();

     doCorrOnTheFly_ = true;

  }

  
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;

    outTree_ = fs->make<TTree>("EDBRCandidates","EDBR Candidates");
    outTreew_ = fs->make<TTree>("EDBRCandidatesw","EDBR Candidates");
    outTree_->Branch("ak8sj11"    ,&ak8sj11 ) ;
    outTree_->Branch("ak8sj21"    ,&ak8sj21  );
    outTree_->Branch("ak8sj31"    ,&ak8sj31   );
    outTree_->Branch("ak8sj12"    ,&ak8sj12 ) ;
    outTree_->Branch("ak8sj22"    ,&ak8sj22  );
    outTree_->Branch("ak8sj32"    ,&ak8sj32   );
    outTree_->Branch("ak8sj13"    ,&ak8sj13 ) ;
    outTree_->Branch("ak8sj23"    ,&ak8sj23  );
    outTree_->Branch("ak8sj33"    ,&ak8sj33   );
    outTree_->Branch("ak8sj14"    ,&ak8sj14 ) ;
    outTree_->Branch("ak8sj24"    ,&ak8sj24  );
    outTree_->Branch("ak8sj34"    ,&ak8sj34   );
    outTree_->Branch("ak8sj15"    ,&ak8sj14 ) ;
    outTree_->Branch("ak8sj25"    ,&ak8sj24  );
    outTree_->Branch("ak8sj35"    ,&ak8sj34   );
  
/// Basic event quantities
  outTree_->Branch("run"             ,&run            ,"run/I");//
  outTree_->Branch("ls"              ,&ls             ,"ls/I"             );//Synch
  outTree_->Branch("event"           ,&nevent         ,"event/I"          );
  outTree_->Branch("nVtx"            ,&nVtx           ,"nVtx/I"           );
  outTree_->Branch("jetAK8puppi_ptJEC"          ,&jetAK8puppi_ptJEC         ,"jetAK8puppi_ptJEC/D"      );
  outTree_->Branch("jetAK8puppi_eta"          ,&jetAK8puppi_eta         ,"jetAK8puppi_eta/D"         );
  outTree_->Branch("jetAK8puppi_phi"          ,&jetAK8puppi_phi         ,"jetAK8puppi_phi/D"         );
  outTree_->Branch("jetAK8puppi_tau1"          ,&jetAK8puppi_tau1         ,"jetAK8puppi_tau1/D"         );
  outTree_->Branch("jetAK8puppi_tau2"          ,&jetAK8puppi_tau2         ,"jetAK8puppi_tau2/D"         );
  outTree_->Branch("jetAK8puppi_tau3"          ,&jetAK8puppi_tau3         ,"jetAK8puppi_tau3/D"         );
  outTree_->Branch("jetAK8puppi_tau21"          ,&jetAK8puppi_tau21         ,"jetAK8puppi_tau21/D"         );
  outTree_->Branch("jetAK8puppi_tau4"          ,&jetAK8puppi_tau4         ,"jetAK8puppi_tau4/D"         );
  outTree_->Branch("jetAK8puppi_tau42"          ,&jetAK8puppi_tau42         ,"jetAK8puppi_tau42/D"         );
  outTree_->Branch("jetAK8puppi_sd"          ,&jetAK8puppi_sd         ,"jetAK8puppi_sd/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC"          ,&jetAK8puppi_sdJEC         ,"jetAK8puppi_sdJEC/D"         );
  outTree_->Branch("jetAK8puppi_sdcorr"          ,&jetAK8puppi_sdcorr         ,"jetAK8puppi_sdcorr/D"         );

  outTree_->Branch("jetAK8puppi_ptJEC_2"          ,&jetAK8puppi_ptJEC_2         ,"jetAK8puppi_ptJEC_2/D"    );
  outTree_->Branch("jetAK8puppi_eta_2"          ,&jetAK8puppi_eta_2         ,"jetAK8puppi_eta_2/D"         );
  outTree_->Branch("jetAK8puppi_phi_2"          ,&jetAK8puppi_phi_2         ,"jetAK8puppi_phi_2/D"         );
  outTree_->Branch("jetAK8puppi_tau1_2"          ,&jetAK8puppi_tau1_2         ,"jetAK8puppi_tau1_2/D"         );
  outTree_->Branch("jetAK8puppi_tau2_2"          ,&jetAK8puppi_tau2_2         ,"jetAK8puppi_tau2_2/D"         );
  outTree_->Branch("jetAK8puppi_tau3_2"          ,&jetAK8puppi_tau3_2         ,"jetAK8puppi_tau3_2/D"         );
  outTree_->Branch("jetAK8puppi_tau21_2"          ,&jetAK8puppi_tau21_2         ,"jetAK8puppi_tau21_2/D"    );
  outTree_->Branch("jetAK8puppi_tau4_2"          ,&jetAK8puppi_tau4_2         ,"jetAK8puppi_tau4_2/D"       );
  outTree_->Branch("jetAK8puppi_tau42_2"          ,&jetAK8puppi_tau42_2         ,"jetAK8puppi_tau42_2/D"  );
  outTree_->Branch("jetAK8puppi_sd_2"          ,&jetAK8puppi_sd_2         ,"jetAK8puppi_sd_2/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC_2"          ,&jetAK8puppi_sdJEC_2         ,"jetAK8puppi_sdJEC_2/D"   );
  outTree_->Branch("jetAK8puppi_sdcorr_2"          ,&jetAK8puppi_sdcorr_2         ,"jetAK8puppi_sdcorr_2/D"  );

  outTree_->Branch("jetAK8puppi_ptJEC_3"          ,&jetAK8puppi_ptJEC_3         ,"jetAK8puppi_ptJEC_3/D"    );
  outTree_->Branch("jetAK8puppi_eta_3"          ,&jetAK8puppi_eta_3         ,"jetAK8puppi_eta_3/D"         );
  outTree_->Branch("jetAK8puppi_phi_3"          ,&jetAK8puppi_phi_3         ,"jetAK8puppi_phi_3/D"         );
  outTree_->Branch("jetAK8puppi_tau1_3"          ,&jetAK8puppi_tau1_3         ,"jetAK8puppi_tau1_3/D"         );
  outTree_->Branch("jetAK8puppi_tau2_3"          ,&jetAK8puppi_tau2_3         ,"jetAK8puppi_tau2_3/D"         );
  outTree_->Branch("jetAK8puppi_tau3_3"          ,&jetAK8puppi_tau3_3         ,"jetAK8puppi_tau3_3/D"         );
  outTree_->Branch("jetAK8puppi_tau21_3"          ,&jetAK8puppi_tau21_3         ,"jetAK8puppi_tau21_3/D"    );
  outTree_->Branch("jetAK8puppi_tau4_3"          ,&jetAK8puppi_tau4_3         ,"jetAK8puppi_tau4_3/D"       );
  outTree_->Branch("jetAK8puppi_tau42_3"          ,&jetAK8puppi_tau42_3         ,"jetAK8puppi_tau42_3/D"  );
  outTree_->Branch("jetAK8puppi_sd_3"          ,&jetAK8puppi_sd_3         ,"jetAK8puppi_sd_3/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC_3"          ,&jetAK8puppi_sdJEC_3         ,"jetAK8puppi_sdJEC_3/D"   );
  outTree_->Branch("jetAK8puppi_sdcorr_3"          ,&jetAK8puppi_sdcorr_3         ,"jetAK8puppi_sdcorr_3/D"  );

  outTree_->Branch("jetAK8puppi_ptJEC_4"          ,&jetAK8puppi_ptJEC_4         ,"jetAK8puppi_ptJEC_4/D"    );
  outTree_->Branch("jetAK8puppi_eta_4"          ,&jetAK8puppi_eta_4         ,"jetAK8puppi_eta_4/D"         );
  outTree_->Branch("jetAK8puppi_phi_4"          ,&jetAK8puppi_phi_4         ,"jetAK8puppi_phi_4/D"         );
  outTree_->Branch("jetAK8puppi_tau1_4"          ,&jetAK8puppi_tau1_4         ,"jetAK8puppi_tau1_4/D"         );
  outTree_->Branch("jetAK8puppi_tau2_4"          ,&jetAK8puppi_tau2_4         ,"jetAK8puppi_tau2_4/D"         );
  outTree_->Branch("jetAK8puppi_tau3_4"          ,&jetAK8puppi_tau3_4         ,"jetAK8puppi_tau3_4/D"         );
  outTree_->Branch("jetAK8puppi_tau21_4"          ,&jetAK8puppi_tau21_4         ,"jetAK8puppi_tau21_4/D"    );
  outTree_->Branch("jetAK8puppi_tau4_4"          ,&jetAK8puppi_tau4_4         ,"jetAK8puppi_tau4_4/D"       );
  outTree_->Branch("jetAK8puppi_tau42_4"          ,&jetAK8puppi_tau42_4         ,"jetAK8puppi_tau42_4/D"  );
  outTree_->Branch("jetAK8puppi_sd_4"          ,&jetAK8puppi_sd_4         ,"jetAK8puppi_sd_4/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC_4"          ,&jetAK8puppi_sdJEC_4         ,"jetAK8puppi_sdJEC_4/D"   );
  outTree_->Branch("jetAK8puppi_sdcorr_4"          ,&jetAK8puppi_sdcorr_4         ,"jetAK8puppi_sdcorr_4/D"  );


  outTree_->Branch("jetAK8puppi_ptJEC_5"          ,&jetAK8puppi_ptJEC_5         ,"jetAK8puppi_ptJEC_5/D"    );
  outTree_->Branch("jetAK8puppi_eta_5"          ,&jetAK8puppi_eta_5         ,"jetAK8puppi_eta_5/D"         );
  outTree_->Branch("jetAK8puppi_phi_5"          ,&jetAK8puppi_phi_5         ,"jetAK8puppi_phi_5/D"         );
  outTree_->Branch("jetAK8puppi_tau1_5"          ,&jetAK8puppi_tau1_5         ,"jetAK8puppi_tau1_5/D"         );
  outTree_->Branch("jetAK8puppi_tau2_5"          ,&jetAK8puppi_tau2_5         ,"jetAK8puppi_tau2_5/D"         );
  outTree_->Branch("jetAK8puppi_tau3_5"          ,&jetAK8puppi_tau3_5         ,"jetAK8puppi_tau3_5/D"         );
  outTree_->Branch("jetAK8puppi_tau21_5"          ,&jetAK8puppi_tau21_5         ,"jetAK8puppi_tau21_5/D"    );
  outTree_->Branch("jetAK8puppi_tau4_5"          ,&jetAK8puppi_tau4_5         ,"jetAK8puppi_tau4_5/D"       );
  outTree_->Branch("jetAK8puppi_tau42_5"          ,&jetAK8puppi_tau42_5         ,"jetAK8puppi_tau42_5/D"  );
  outTree_->Branch("jetAK8puppi_sd_5"          ,&jetAK8puppi_sd_5         ,"jetAK8puppi_sd_5/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC_5"          ,&jetAK8puppi_sdJEC_5         ,"jetAK8puppi_sdJEC_5/D"   );
  outTree_->Branch("jetAK8puppi_sdcorr_5"          ,&jetAK8puppi_sdcorr_5         ,"jetAK8puppi_sdcorr_5/D"  );


  outTree_->Branch("jetAK8puppi_ptJEC_6"          ,&jetAK8puppi_ptJEC_6         ,"jetAK8puppi_ptJEC_6/D"    );
  outTree_->Branch("jetAK8puppi_eta_6"          ,&jetAK8puppi_eta_6         ,"jetAK8puppi_eta_6/D"         );
  outTree_->Branch("jetAK8puppi_phi_6"          ,&jetAK8puppi_phi_6         ,"jetAK8puppi_phi_6/D"         );
  outTree_->Branch("jetAK8puppi_tau1_6"          ,&jetAK8puppi_tau1_6         ,"jetAK8puppi_tau1_6/D"         );
  outTree_->Branch("jetAK8puppi_tau2_6"          ,&jetAK8puppi_tau2_6         ,"jetAK8puppi_tau2_6/D"         );
  outTree_->Branch("jetAK8puppi_tau3_6"          ,&jetAK8puppi_tau3_6         ,"jetAK8puppi_tau3_6/D"         );
  outTree_->Branch("jetAK8puppi_tau21_6"          ,&jetAK8puppi_tau21_6         ,"jetAK8puppi_tau21_6/D"    );
  outTree_->Branch("jetAK8puppi_tau4_6"          ,&jetAK8puppi_tau4_6         ,"jetAK8puppi_tau4_6/D"       );
  outTree_->Branch("jetAK8puppi_tau42_6"          ,&jetAK8puppi_tau42_6         ,"jetAK8puppi_tau42_6/D"  );
  outTree_->Branch("jetAK8puppi_sd_6"          ,&jetAK8puppi_sd_6         ,"jetAK8puppi_sd_6/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC_6"          ,&jetAK8puppi_sdJEC_6         ,"jetAK8puppi_sdJEC_6/D"   );
  outTree_->Branch("jetAK8puppi_sdcorr_6"          ,&jetAK8puppi_sdcorr_6         ,"jetAK8puppi_sdcorr_6/D"  );


  outTree_->Branch("jetAK8puppi_ptJEC_7"          ,&jetAK8puppi_ptJEC_7         ,"jetAK8puppi_ptJEC_7/D"    );
  outTree_->Branch("jetAK8puppi_eta_7"          ,&jetAK8puppi_eta_7         ,"jetAK8puppi_eta_7/D"         );
  outTree_->Branch("jetAK8puppi_phi_7"          ,&jetAK8puppi_phi_7         ,"jetAK8puppi_phi_7/D"         );
  outTree_->Branch("jetAK8puppi_tau1_7"          ,&jetAK8puppi_tau1_7         ,"jetAK8puppi_tau1_7/D"         );
  outTree_->Branch("jetAK8puppi_tau2_7"          ,&jetAK8puppi_tau2_7         ,"jetAK8puppi_tau2_7/D"         );
  outTree_->Branch("jetAK8puppi_tau3_7"          ,&jetAK8puppi_tau3_7         ,"jetAK8puppi_tau3_7/D"         );
  outTree_->Branch("jetAK8puppi_tau21_7"          ,&jetAK8puppi_tau21_7         ,"jetAK8puppi_tau21_7/D"    );
  outTree_->Branch("jetAK8puppi_tau4_7"          ,&jetAK8puppi_tau4_7         ,"jetAK8puppi_tau4_7/D"       );
  outTree_->Branch("jetAK8puppi_tau42_7"          ,&jetAK8puppi_tau42_7         ,"jetAK8puppi_tau42_7/D"  );
  outTree_->Branch("jetAK8puppi_sd_7"          ,&jetAK8puppi_sd_7         ,"jetAK8puppi_sd_7/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC_7"          ,&jetAK8puppi_sdJEC_7         ,"jetAK8puppi_sdJEC_7/D"   );
  outTree_->Branch("jetAK8puppi_sdcorr_7"          ,&jetAK8puppi_sdcorr_7         ,"jetAK8puppi_sdcorr_7/D"  );


  outTree_->Branch("jetAK8puppi_ptJEC_8"          ,&jetAK8puppi_ptJEC_8         ,"jetAK8puppi_ptJEC_8/D"    );
  outTree_->Branch("jetAK8puppi_eta_8"          ,&jetAK8puppi_eta_8         ,"jetAK8puppi_eta_8/D"         );
  outTree_->Branch("jetAK8puppi_phi_8"          ,&jetAK8puppi_phi_8         ,"jetAK8puppi_phi_8/D"         );
  outTree_->Branch("jetAK8puppi_tau1_8"          ,&jetAK8puppi_tau1_8         ,"jetAK8puppi_tau1_8/D"         );
  outTree_->Branch("jetAK8puppi_tau2_8"          ,&jetAK8puppi_tau2_8         ,"jetAK8puppi_tau2_8/D"         );
  outTree_->Branch("jetAK8puppi_tau3_8"          ,&jetAK8puppi_tau3_8         ,"jetAK8puppi_tau3_8/D"         );
  outTree_->Branch("jetAK8puppi_tau21_8"          ,&jetAK8puppi_tau21_8         ,"jetAK8puppi_tau21_8/D"    );
  outTree_->Branch("jetAK8puppi_tau4_8"          ,&jetAK8puppi_tau4_8         ,"jetAK8puppi_tau4_8/D"       );
  outTree_->Branch("jetAK8puppi_tau42_8"          ,&jetAK8puppi_tau42_8         ,"jetAK8puppi_tau42_8/D"  );
  outTree_->Branch("jetAK8puppi_sd_8"          ,&jetAK8puppi_sd_8         ,"jetAK8puppi_sd_8/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC_8"          ,&jetAK8puppi_sdJEC_8         ,"jetAK8puppi_sdJEC_8/D"   );
  outTree_->Branch("jetAK8puppi_sdcorr_8"          ,&jetAK8puppi_sdcorr_8         ,"jetAK8puppi_sdcorr_8/D"  );

  /// Other quantities
  outTree_->Branch("theWeight", &theWeight, "theWeight/D");  
  outTreew_->Branch("theWeight", &theWeight, "theWeight/D");
  outTree_->Branch("nump", &nump, "nump/D");  
  outTree_->Branch("numm", &numm, "numm/D");  
  outTree_->Branch("npT"           ,&npT         ,"npT/D"          );
  outTree_->Branch("npIT"           ,&npIT         ,"npIT/D"          );
  outTree_->Branch("nBX"           ,&nBX         ,"nBX/I"          );
  outTree_->Branch("triggerWeight"   ,&triggerWeight  ,"triggerWeight/D"  );
  outTree_->Branch("lumiWeight"      ,&lumiWeight     ,"lumiWeight/D"     );
  outTree_->Branch("pileupWeight"    ,&pileupWeight   ,"pileupWeight/D"   );
//after JEC varible
  outTree_->Branch("met"             ,&met            ,"met/D"            );
  outTree_->Branch("metPhi"          ,&metPhi         ,"metPhi/D"         );
  outTree_->Branch("METraw_et",&METraw_et,"METraw_et/D");
  outTree_->Branch("METraw_phi",&METraw_phi,"METraw_phi/D");
  outTree_->Branch("METraw_sumEt",&METraw_sumEt,"METraw_sumEt/D");
  outTree_->Branch("MET_et",&MET_et,"MET_et/D");
  outTree_->Branch("MET_phi",&MET_phi,"MET_phi/D");
  outTree_->Branch("MET_sumEt",&MET_sumEt,"MET_sumEt/D");

  outTree_->Branch("jetAK8_pt",&jetAK8_pt,"jetAK8_pt/D");
  outTree_->Branch("jetAK8_mass",&jetAK8_mass,"jetAK8_mass/D");
  outTree_->Branch("jetAK8_jec",&jetAK8_jec,"jetAK8_jec/D");
  outTree_->Branch("jetAK8_pt1",&jetAK8_pt1,"jetAK8_pt1[3]/D");
  outTree_->Branch("jetAK8_eta1",&jetAK8_eta1,"jetAK8_eta1[3]/D");
  outTree_->Branch("jetAK8_mass1",&jetAK8_mass1,"jetAK8_mass1[3]/D");
  outTree_->Branch("jetAK8_SF_mass1",&jetAK8_SF_mass1,"jetAK8_SF_mass1[3]/D");
  outTree_->Branch("jetAK8_SF_mass2",&jetAK8_SF_mass2,"jetAK8_SF_mass2[3]/D");
  outTree_->Branch("jetAK8_jec1",&jetAK8_jec1,"jetAK8_jec1[3]/D");
  outTree_->Branch("jetAK8_eta",&jetAK8_eta,"jetAK8_eta/D");
  outTree_->Branch("jetAK8_phi",&jetAK8_phi,"jetAK8_phi/D");

  outTree_->Branch("candMasspuppiJEC",&candMasspuppiJEC,"candMasspuppiJEC/D");
  outTree_->Branch("massww",&massww,"massww[3]/D");

 
  ///HLT bits
  outTree_->Branch("HLT_Mu1"   ,&HLT_Mu1  ,"HLT_Mu1/I"  );
  outTree_->Branch("HLT_Mu2"   ,&HLT_Mu2  ,"HLT_Mu2/I"  );
  outTree_->Branch("HLT_Mu3"   ,&HLT_Mu3  ,"HLT_Mu3/I"  );
  outTree_->Branch("HLT_Mu4"   ,&HLT_Mu4  ,"HLT_Mu4/I"  );
  outTree_->Branch("HLT_Mu5"   ,&HLT_Mu5  ,"HLT_Mu5/I"  );
  outTree_->Branch("HLT_Mu6"   ,&HLT_Mu6  ,"HLT_Mu6/I"  );
  outTree_->Branch("HLT_Mu7"   ,&HLT_Mu7  ,"HLT_Mu7/I"  );
  outTree_->Branch("HLT_Mu8"   ,&HLT_Mu8  ,"HLT_Mu8/I"  );
  outTree_->Branch("HLT_Mu9"   ,&HLT_Mu9  ,"HLT_Mu9/I"  );
  outTree_->Branch("HLT_Mu10"   ,&HLT_Mu10  ,"HLT_Mu10/I"  );
  outTree_->Branch("HLT_Mu11"   ,&HLT_Mu11  ,"HLT_Mu11/I"  );
  outTree_->Branch("HLT_Mu12"   ,&HLT_Mu12  ,"HLT_Mu12/I"  );
  outTree_->Branch("HLT_Mu13"   ,&HLT_Mu13  ,"HLT_Mu13/I"  );
  outTree_->Branch("HLT_Mu14"   ,&HLT_Mu14  ,"HLT_Mu14/I"  );
  outTree_->Branch("HLT_Mu15"   ,&HLT_Mu15  ,"HLT_Mu15/I"  );
  outTree_->Branch("HLT_Mu16"   ,&HLT_Mu16  ,"HLT_Mu16/I"  );
// filter
  outTree_->Branch("passFilter_HBHE"                 ,&passFilter_HBHE_                ,"passFilter_HBHE_/O");
  outTree_->Branch("passFilter_HBHEIso"                 ,&passFilter_HBHEIso_                ,"passFilter_HBHEIso_/O");
  outTree_->Branch("passFilter_GlobalHalo"              ,&passFilter_GlobalHalo_             ,"passFilter_GlobalHalo_/O");
  outTree_->Branch("passFilter_ECALDeadCell"         ,&passFilter_ECALDeadCell_        ,"passFilter_ECALDeadCell_/O");
  outTree_->Branch("passFilter_GoodVtx"              ,&passFilter_GoodVtx_             ,"passFilter_GoodVtx_/O");
  outTree_->Branch("passFilter_EEBadSc"              ,&passFilter_EEBadSc_             ,"passFilter_EEBadSc_/O");
  outTree_->Branch("passFilter_badMuon"                 ,&passFilter_badMuon_                ,"passFilter_badMuon_/O");
  outTree_->Branch("passFilter_badChargedHadron"                 ,&passFilter_badChargedHadron_                ,"passFilter_badChargedHadron_/O");

  /// AK4 Jets Info
  outTree_->Branch("ak4jet_hf"        , ak4jet_hf       ,"ak4jet_hf[8]/I"       );
  outTree_->Branch("ak4jet_pf"        , ak4jet_pf       ,"ak4jet_pf[8]/I"       );
  outTree_->Branch("ak4jet_pt"        , ak4jet_pt       ,"ak4jet_pt[8]/D"       );
  outTree_->Branch("ak4jet_pt_uncorr"        , ak4jet_pt_uncorr       ,"ak4jet_pt_uncorr[8]/D"       );
  outTree_->Branch("ak4jet_eta"        , ak4jet_eta       ,"ak4jet_eta[8]/D"       );
  outTree_->Branch("ak4jet_phi"        , ak4jet_phi       ,"ak4jet_phi[8]/D"       );
  outTree_->Branch("ak4jet_e"        , ak4jet_e       ,"ak4jet_e[8]/D"       );
  outTree_->Branch("ak4jet_dr"        , ak4jet_dr       ,"ak4jet_dr[8]/D"       );
  outTree_->Branch("ak4jet_csv"        , ak4jet_csv       ,"ak4jet_csv[8]/D"       );
  outTree_->Branch("ak4jet_icsv"        , ak4jet_icsv       ,"ak4jet_icsv[8]/D"       );
  outTree_->Branch("ak4jet_IDLoose"        , ak4jet_IDLoose       ,"ak4jet_IDLoose[8]/D"       );
  outTree_->Branch("ak4jet_IDTight"        , ak4jet_IDTight       ,"ak4jet_IDTight[8]/D"       );

  /// Gen Level quantities
  outTree_->Branch("gen_gra_m"        ,&gen_gra_m       ,"gen_gra_m/D"       );
  outTree_->Branch("gen_gra_pt"        ,&gen_gra_pt       ,"gen_gra_pt/D"       );
  outTree_->Branch("gen_gra_eta"        ,&gen_gra_eta       ,"gen_gra_eta/D"       );
  outTree_->Branch("gen_rad_m"        ,&gen_rad_m       ,"gen_rad_m/D"       );
  outTree_->Branch("gen_rad_pt"        ,&gen_rad_pt       ,"gen_rad_pt/D"       );
  outTree_->Branch("gen_rad_eta"        ,&gen_rad_eta       ,"gen_rad_eta/D"       );
  outTree_->Branch("gen_tau_pt"        ,&gen_tau_pt       ,"gen_tau_pt/D"       );
  outTree_->Branch("gen_tau_eta"        ,&gen_tau_eta       ,"gen_tau_eta/D"       );
  outTree_->Branch("gen_tau_phi"        ,&gen_tau_phi       ,"gen_tau_phi/D"       );
  outTree_->Branch("gen_tau_e"        ,&gen_tau_e       ,"gen_tau_e/D"       );
  outTree_->Branch("gen_tau_pt_2"        ,&gen_tau_pt_2       ,"gen_tau_pt_2/D"       );
  outTree_->Branch("gen_tau_eta_2"        ,&gen_tau_eta_2       ,"gen_tau_eta_2/D"       );
  outTree_->Branch("gen_tau_phi_2"        ,&gen_tau_phi_2       ,"gen_tau_phi_2/D"       );
  outTree_->Branch("gen_tau_e_2"        ,&gen_tau_e_2       ,"gen_tau_e_2/D"       );
  outTree_->Branch("gen_tau_pt_3"        ,&gen_tau_pt_3       ,"gen_tau_pt_3/D"       );
  outTree_->Branch("gen_tau_eta_3"        ,&gen_tau_eta_3       ,"gen_tau_eta_3/D"       );
  outTree_->Branch("gen_tau_phi_3"        ,&gen_tau_phi_3       ,"gen_tau_phi_3/D"       );
  outTree_->Branch("gen_tau_e_3"        ,&gen_tau_e_3       ,"gen_tau_e_3/D"       );
 
  outTree_->Branch("gentop_pt"        ,&gentop_pt       ,"gentop_pt/D"       );
  outTree_->Branch("gentop_eta"        ,&gentop_eta       ,"gentop_eta/D"       );
  outTree_->Branch("gentop_phi"        ,&gentop_phi       ,"gentop_phi/D"       );
  outTree_->Branch("gentop_mass"        ,&gentop_mass       ,"gentop_mass/D"       );
  outTree_->Branch("genantitop_pt"        ,&genantitop_pt       ,"genantitop_pt/D"       );
  outTree_->Branch("genantitop_eta"        ,&genantitop_eta       ,"genantitop_eta/D"       );
  outTree_->Branch("genantitop_phi"        ,&genantitop_phi       ,"genantitop_phi/D"       );
  outTree_->Branch("genantitop_mass"        ,&genantitop_mass       ,"genantitop_mass/D"       );
  outTree_->Branch("ptq11"        ,&ptq11       ,"ptq11/D"       );
  outTree_->Branch("etaq11"        ,&etaq11       ,"etaq11/D"       );
  outTree_->Branch("phiq11"        ,&phiq11       ,"phiq11/D"       );
  outTree_->Branch("massq11"        ,&massq11       ,"massq11/D"       );
  outTree_->Branch("ptq21"        ,&ptq21       ,"ptq21/D"       );
  outTree_->Branch("etaq21"        ,&etaq21       ,"etaq21/D"       );
  outTree_->Branch("phiq21"        ,&phiq21       ,"phiq21/D"       );
  outTree_->Branch("massq21"        ,&massq21       ,"massq21/D"       );
  outTree_->Branch("ptq31"        ,&ptq31       ,"ptq31/D"       );
  outTree_->Branch("etaq31"        ,&etaq31       ,"etaq31/D"       );
  outTree_->Branch("phiq31"        ,&phiq31       ,"phiq31/D"       );
  outTree_->Branch("massq31"        ,&massq31       ,"massq31/D"       );

  outTree_->Branch("ptq12"        ,&ptq12       ,"ptq12/D"       );
  outTree_->Branch("etaq12"        ,&etaq12       ,"etaq12/D"       );
  outTree_->Branch("phiq12"        ,&phiq12       ,"phiq12/D"       );
  outTree_->Branch("massq12"        ,&massq12       ,"massq12/D"       );
  outTree_->Branch("ptq22"        ,&ptq22       ,"ptq22/D"       );
  outTree_->Branch("etaq22"        ,&etaq22       ,"etaq22/D"       );
  outTree_->Branch("phiq22"        ,&phiq22       ,"phiq22/D"       );
  outTree_->Branch("massq22"        ,&massq22       ,"massq22/D"       );
  outTree_->Branch("ptq32"        ,&ptq32       ,"ptq32/D"       );
  outTree_->Branch("etaq32"        ,&etaq32       ,"etaq32/D"       );
  outTree_->Branch("phiq32"        ,&phiq32       ,"phiq32/D"       );
  outTree_->Branch("massq32"        ,&massq32       ,"massq32/D"       );

  outTree_->Branch("ptGenVhad"        ,&ptGenVhad       ,"ptGenVhad/D"       );
  outTree_->Branch("etaGenVhad"        ,&etaGenVhad       ,"etaGenVhad/D"       );
  outTree_->Branch("phiGenVhad"        ,&phiGenVhad       ,"phiGenVhad/D"       );
  outTree_->Branch("massGenVhad"        ,&massGenVhad       ,"massGenVhad/D"       );
  outTree_->Branch("ptGenV_2"        ,&ptGenV_2       ,"ptGenV_2/D"       );
  outTree_->Branch("etaGenV_2"        ,&etaGenV_2       ,"etaGenV_2/D"       );
  outTree_->Branch("phiGenV_2"        ,&phiGenV_2       ,"phiGenV_2/D"       );
  outTree_->Branch("massGenV_2"        ,&massGenV_2       ,"massGenV_2/D"       );
  outTree_->Branch("ptGenV_3"        ,&ptGenV_3       ,"ptGenV_3/D"       );
  outTree_->Branch("etaGenV_3"        ,&etaGenV_3       ,"etaGenV_3/D"       );
  outTree_->Branch("phiGenV_3"        ,&phiGenV_3       ,"phiGenV_3/D"       );
  outTree_->Branch("massGenV_3"        ,&massGenV_3       ,"massGenV_3/D"       );

  outTree_->Branch("ptGenVlep"        ,&ptGenVlep       ,"ptGenVlep/D"       );
  outTree_->Branch("etaGenVlep"        ,&etaGenVlep       ,"etaGenVlep/D"       );
  outTree_->Branch("phiGenVlep"        ,&phiGenVlep       ,"phiGenVlep/D"       );
  outTree_->Branch("massGenVlep"        ,&massGenVlep       ,"massGenVlep/D"       );
  outTree_->Branch("ptGenVlep_2"        ,&ptGenVlep_2       ,"ptGenVlep_2/D"       );
  outTree_->Branch("etaGenVlep_2"        ,&etaGenVlep_2       ,"etaGenVlep_2/D"       );
  outTree_->Branch("phiGenVlep_2"        ,&phiGenVlep_2       ,"phiGenVlep_2/D"       );
  outTree_->Branch("massGenVlep_2"        ,&massGenVlep_2       ,"massGenVlep_2/D"       );
  outTree_->Branch("ptGenVlep_3"        ,&ptGenVlep_3       ,"ptGenVlep_3/D"       );
  outTree_->Branch("etaGenVlep_3"        ,&etaGenVlep_3       ,"etaGenVlep_3/D"       );
  outTree_->Branch("phiGenVlep_3"        ,&phiGenVlep_3       ,"phiGenVlep_3/D"       );
  outTree_->Branch("massGenVlep_3"        ,&massGenVlep_3       ,"massGenVlep_3/D"       );

  outTree_->Branch("status_1"           ,&status_1         ,"status_1/I"          );
  outTree_->Branch("status_2"           ,&status_2         ,"status_2/I"          );
  outTree_->Branch("status_3"           ,&status_3         ,"status_3/I"          );
  outTree_->Branch("pweight" ,pweight ,"pweight[211]/D" );
}

const reco::Candidate*  EDBRTreeMaker::findLastW(const reco::Candidate *particle,int IDpdg){
    int iw=0;
    int pidw=0;
    const reco::Candidate* pw=particle;
    for(int ii=0;particle->daughter(ii)!=NULL;ii++){
        if(abs(particle->daughter(ii)->pdgId())>pidw) {
            iw=ii;
            pidw=abs(particle->daughter(ii)->pdgId());
        }
    }
    if( abs(pidw) == IDpdg ){
        pw = particle->daughter(iw);
        return (findLastW(pw,IDpdg));
    }
    return pw;
    }
const reco::Candidate*  EDBRTreeMaker::findFirstW(const reco::Candidate *particle){
    if (particle->mother(0)!=NULL){
        if(abs(particle->mother(0)->pdgId()) == 24 )
              return (findFirstW(particle->mother(0)));
    }
    return particle;
}

EDBRTreeMaker::~EDBRTreeMaker()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


bool
EDBRTreeMaker::looseJetID( const pat::Jet& j ) {
	double NHF = j.neutralHadronEnergyFraction();
	double NEMF = j.neutralEmEnergyFraction();
	double CHF = j.chargedHadronEnergyFraction();
	double CEMF = j.chargedEmEnergyFraction();
	int NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
	int NumNeutralParticle =j.neutralMultiplicity();
	int CHM = j.chargedMultiplicity(); 
	double eta = j.eta();
          return ((  (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=2.7 ) || (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2 && abs(eta)>2.7 && abs(eta)<=3.0 ) || (NEMF<0.90 && NumNeutralParticle>10 && abs(eta)>3.0) ) ;

}

bool
EDBRTreeMaker::tightJetID( const pat::Jet& j ) {
// refer to https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
        double NHF = j.neutralHadronEnergyFraction();
        double NEMF = j.neutralEmEnergyFraction();
        double CHF = j.chargedHadronEnergyFraction();
        double CEMF = j.chargedEmEnergyFraction();
        int NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
        int NumNeutralParticle =j.neutralMultiplicity();
        int CHM = j.chargedMultiplicity();
        double eta = j.eta();
        return ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=3.0 ) || (NEMF<0.90 && NumNeutralParticle>10 && abs(eta)>3.0 )  ;
}

float
EDBRTreeMaker::dEtaInSeed( const pat::Electron*  ele ){
  return ele->superCluster().isNonnull() && ele->superCluster()->seed().isNonnull() ?
    ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta() : std::numeric_limits<float>::max();
}

void EDBRTreeMaker::initJetCorrFactors( void ){
   std::vector<JetCorrectorParameters> vPar;
   for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNames_.begin(), payloadEnd = jecAK8PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
   }
  // Make the FactorizedJetCorrector
  jecAK8_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNamesGroomed_.begin(), payloadEnd = jecAK8PayloadNamesGroomed_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
  // Make the FactorizedJetCorrector
  jecAK8Groomed_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNamesGroomed_.begin(), payloadEnd = jecAK8PayloadNamesGroomed_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
  jecAK8GroomedSD_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8puppiPayloadNames_.begin(), payloadEnd = jecAK8puppiPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
 jecAK8puppi_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8puppiPayloadNamesGroomed_.begin(), payloadEnd = jecAK8puppiPayloadNamesGroomed_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
 jecAK8puppiGroomed_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK4PayloadNames_.begin(), payloadEnd = jecAK4PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
  // Make the FactorizedJetCorrector
  jecAK4_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = offsetCorrLabel_.begin(), payloadEnd = offsetCorrLabel_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
  jecOffset_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
}


double EDBRTreeMaker::getJEC( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){
   double jetCorrFactor = 1.;
   if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
      jecAK4_->setJetEta( rawJetP4.eta() );
      jecAK4_->setJetPt ( rawJetP4.pt() );
      jecAK4_->setJetE  ( rawJetP4.energy() );
      jecAK4_->setJetPhi( rawJetP4.phi()    );
      jecAK4_->setJetA  ( jet.jetArea() );
      jecAK4_->setRho   ( *(rho_.product()) );
      jecAK4_->setNPV   ( nVtx );
      jetCorrFactor = jecAK4_->getCorrection();
   }
   reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
   corrJetP4 *= jetCorrFactor;
   return jetCorrFactor;
}

double EDBRTreeMaker::getJECOffset( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){
   double jetCorrFactor = 1.;
   if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
      jecOffset_->setJetEta( rawJetP4.eta()     );
      jecOffset_->setJetPt ( rawJetP4.pt()      );
      jecOffset_->setJetE  ( rawJetP4.energy()  );
      jecOffset_->setJetPhi( rawJetP4.phi()     );
      jecOffset_->setJetA  ( jet.jetArea()      );
      jecOffset_->setRho   ( *(rho_.product())  );
      jecOffset_->setNPV   ( nVtx  );
      jetCorrFactor = jecOffset_->getCorrection();
   }
   reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
   corrJetP4 *= jetCorrFactor;
   return jetCorrFactor;
}

//-------------------------------------------------------------------------------------------------------------------------------------//
//
// member functions
//
void EDBRTreeMaker::addTypeICorr( edm::Event const & event ){
   TypeICorrMap_.clear();

   event.getByToken(jetToken_      , jets_    );
   event.getByToken(rhoToken_      , rho_     );
   edm::Handle<reco::VertexCollection> vertices_;
   event.getByToken(vtxToken_, vertices_);
   edm::Handle<edm::View<pat::Muon>> muons_;
   event.getByToken(t1muSrc_,muons_);
   bool skipEM_                    = true;
   double skipEMfractionThreshold_ = 0.9;
   bool skipMuons_                 = true;
   std::string skipMuonSelection_string = "isGlobalMuon | isStandAloneMuon";
   StringCutObjectSelector<reco::Candidate>* skipMuonSelection_ = new StringCutObjectSelector<reco::Candidate>(skipMuonSelection_string,true);
   double jetCorrEtaMax_           = 9.9;
   double type1JetPtThreshold_     = 15.0; //10.0;
   double corrEx    = 0;
   double corrEy    = 0;
   double corrSumEt = 0;
   for (const pat::Jet &jet : *jets_) {
     double emEnergyFraction = jet.chargedEmEnergyFraction() + jet.neutralEmEnergyFraction();
     if ( skipEM_ && emEnergyFraction > skipEMfractionThreshold_ ) continue;
     reco::Candidate::LorentzVector rawJetP4 = jet.correctedP4(0);
     double corr = getJEC(rawJetP4, jet, jetCorrEtaMax_, jetCorrLabel_);
     if ( skipMuons_ ) {
          const std::vector<reco::CandidatePtr> & cands = jet.daughterPtrVector();
          for ( std::vector<reco::CandidatePtr>::const_iterator cand = cands.begin();
             cand != cands.end(); ++cand ) {
     	     const reco::PFCandidate *pfcand = dynamic_cast<const reco::PFCandidate *>(cand->get());
     	     const reco::Candidate *mu = (pfcand != 0 ? ( pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand->get());
               if ( mu != 0 && (*skipMuonSelection_)(*mu) ) {
                  reco::Candidate::LorentzVector muonP4 = (*cand)->p4();
                  rawJetP4 -= muonP4;
               }
           }
      }
     reco::Candidate::LorentzVector corrJetP4 = corr*rawJetP4;
     if ( corrJetP4.pt() > type1JetPtThreshold_ ) {
                 reco::Candidate::LorentzVector tmpP4 = jet.correctedP4(0);
                 corr = getJECOffset(tmpP4, jet, jetCorrEtaMax_, offsetCorrLabel_);
                 reco::Candidate::LorentzVector rawJetP4offsetCorr = corr*rawJetP4;
                 corrEx    -= (corrJetP4.px() - rawJetP4offsetCorr.px());
                 corrEy    -= (corrJetP4.py() - rawJetP4offsetCorr.py());
                 corrSumEt += (corrJetP4.Et() - rawJetP4offsetCorr.Et());
      }
  }
 TypeICorrMap_["corrEx"]    = corrEx;
 TypeICorrMap_["corrEy"]    = corrEy;
 TypeICorrMap_["corrSumEt"] = corrSumEt;
}

//-------------------------------------------------------------------------------------------------------------------------------------//
math::XYZTLorentzVector
EDBRTreeMaker::getNeutrinoP4(double& MetPt, double& MetPhi, TLorentzVector& lep, int lepType){
    double leppt = lep.Pt();
    double lepphi = lep.Phi();
    double lepeta = lep.Eta();
    double lepenergy = lep.Energy();
    double metpt = MetPt;
    double metphi = MetPhi;
    double  px = metpt*cos(metphi);
    double  py = metpt*sin(metphi);
    double  pz = 0;
    double  pxl= leppt*cos(lepphi);
    double  pyl= leppt*sin(lepphi);
    double  pzl= leppt*sinh(lepeta);
    double  El = lepenergy;
    double  a = pow(MW_,2) + pow(px+pxl,2) + pow(py+pyl,2) - px*px - py*py - El*El + pzl*pzl;
    double  b = 2.*pzl;
    double  A = b*b -4.*El*El;
    double  B = 2.*a*b;
    double  C = a*a-4.*(px*px+py*py)*El*El;
    double M_mu =  0;
    int type=2; // use the small abs real root
    a = MW_*MW_ - M_mu*M_mu + 2.0*pxl*px + 2.0*pyl*py;
    A = 4.0*(El*El - pzl*pzl);
    B = -4.0*a*pzl;
    C = 4.0*El*El*(px*px + py*py) - a*a;
    double tmproot = B*B - 4.0*A*C;
    if (tmproot<0) {
        pz = - B/(2*A); // take real part of complex roots
    }
    else {
        double tmpsol1 = (-B + sqrt(tmproot))/(2.0*A);
        double tmpsol2 = (-B - sqrt(tmproot))/(2.0*A);
        if (type == 0 ) {
            if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
            else { pz = tmpsol1; }
            if ( abs(pz) > 300. ) {
                if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
                else { pz = tmpsol2; }
            }
        }
        if (type == 1 ) {
            // two real roots, pick the one closest to pz of muon
            if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
            else {pz = tmpsol1; }
        }
        if (type == 2 ) {
            // pick the most central root.
            if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
            else { pz = tmpsol2; }
        }
    }//endl of if real root
    //dont correct pt neutrino
    math::XYZTLorentzVector outP4(px,py,pz,sqrt(px*px+py*py+pz*pz));
    return outP4;
}//end neutrinoP4


// ------------ method called for each event  ------------
void
EDBRTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   setDummyValues(); //Initalize variables with dummy values
   nevent = iEvent.eventAuxiliary().event();
   run    = iEvent.eventAuxiliary().run();
   ls     = iEvent.eventAuxiliary().luminosityBlock();
   Handle<TriggerResults> trigRes;
   iEvent.getByToken(hltToken_, trigRes);
 

   int mtemp1=0;
   for (size_t i=0; i<muPaths1.size();i++) {
      mtemp1 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths1[i]));
      if(HLT_Mu1<mtemp1) HLT_Mu1=mtemp1;
   }
   int mtemp2=0;
   for (size_t i=0; i<muPaths2.size();i++) {
      mtemp2 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths2[i]));
      if(HLT_Mu2<mtemp2) HLT_Mu2=mtemp2;
   }
   int mtemp3=0;
   for (size_t i=0; i<muPaths3.size();i++) {
      mtemp3 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths3[i]));
      if(HLT_Mu3<mtemp3) HLT_Mu3=mtemp3;
   }
   int mtemp4=0;
   for (size_t i=0; i<muPaths4.size();i++) {
      mtemp4 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths4[i]));
      if(HLT_Mu4<mtemp4) HLT_Mu4=mtemp4;
   }
   int mtemp5=0;
   for (size_t i=0; i<muPaths5.size();i++) {
      mtemp5 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths5[i]));
      if(HLT_Mu5<mtemp5) HLT_Mu5=mtemp5;
   }
   int mtemp6=0;
   for (size_t i=0; i<muPaths6.size();i++) {
      mtemp6 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths6[i]));
      if(HLT_Mu6<mtemp6) HLT_Mu6=mtemp6;
   }
   int mtemp7=0;
   for (size_t i=0; i<muPaths7.size();i++) {
      mtemp7 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths7[i]));
      if(HLT_Mu7<mtemp7) HLT_Mu7=mtemp7;
   }
   int mtemp8=0;
   for (size_t i=0; i<muPaths8.size();i++) {
      mtemp8 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths8[i]));
      if(HLT_Mu8<mtemp8) HLT_Mu8=mtemp8;
   }
   int mtemp9=0;
   for (size_t i=0; i<muPaths9.size();i++) {
      mtemp9 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths9[i]));
      if(HLT_Mu9<mtemp9) HLT_Mu9=mtemp9;
   }
   int mtemp10=0;
   for (size_t i=0; i<el1.size();i++) {
      mtemp10 = (int)trigRes->accept(hltConfig.triggerIndex(el1[i]));
      if(HLT_Mu10<mtemp10) HLT_Mu10=mtemp10;
   }
   int mtemp11=0;
   for (size_t i=0; i<el2.size();i++) {
      mtemp11 = (int)trigRes->accept(hltConfig.triggerIndex(el2[i]));
      if(HLT_Mu11<mtemp11) HLT_Mu11=mtemp11;
   }

   int mtemp12=0;
   for (size_t i=0; i<el3.size();i++) {
      mtemp12 = (int)trigRes->accept(hltConfig.triggerIndex(el3[i]));
      if(HLT_Mu12<mtemp12) HLT_Mu12=mtemp12;
   }

   int mtemp13=0;
   for (size_t i=0; i<mu1.size();i++) {
      mtemp13 = (int)trigRes->accept(hltConfig.triggerIndex(mu1[i]));
      if(HLT_Mu13<mtemp13) HLT_Mu13=mtemp13;
   }

   int mtemp14=0;
   for (size_t i=0; i<mu2.size();i++) {
      mtemp14 = (int)trigRes->accept(hltConfig.triggerIndex(mu2[i]));
      if(HLT_Mu14<mtemp14) HLT_Mu14=mtemp14;
   }  

   int mtemp15=0;
   for (size_t i=0; i<mu3.size();i++) {
      mtemp15 = (int)trigRes->accept(hltConfig.triggerIndex(mu3[i]));
      if(HLT_Mu15<mtemp15) HLT_Mu15=mtemp15;
   }  

   int mtemp16=0;
   for (size_t i=0; i<mu4.size();i++) {
      mtemp16 = (int)trigRes->accept(hltConfig.triggerIndex(mu4[i]));
      if(HLT_Mu16<mtemp16) HLT_Mu16=mtemp16;
   }  

   edm::Handle<edm::View<pat::Jet> > hadronicVs;

   iEvent.getByToken(hadronicVSrc_, hadronicVs);
   edm::Handle<edm::View<reco::Candidate> > leptonicVs;
   iEvent.getByToken(leptonicVSrc_, leptonicVs);
   edm::Handle<double> rho;

   iEvent.getByToken(rhoToken_      , rho     );
   double fastJetRho = *(rho.product());
   useless = fastJetRho;
   edm::Handle<edm::View<pat::Jet> > ak4jets;
   iEvent.getByToken(ak4jetsSrc_, ak4jets);
   edm::Handle<edm::View<reco::Candidate> > gravitons;

//   iEvent.getByLabel(gravitonSrc_.c_str(), gravitons);
   iEvent.getByToken(gravitonSrc_, gravitons);
   edm::Handle<edm::View<reco::Candidate> > metHandle;
   iEvent.getByToken(metSrc_, metHandle);
   edm::Handle<edm::View<pat::Muon>> loosemus;
   iEvent.getByToken(loosemuonToken_,loosemus);

   edm::Handle<edm::View<pat::Muon>> vetomus;
   iEvent.getByToken(vetomuonToken_,vetomus);
   edm::Handle<edm::View<pat::Electron>> looseels;

   iEvent.getByToken(looseelectronToken_, looseels);
   edm::Handle<edm::View<pat::Electron>> vetoels;

   iEvent.getByToken(vetoelectronToken_, vetoels);
   edm::Handle<edm::View<reco::GenParticle> > genParticles;//define genParticle
   iEvent.getByToken(genSrc_, genParticles);

   if (RunOnMC_){
     edm::Handle<LHEEventProduct> wgtsource;
     iEvent.getByToken(LheToken_, wgtsource);
     for ( int i=0; i<211; ++i) {
       pweight[i]= wgtsource->weights()[i].wgt/wgtsource->originalXWGTUP();
                              }
     edm::Handle<GenEventInfoProduct> genEvtInfo;
     iEvent.getByToken(GenToken_,genEvtInfo);
     theWeight = genEvtInfo->weight();
     if(theWeight>0) nump = nump+1;
     if(theWeight<0) numm = numm+1;

     edm::Handle<std::vector<PileupSummaryInfo>>  PupInfo;
     iEvent.getByToken(PUToken_, PupInfo);
     std::vector<PileupSummaryInfo>::const_iterator PVI;
     for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
         nBX = PVI->getBunchCrossing();
         if(nBX == 0) { 
            npT = PVI->getTrueNumInteractions();
            npIT = PVI->getPU_NumInteractions();
         }
      }
   }
//filter
    iEvent.getByToken(noiseFilterToken_, noiseFilterBits_);
    const edm::TriggerNames &names = iEvent.triggerNames(*noiseFilterBits_);
    for (unsigned int i = 0, n = noiseFilterBits_->size(); i < n; ++i) {
      if (names.triggerName(i) == HBHENoiseFilter_Selector_)
        passFilter_HBHE_ = noiseFilterBits_->accept(i); // TO BE USED
      if (names.triggerName(i) == HBHENoiseIsoFilter_Selector_)
        passFilter_HBHEIso_ = noiseFilterBits_->accept(i); // TO BE USED
      if (names.triggerName(i) == GlobalHaloNoiseFilter_Selector_)
        passFilter_GlobalHalo_ = noiseFilterBits_->accept(i); // TO BE USED
      if (names.triggerName(i) == ECALDeadCellNoiseFilter_Selector_)
        passFilter_ECALDeadCell_ = noiseFilterBits_->accept(i); // under scrutiny
      if (names.triggerName(i) == GoodVtxNoiseFilter_Selector_)
        passFilter_GoodVtx_ = noiseFilterBits_->accept(i); // TO BE USED
      if (names.triggerName(i) == EEBadScNoiseFilter_Selector_)
        passFilter_EEBadSc_ = noiseFilterBits_->accept(i); // under scrutiny  
    }
     edm::Handle<bool> badMuonResultHandle;
     edm::Handle<bool> badChargedHadronResultHandle;
     iEvent.getByToken(badMuonNoiseFilter_Selector_, badMuonResultHandle);
     iEvent.getByToken(badChargedHadronNoiseFilter_Selector_, badChargedHadronResultHandle);
     passFilter_badMuon_ = *badMuonResultHandle;
     passFilter_badChargedHadron_ = *badChargedHadronResultHandle;

   numCands = gravitons->size();
 
// ************************* Gen Level Information******************//
   if(RunOnMC_)
   {//MC Info
       Int_t havegra=0;
	  for(size_t ik=0; ik<genParticles->size();ik++)
	{// loop on gen
              const reco::Candidate* ptop = &(*genParticles)[ik];
                if(ptop->pdgId()== 6) {
                   gentop_pt = ptop->pt();
                   gentop_eta = ptop->eta();
                   gentop_phi = ptop->phi();
                   gentop_mass = ptop->mass(); }
                if(ptop->pdgId()== -6) {
                   genantitop_pt = ptop->pt();
                   genantitop_eta = ptop->eta();
                   genantitop_phi = ptop->phi();
                   genantitop_mass = ptop->mass();}
        if( abs((*genParticles)[ik].pdgId())==9000024) havegra=1;
		if( abs((*genParticles)[ik].pdgId())==9000024 || abs((*genParticles)[ik].pdgId())==6 )
		{//if Wkk
            const reco::Candidate* pwkk0 = &(*genParticles)[ik];
             const reco::Candidate* pwkk=findLastW(pwkk0,9000024);
                   gen_gra_m=pwkk->mass();
		   gen_gra_pt=pwkk->pt();
		   gen_gra_eta=pwkk->eta();
                     qnumber=1;
                for(int i=0;pwkk->daughter(i)!=NULL;i++)
                   {//loop on Wkk daughter
                    if(abs(pwkk->daughter(i)->pdgId())==24)
                       {//if w
                         const reco::Candidate* pw0 = pwkk->daughter(i);
                           const reco::Candidate* pw= findLastW(pw0,24);  
                         if(pw->daughter(0)!=NULL)
                         {//loop on w daughter
                            const reco::Candidate* pl = pw->daughter(0);
                            if( (abs(pl->pdgId())==11) || (abs(pl->pdgId())==13)|| (abs(pl->pdgId())==15)||(abs(pl->pdgId())==12) || (abs(pl->pdgId())==14)|| (abs(pl->pdgId())==16))
                            {//beign of lep-w
                                ptGenVlep = pw->pt();
                                etaGenVlep = pw->eta();
                                phiGenVlep = pw->phi();
                                massGenVlep = pw->mass();
                                status_1=0;
                            for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                                const reco::Candidate* pl = pw->daughter(ii);
                               if(abs(pl->pdgId())==11)
                               {
                                 gen_ele_pt=pl->pt();
                                 gen_ele_eta=pl->eta();
                                 gen_ele_phi=pl->phi();
                                 gen_ele_e=pl->energy();
                                   status_1=1;
                               }
                               if(abs(pl->pdgId())==13)
                               {
                                 gen_mu_pt=pl->pt();
                                 gen_mu_eta=pl->eta();
                                 gen_mu_phi=pl->phi();
                                 gen_mu_e=pl->energy();
                                   status_1=2;
                               }
                                if(abs(pl->pdgId())==15)
                                    {
                            gen_tau_pt=pl->pt();
                            gen_tau_eta=pl->eta();
                            gen_tau_phi=pl->phi();
                            gen_tau_e=pl->energy();
                                        status_1=3;
                                    }
                                }
                             }//end of if lep-w

                             for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                                const reco::Candidate* pl = pw->daughter(ii);
                               if(abs(pl->pdgId())<6)
                               {
                              if(qnumber==1)
                                 {
                                 ptq11 = pl->pt();
                                 etaq11 = pl->eta();
                                 phiq11 = pl->phi();
                                 massq11 = pl->mass();
                                 }
                                 if(qnumber==2)
                                 {
                                 ptq12 = pl->pt();
                                 etaq12 = pl->eta();
                                 phiq12 = pl->phi();
                                 massq12 = pl->mass();
                                 }
                                 qnumber=qnumber+1;
                               }                                          
                                                                     }
                             pl = pw->daughter(0);

        		     if(abs(pl->pdgId())<6) 
                             {
			         ptGenVhad = pw->pt();
			         etaGenVhad = pw->eta();
			         phiGenVhad = pw->phi();
			         massGenVhad = pw->mass();
                                 status_1=4;
                             }
        		     if(abs(pl->pdgId())==24) 
                             {
                                 status_1=5;
                             }
			   }//end of loop on w daughter
                       }//end of if w				 
                   }//end of loop on Wkk daughter
 		}//end of if Wkk
        if( havegra==0&&abs((*genParticles)[ik].pdgId())==24 )
        {//if W
             qnumber=1;
            const reco::Candidate* pw0 = &(*genParticles)[ik];
            const reco::Candidate* pw=findFirstW(pw0);
            if(pw->mother(0)->pdgId()!=9000025){
                const reco::Candidate* pw1= findLastW(pw0,24);
                if(pw1->daughter(0)!=NULL)
                {//loop on w daughter
                const reco::Candidate* pl = pw1->daughter(0);
                if( (abs(pl->pdgId())==11) || (abs(pl->pdgId())==13)|| (abs(pl->pdgId())==15)||(abs(pl->pdgId())==12) || (abs(pl->pdgId())==14)|| (abs(pl->pdgId())==16))
                {//beign of lep-w
                    ptGenVlep = pw1->pt();
                    etaGenVlep = pw1->eta();
                    phiGenVlep = pw1->phi();
                    massGenVlep = pw1->mass();
                    status_1=0;
                    for(int ii=0;pw1->daughter(ii)!=NULL;ii++){
                        const reco::Candidate* pl = pw1->daughter(ii);
                        if(abs(pl->pdgId())==11)
                        {
                            gen_ele_pt=pl->pt();
                            gen_ele_eta=pl->eta();
                            gen_ele_phi=pl->phi();
                            gen_ele_e=pl->energy();
                            status_1=1;
                        }
                        if(abs(pl->pdgId())==13)
                        {
                            gen_mu_pt=pl->pt();
                            gen_mu_eta=pl->eta();
                            gen_mu_phi=pl->phi();
                            gen_mu_e=pl->energy();
                            status_1=2;
                        }
                        if(abs(pl->pdgId())==15)
                        {
                            gen_tau_pt=pl->pt();
                            gen_tau_eta=pl->eta();
                            gen_tau_phi=pl->phi();
                            gen_tau_e=pl->energy();
                            status_1=3;
                        }
                    }
                }//end of if lep-w

                             for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                                const reco::Candidate* pl = pw->daughter(ii);
                               if(abs(pl->pdgId())<6)
                               {
                              if(qnumber==1)
                                 {
                                 ptq11 = pl->pt();
                                 etaq11 = pl->eta();
                                 phiq11 = pl->phi();
                                 massq11 = pl->mass();
                                 }
                                 if(qnumber==2)
                                 {
                                 ptq12 = pl->pt();
                                 etaq12 = pl->eta();
                                 phiq12 = pl->phi();
                                 massq12 = pl->mass();
                                 }
                                 qnumber=qnumber+1;
                               }
                                                                     }
                             pl = pw->daughter(0);


                if(abs(pl->pdgId())<6)
                {
                    ptGenVhad = pw1->pt();
                    etaGenVhad = pw1->eta();
                    phiGenVhad = pw1->phi();
                    massGenVhad = pw1->mass();
                    status_1=4;
                }
                if(abs(pl->pdgId())==24)
                {
                    status_1=5;
                }
                }
            }
        }
		if( abs((*genParticles)[ik].pdgId())==9000025 )
		{//if Radion
                   gen_rad_m=(*genParticles)[ik].mass();
		   gen_rad_pt=(*genParticles)[ik].pt();
	           gen_rad_eta=(*genParticles)[ik].eta();
                   qnumber2=1;
                   qnumber3=1;
                   for(int i=0;(*genParticles)[ik].daughter(i)!=NULL;i++)
                   {//loop on Radion daughter
                                 
                      if(((*genParticles)[ik].daughter(i)->pdgId())==24)
                       {//if w-
                         const reco::Candidate* pw0 = (*genParticles)[ik].daughter(i);
                           const reco::Candidate* pw= findLastW(pw0,24);
                         if(pw->daughter(0)!=NULL)
                         {//loop on w daughter
                            const reco::Candidate* pl = pw->daughter(0);
                            if( (abs(pl->pdgId())==11) || (abs(pl->pdgId())==13)||(abs(pl->pdgId())==15)||(abs(pl->pdgId())==12) || (abs(pl->pdgId())==14)|| (abs(pl->pdgId())==16))
                            {//beign of lep-w
                               ptGenVlep_2 = pw->pt();
                               etaGenVlep_2 = pw->eta();
                               phiGenVlep_2 = pw->phi();
                               massGenVlep_2 = pw->mass();
                               status_2=0;
                            for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                               const reco::Candidate* pl = pw->daughter(ii);
                               if(abs(pl->pdgId())==11)
                               {
                                 gen_ele_pt_2=pl->pt();
                                 gen_ele_eta_2=pl->eta();
                                 gen_ele_phi_2=pl->phi();
                                 gen_ele_e_2=pl->energy();
                                 status_2=1;  
                               }
                               if(abs(pl->pdgId())==13)
                               {
                                 gen_mu_pt_2=pl->pt();
                                 gen_mu_eta_2=pl->eta();
                                 gen_mu_phi_2=pl->phi();
                                 gen_mu_e_2=pl->energy();
                                 status_2=2;  
                               }
				if(abs(pl->pdgId())==15)
                               {
                                 gen_tau_pt_2=pl->pt();
                                 gen_tau_eta_2=pl->eta();
                                 gen_tau_phi_2=pl->phi();
                                 gen_tau_e_2=pl->energy();
                                 status_2=3;
                               }
                            }
                             }//end of if lep-w

                             for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                                const reco::Candidate* pl = pw->daughter(ii);
                               if(abs(pl->pdgId())<6)
                               {
                                if(qnumber2==1)
                                 {
                                 ptq21 = pl->pt();
                                 etaq21 = pl->eta();
                                 phiq21 = pl->phi();
                                 massq21 = pl->mass();
                                 }
                                 if(qnumber2==2)
                                 {
                                 ptq22 = pl->pt();
                                 etaq22 = pl->eta();
                                 phiq22 = pl->phi();
                                 massq22 = pl->mass();
                                 }
                                 qnumber2=qnumber2+1;                               
                               }
                                                                     }
                             pl = pw->daughter(0);

        		     if(abs(pl->pdgId())<6) 
                             {
                                 ptGenV_2 = pw->pt();
                                 etaGenV_2 = pw->eta();
                                 phiGenV_2 = pw->phi();
                                 massGenV_2 = pw->mass();
                                 status_2=4;  
                             }
        		     if(abs(pl->pdgId())==24) 
                             {
                                 status_2=5;  
                             }
			   }//end of loop on w daughter
                       }//end of if w-	
                      if(((*genParticles)[ik].daughter(i)->pdgId())==-24)
                       {//if w+
                         const reco::Candidate* pw0 = (*genParticles)[ik].daughter(i);
                           const reco::Candidate* pw= findLastW(pw0,24);
                         if(pw->daughter(0)!=NULL)
                         {
                            const reco::Candidate* pl = pw->daughter(0);
                            if( (abs(pl->pdgId())==11) || (abs(pl->pdgId())==13)|| (abs(pl->pdgId())==15)||(abs(pl->pdgId())==12) || (abs(pl->pdgId())==14)|| (abs(pl->pdgId())==16))
                            {//beign of lep-w
                               ptGenVlep_3 = pw->pt();
                               etaGenVlep_3 = pw->eta();
                               phiGenVlep_3 = pw->phi();
                               massGenVlep_3 = pw->mass();
                               status_3=0;
                            for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                              const reco::Candidate* pl = pw->daughter(ii);
                               if(abs(pl->pdgId())==11)
                               {
                                 gen_ele_pt_3=pl->pt();
                                 gen_ele_eta_3=pl->eta();
                                 gen_ele_phi_3=pl->phi();
                                 gen_ele_e_3=pl->energy();
                                 status_3=1;  
                               }
                               if(abs(pl->pdgId())==13)
                               {
                                 gen_mu_pt_3=pl->pt();
                                 gen_mu_eta_3=pl->eta();
                                 gen_mu_phi_3=pl->phi();
                                 gen_mu_e_3=pl->energy();
                                 status_3=2;  
                               }
                               if(abs(pl->pdgId())==15)
                               {
                                 gen_tau_pt_3=pl->pt();
                                 gen_tau_eta_3=pl->eta();
                                 gen_tau_phi_3=pl->phi();
                                 gen_tau_e_3=pl->energy();
                                 status_3=3;  
                               }
                            }
                             }//end of if lep-w

                             for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                                const reco::Candidate* pl = pw->daughter(ii);
                               if(abs(pl->pdgId())<6)
                               {
                                if(qnumber3==1)
                                 {
                                 ptq31 = pl->pt();
                                 etaq31 = pl->eta();
                                 phiq31 = pl->phi();
                                 massq31 = pl->mass();
                                 }
                                 if(qnumber3==2)
                                 {
                                 ptq32 = pl->pt();
                                 etaq32 = pl->eta();
                                 phiq32 = pl->phi();
                                 massq32 = pl->mass();
                                 }
                                 qnumber3=qnumber3+1;
                               }
                                                                      }
                             pl = pw->daughter(0);


        		     if(abs(pl->pdgId())<6) 
                             {
                                 ptGenV_3 = pw->pt();
                                 etaGenV_3 = pw->eta();
                                 phiGenV_3 = pw->phi();
                                 massGenV_3 = pw->mass();
                                 status_3=4;  
                             }
        		     if(abs(pl->pdgId())==24) 
                             {
                                 status_3=5;  
                             }
			   }//end of loop on w daughter
                       }//end of if w+	
			 
                   }//end of loop on Radion daughter
 		}//end of if Radion

        }//end of loop on gen
   }//end of MC Info
// *************************End of Gen Level Information******************//


    if(hadronicVs->size()!= 0){
       const reco::Candidate& metCand = metHandle->at(0);
       nLooseMu = loosemus->size();
       nLooseEle = looseels->size();
       nVetoMu = vetomus->size();
       nVetoEle = vetoels->size();

       edm::Handle<reco::VertexCollection> vertices;
       iEvent.getByToken(vtxToken_, vertices);
       if (nVetoMu != 0) return;
       if (hadronicVs->size()<2) return;
       if (vertices->empty()) return; // skip the event if no PV found
       nVtx = vertices->size();
       reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
       for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx) {
               if (  /// !vtx->isFake() &&
                     !(vtx->chi2()==0 && vtx->ndof()==0) 
	             &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
	             && fabs(vtx->position().Z())<=24.0) {
                     firstGoodVertex = vtx;
                     break;
                    }           
       }
       if ( firstGoodVertex==vertices->end() ) return; // skip event if there are no good PVs


           // ************************* MET ********************** //
              iEvent.getByToken(metInputToken_ , METs_ );
                addTypeICorr(iEvent);
                for (const pat::MET &met : *METs_) {
                        const float rawPt	 = met.uncorPt();
		    const float rawPhi   = met.uncorPhi();
		    const float rawSumEt = met.uncorSumEt();
                        TVector2 rawMET_;
                        rawMET_.SetMagPhi (rawPt, rawPhi );
                        Double_t rawPx = rawMET_.Px();
                        Double_t rawPy = rawMET_.Py();
                        Double_t rawEt = std::hypot(rawPx,rawPy);
            	    METraw_et = rawEt;
        	   	    METraw_phi = rawPhi;
        	        	    METraw_sumEt = rawSumEt;
                        double pxcorr = rawPx+TypeICorrMap_["corrEx"];
                        double pycorr = rawPy+TypeICorrMap_["corrEy"];
                        double et     = std::hypot(pxcorr,pycorr);
                        double sumEtcorr = rawSumEt+TypeICorrMap_["corrSumEt"];
                        TLorentzVector corrmet; corrmet.SetPxPyPzE(pxcorr,pycorr,0.,et);
            	    useless = sumEtcorr;
            	    useless = rawEt;
            	    MET_et = et;
            	    MET_phi = corrmet.Phi();
            	    MET_sumEt = sumEtcorr;
            	    MET_corrPx = TypeICorrMap_["corrEx"];
            	    MET_corrPy = TypeICorrMap_["corrEy"]; 
                }
           // ***************************************************************** //  
 
       /// For the time being, set these to 1
        triggerWeight=1.0;
        pileupWeight=1.0;
        double targetEvents = targetLumiInvPb_*crossSectionPb_;
        lumiWeight = targetEvents/originalNEvents_;
        met          = metCand.pt();
        metPhi       = metCand.phi();

                              if(leptonicVs->size()!= 0){  
        const reco::Candidate& leptonicV = leptonicVs->at(0);
 
        if( leptonicV.daughter(0)->isElectron()||leptonicV.daughter(1)->isElectron() ) {
                       const pat::Electron *el1 = leptonicV.daughter(0)->isElectron() ?
                                                  (pat::Electron*)leptonicV.daughter(0):
                                                  (pat::Electron*)leptonicV.daughter(1);
                double etaSC1         = el1->superCluster()->eta();
                double d01            = (-1)*el1->gsfTrack()->dxy(firstGoodVertex->position());
                isHEEP = false;
                et = el1->energy()!=0. ? el1->et()/el1->energy()*el1->caloEnergy() : 0.;
                if( et > 35. ) {
                     if( fabs(etaSC1) < 1.4442 ){
                        iso = el1->dr03EcalRecHitSumEt() + el1->dr03HcalDepth1TowerSumEt();
                        isoCut = 2 + 0.03*et + 0.28*fastJetRho;
                        if( el1->ecalDriven() == 1 && dEtaInSeed( el1 ) < 0.004 && el1->deltaPhiSuperClusterTrackAtVtx() < 0.06 &&
                         el1->hadronicOverEm() < (1./el1->superCluster()->energy()+0.05) &&
                         (el1->full5x5_e2x5Max()/el1->full5x5_e5x5() > 0.94 || el1->full5x5_e1x5()/el1->full5x5_e5x5() > 0.83) &&
                         el1->dr03TkSumPt() < 5. && el1->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 && iso < isoCut && fabs(d01) < 0.02 ) isHEEP = true;
                     }
                     if( fabs(etaSC1) > 1.566 && fabs(etaSC1) < 2.5 ){
                        iso = el1->dr03EcalRecHitSumEt() + el1->dr03HcalDepth1TowerSumEt();
                        if( et <= 50 )
                                isoCut = 2.5 + 0.28*fastJetRho;
                        else
                                isoCut = 2.5+0.03*(et-50.) + 0.28*fastJetRho;
                        if( el1->ecalDriven() == 1 && dEtaInSeed( el1 ) < 0.006 && el1->deltaPhiSuperClusterTrackAtVtx() < 0.06 &&
                         el1->hadronicOverEm() < (5./el1->superCluster()->energy()+0.05) && el1->full5x5_sigmaIetaIeta() < 0.03 &&
                         el1->dr03TkSumPt() < 5. && el1->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&
                         iso < isoCut && fabs(d01) < 0.05 ) isHEEP = true;
                     }
        }
        }
                                     }
                                else{
                              isHEEP = false;
                                    }
////////////////////////JEC for AK8/////////////////////////////////
        reco::Candidate::LorentzVector uncorrPrunedJet;
        bool doPuppi  = iEvent.getByToken(puppijetInputToken_, puppijets_ );
        if( doPuppi ){//1
         for(size_t ij=0; ij<puppijets_->size();ij++){
           corr_AK8puppi[ij] = 1;
           corr_AK8puppiSD[ij] = 1;
           const pat::Jet& hadronicVa = puppijets_->at(ij);
           reco::Candidate::LorentzVector uncorrJet;
           if(not isJEC_) doCorrOnTheFly_ = false;
           if( doCorrOnTheFly_ ){
              uncorrJet = hadronicVa.correctedP4(0);
              jecAK8puppi_->setJetEta( uncorrJet.eta()          );
              jecAK8puppi_->setJetPt ( uncorrJet.pt()           );
              jecAK8puppi_->setJetE  ( uncorrJet.energy()       );
              jecAK8puppi_->setRho   (fastJetRho);
              jecAK8puppi_->setNPV   (nVtx);
              jecAK8puppi_->setJetA  (hadronicVa.jetArea());
              corr_AK8puppi[ij] = jecAK8puppi_->getCorrection();
              jecAK8puppiGroomed_->setJetEta( uncorrJet.eta()          );
              jecAK8puppiGroomed_->setJetPt ( uncorrJet.pt()           );
              jecAK8puppiGroomed_->setJetE  ( uncorrJet.energy()       );
              jecAK8puppiGroomed_->setRho   (fastJetRho);
              jecAK8puppiGroomed_->setNPV   (nVtx);
              jecAK8puppiGroomed_->setJetA  (hadronicVa.jetArea());
              corr_AK8puppiSD[ij] = jecAK8puppiGroomed_->getCorrection();
           }
           else{uncorrJet = hadronicVa.p4();}

           if(ij<8){
              jetAK8puppi_pt1[ij] = corr_AK8puppi[ij]*uncorrJet.pt();
              jetAK8puppi_mass1[ij] = corr_AK8puppi[ij]*uncorrJet.mass();
              jetAK8puppi_eta1[ij] = uncorrJet.eta();
              jetAK8puppi_jec1[ij] = corr_AK8puppi[ij];
              jetAK8puppiSD_jec1[ij] = corr_AK8puppiSD[ij];
           }
         }

         int usenumber3 = -1; double pt_larger=0;
         int numvhad = puppijets_->size();
         for( int inum = 0; inum< numvhad; inum++){
           const pat::Jet& Vpuppi = puppijets_->at(inum);
           if(looseJetID(Vpuppi)<1) continue;     
           if(jetAK8puppi_pt1[inum] > pt_larger && fabs(jetAK8puppi_eta1[inum])<2.4 && inum<8) {pt_larger = jetAK8puppi_pt1[inum]; usenumber3 = inum; continue;}
        }
       if (usenumber3>-1) {//2
        const pat::Jet& hadronicVpuppi = puppijets_->at(usenumber3);
                IDLoose = looseJetID(hadronicVpuppi);
                IDTight = tightJetID(hadronicVpuppi);
                jetAK8puppi_ptJEC       = jetAK8puppi_pt1[usenumber3]; // unpruned corrected jet pt
                jetAK8puppi_eta     = jetAK8puppi_eta1[usenumber3]; // unpruned (w/o jec) jet eta
                jetAK8puppi_phi      = hadronicVpuppi.phi(); // unpruned (w/o jec) jet phi
                jetAK8puppi_tau1         = hadronicVpuppi.userFloat("NjettinessAK8:tau1");
                jetAK8puppi_tau2         = hadronicVpuppi.userFloat("NjettinessAK8:tau2");
                jetAK8puppi_tau3         = hadronicVpuppi.userFloat("NjettinessAK8:tau3");
                jetAK8puppi_tau21        = jetAK8puppi_tau2/jetAK8puppi_tau1;
                jetAK8puppi_tau4         = hadronicVpuppi.userFloat("NjettinessAK8:tau4");
                jetAK8puppi_tau42        = jetAK8puppi_tau4/jetAK8puppi_tau2;
                jetAK8puppi_sd       =  hadronicVpuppi.userFloat("ak8PFJetsCHSSoftDropMass"); // uncorrected pruned mass
                jetAK8puppi_sdJEC  =corr_AK8puppiSD[usenumber3]*jetAK8puppi_sd;
                Double_t gencorrect=1.0;
                Double_t recocorrect_0eta1p3=1.0;
                Double_t recocorrect_1p3eta2p5=1.0;
                gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC*0.08,-1.2);
                recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC+3.449e-07*pow(jetAK8puppi_ptJEC,2)-2.681e-10*pow(jetAK8puppi_ptJEC,3)+8.674e-14*pow(jetAK8puppi_ptJEC,4)-1.001e-17*pow(jetAK8puppi_ptJEC,5);
                recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC+8.37e-07*pow(jetAK8puppi_ptJEC,2)-5.204e-10*pow(jetAK8puppi_ptJEC,3)+1.454e-13*pow(jetAK8puppi_ptJEC,4)-1.504e-17*pow(jetAK8puppi_ptJEC,5);
                if (fabs(jetAK8puppi_eta)<=1.3){jetAK8puppi_sdcorr=jetAK8puppi_sd*gencorrect*recocorrect_0eta1p3;}
                else if (fabs(jetAK8puppi_eta)<2.5 && fabs(jetAK8puppi_eta)>1.3){jetAK8puppi_sdcorr=jetAK8puppi_sd*gencorrect*recocorrect_1p3eta2p5;}
           
           int usenumber2 = -1; double pt_larger2=0;
           for( int inum = 0; inum< numvhad; inum++){
               const pat::Jet& Vpuppi = puppijets_->at(inum);
               if(looseJetID(Vpuppi)<1) continue;
               if(jetAK8puppi_pt1[inum] > pt_larger2 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber3 && inum<8) {pt_larger2 = jetAK8puppi_pt1[inum]; usenumber2 = inum; continue;}
           }
           
           
          if(usenumber2>-1)  {
               const pat::Jet& hadronicVpuppi_2 = puppijets_->at(usenumber2);
                IDLoose_2 = looseJetID(hadronicVpuppi_2);
                IDTight_2 = tightJetID(hadronicVpuppi_2);
               jetAK8puppi_ptJEC_2       = jetAK8puppi_pt1[usenumber2]; // unpruned corrected jet pt
               jetAK8puppi_eta_2     = jetAK8puppi_eta1[usenumber2]; // unpruned (w/o jec) jet eta
               jetAK8puppi_phi_2      = hadronicVpuppi_2.phi(); // unpruned (w/o jec) jet phi
               jetAK8puppi_tau1_2         = hadronicVpuppi_2.userFloat("NjettinessAK8:tau1");
               jetAK8puppi_tau2_2         = hadronicVpuppi_2.userFloat("NjettinessAK8:tau2");
               jetAK8puppi_tau3_2         = hadronicVpuppi_2.userFloat("NjettinessAK8:tau3");
               jetAK8puppi_tau21_2        = jetAK8puppi_tau2_2/jetAK8puppi_tau1_2;
               jetAK8puppi_tau4_2         = hadronicVpuppi_2.userFloat("NjettinessAK8:tau4");
               jetAK8puppi_tau42_2        = jetAK8puppi_tau4_2/jetAK8puppi_tau2_2;
               
               jetAK8puppi_sd_2       =  hadronicVpuppi_2.userFloat("ak8PFJetsCHSSoftDropMass"); // uncorrected pruned mass
               jetAK8puppi_sdJEC_2  =corr_AK8puppiSD[usenumber2]*jetAK8puppi_sd_2;
               Double_t gencorrect=1.0;
               Double_t recocorrect_0eta1p3=1.0;
               Double_t recocorrect_1p3eta2p5=1.0;
               gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_2*0.08,-1.2);
               recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_2+3.449e-07*pow(jetAK8puppi_ptJEC_2,2)-2.681e-10*pow(jetAK8puppi_ptJEC_2,3)+8.674e-14*pow(jetAK8puppi_ptJEC_2,4)-1.001e-17*pow(jetAK8puppi_ptJEC_2,5);
               recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_2+8.37e-07*pow(jetAK8puppi_ptJEC_2,2)-5.204e-10*pow(jetAK8puppi_ptJEC_2,3)+1.454e-13*pow(jetAK8puppi_ptJEC_2,4)-1.504e-17*pow(jetAK8puppi_ptJEC_2,5);
               if (fabs(jetAK8puppi_eta_2)<=1.3){jetAK8puppi_sdcorr_2=jetAK8puppi_sd_2*gencorrect*recocorrect_0eta1p3;}
               else if (fabs(jetAK8puppi_eta_2)<2.5 && fabs(jetAK8puppi_eta_2)>1.3){jetAK8puppi_sdcorr_2=jetAK8puppi_sd_2*gencorrect*recocorrect_1p3eta2p5;}
           }
           


           int usenumber1 = -1; double pt_larger3=0;
           for( int inum = 0; inum< numvhad; inum++){
               const pat::Jet& Vpuppi = puppijets_->at(inum);
               if(looseJetID(Vpuppi)<1) continue;
               if(jetAK8puppi_pt1[inum] > pt_larger3 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber3 && inum != usenumber2 && inum<8) {pt_larger3 = jetAK8puppi_pt1[inum]; usenumber1 = inum; continue;}
           }
           
           
          if(usenumber1>-1)  {
               const pat::Jet& hadronicVpuppi_3 = puppijets_->at(usenumber1);
                IDLoose_3 = looseJetID(hadronicVpuppi_3);
                IDTight_3 = tightJetID(hadronicVpuppi_3);
               jetAK8puppi_ptJEC_3       = jetAK8puppi_pt1[usenumber1]; // unpruned corrected jet pt
               jetAK8puppi_eta_3     = jetAK8puppi_eta1[usenumber1]; // unpruned (w/o jec) jet eta
               jetAK8puppi_phi_3      = hadronicVpuppi_3.phi(); // unpruned (w/o jec) jet phi
               jetAK8puppi_tau1_3         = hadronicVpuppi_3.userFloat("NjettinessAK8:tau1");
               jetAK8puppi_tau2_3         = hadronicVpuppi_3.userFloat("NjettinessAK8:tau2");
               jetAK8puppi_tau3_3         = hadronicVpuppi_3.userFloat("NjettinessAK8:tau3");
               jetAK8puppi_tau21_3        = jetAK8puppi_tau2_3/jetAK8puppi_tau1_3;
               jetAK8puppi_tau4_3         = hadronicVpuppi_3.userFloat("NjettinessAK8:tau4");
               jetAK8puppi_tau42_3        = jetAK8puppi_tau4_3/jetAK8puppi_tau2_3;
               
               jetAK8puppi_sd_3       =  hadronicVpuppi_3.userFloat("ak8PFJetsCHSSoftDropMass"); // uncorrected pruned mass
               jetAK8puppi_sdJEC_3  =corr_AK8puppiSD[usenumber1]*jetAK8puppi_sd_3;
               Double_t gencorrect=1.0;
               Double_t recocorrect_0eta1p3=1.0;
               Double_t recocorrect_1p3eta2p5=1.0;
               gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_3*0.08,-1.2);
               recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_3+3.449e-07*pow(jetAK8puppi_ptJEC_3,2)-2.681e-10*pow(jetAK8puppi_ptJEC_3,3)+8.674e-14*pow(jetAK8puppi_ptJEC_3,4)-1.001e-17*pow(jetAK8puppi_ptJEC_3,5);
               recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_3+8.37e-07*pow(jetAK8puppi_ptJEC_3,2)-5.204e-10*pow(jetAK8puppi_ptJEC_3,3)+1.454e-13*pow(jetAK8puppi_ptJEC_3,4)-1.504e-17*pow(jetAK8puppi_ptJEC_3,5);
               if (fabs(jetAK8puppi_eta_3)<=1.3){jetAK8puppi_sdcorr_3=jetAK8puppi_sd_3*gencorrect*recocorrect_0eta1p3;}
               else if (fabs(jetAK8puppi_eta_3)<2.5 && fabs(jetAK8puppi_eta_3)>1.3){jetAK8puppi_sdcorr_3=jetAK8puppi_sd_3*gencorrect*recocorrect_1p3eta2p5;}
           }
           

           int usenumber4 = -1; double pt_larger4=0;
           for( int inum = 0; inum< numvhad; inum++){
               const pat::Jet& Vpuppi = puppijets_->at(inum);
               if(looseJetID(Vpuppi)<1) continue;
               if(jetAK8puppi_pt1[inum] > pt_larger4 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber3 && inum != usenumber2 && inum != usenumber1 && inum<8) {pt_larger4 = jetAK8puppi_pt1[inum]; usenumber4 = inum; continue;}
           }


          if(usenumber4>-1)  {
               const pat::Jet& hadronicVpuppi_4 = puppijets_->at(usenumber4);
                IDLoose_4 = looseJetID(hadronicVpuppi_4);
                IDTight_4 = tightJetID(hadronicVpuppi_4);
               jetAK8puppi_ptJEC_4       = jetAK8puppi_pt1[usenumber4]; // unpruned corrected jet pt
               jetAK8puppi_eta_4     = jetAK8puppi_eta1[usenumber4]; // unpruned (w/o jec) jet eta
               jetAK8puppi_phi_4      = hadronicVpuppi_4.phi(); // unpruned (w/o jec) jet phi
               jetAK8puppi_tau1_4         = hadronicVpuppi_4.userFloat("NjettinessAK8:tau1");
               jetAK8puppi_tau2_4         = hadronicVpuppi_4.userFloat("NjettinessAK8:tau2");
               jetAK8puppi_tau3_4         = hadronicVpuppi_4.userFloat("NjettinessAK8:tau3");
               jetAK8puppi_tau21_4        = jetAK8puppi_tau2_4/jetAK8puppi_tau1_4;
               jetAK8puppi_tau4_4         = hadronicVpuppi_4.userFloat("NjettinessAK8:tau4");
               jetAK8puppi_tau42_4        = jetAK8puppi_tau4_4/jetAK8puppi_tau2_4;

               jetAK8puppi_sd_4       =  hadronicVpuppi_4.userFloat("ak8PFJetsCHSSoftDropMass"); // uncorrected pruned mass
               jetAK8puppi_sdJEC_4  =corr_AK8puppiSD[usenumber4]*jetAK8puppi_sd_4;
               Double_t gencorrect=1.0;
               Double_t recocorrect_0eta1p3=1.0;
               Double_t recocorrect_1p3eta2p5=1.0;
               gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_4*0.08,-1.2);
               recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_4+3.449e-07*pow(jetAK8puppi_ptJEC_4,2)-2.681e-10*pow(jetAK8puppi_ptJEC_4,3)+8.674e-14*pow(jetAK8puppi_ptJEC_4,4)-1.001e-17*pow(jetAK8puppi_ptJEC_4,5);
               recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_4+8.37e-07*pow(jetAK8puppi_ptJEC_4,2)-5.204e-10*pow(jetAK8puppi_ptJEC_4,3)+1.454e-13*pow(jetAK8puppi_ptJEC_4,4)-1.504e-17*pow(jetAK8puppi_ptJEC_4,5);
               if (fabs(jetAK8puppi_eta_4)<=1.3){jetAK8puppi_sdcorr_4=jetAK8puppi_sd_4*gencorrect*recocorrect_0eta1p3;}
               else if (fabs(jetAK8puppi_eta_4)<2.5 && fabs(jetAK8puppi_eta_4)>1.3){jetAK8puppi_sdcorr_4=jetAK8puppi_sd_4*gencorrect*recocorrect_1p3eta2p5;}
           }


           int usenumber5 = -1; double pt_larger5=0;
           for( int inum = 0; inum< numvhad; inum++){
               const pat::Jet& Vpuppi = puppijets_->at(inum);
               if(looseJetID(Vpuppi)<1) continue;
               if(jetAK8puppi_pt1[inum] > pt_larger5 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber4 && inum != usenumber3 && inum != usenumber2 && inum != usenumber1 && inum<8) {pt_larger5 = jetAK8puppi_pt1[inum]; usenumber5 = inum; continue;}
           }


          if(usenumber5>-1)  {
               const pat::Jet& hadronicVpuppi_5 = puppijets_->at(usenumber5);
               jetAK8puppi_ptJEC_5       = jetAK8puppi_pt1[usenumber5]; // unpruned corrected jet pt
               jetAK8puppi_eta_5     = jetAK8puppi_eta1[usenumber5]; // unpruned (w/o jec) jet eta
               jetAK8puppi_phi_5      = hadronicVpuppi_5.phi(); // unpruned (w/o jec) jet phi
               jetAK8puppi_tau1_5         = hadronicVpuppi_5.userFloat("NjettinessAK8:tau1");
               jetAK8puppi_tau2_5         = hadronicVpuppi_5.userFloat("NjettinessAK8:tau2");
               jetAK8puppi_tau3_5         = hadronicVpuppi_5.userFloat("NjettinessAK8:tau3");
               jetAK8puppi_tau21_5        = jetAK8puppi_tau2_5/jetAK8puppi_tau1_5;
               jetAK8puppi_tau4_5         = hadronicVpuppi_5.userFloat("NjettinessAK8:tau4");
               jetAK8puppi_tau42_5        = jetAK8puppi_tau4_5/jetAK8puppi_tau2_5;

               jetAK8puppi_sd_5       =  hadronicVpuppi_5.userFloat("ak8PFJetsCHSSoftDropMass"); // uncorrected pruned mass
               jetAK8puppi_sdJEC_5  =corr_AK8puppiSD[usenumber5]*jetAK8puppi_sd_5;
               Double_t gencorrect=1.0;
               Double_t recocorrect_0eta1p3=1.0;
               Double_t recocorrect_1p3eta2p5=1.0;
               gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_5*0.08,-1.2);
               recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_5+3.449e-07*pow(jetAK8puppi_ptJEC_5,2)-2.681e-10*pow(jetAK8puppi_ptJEC_5,3)+8.674e-14*pow(jetAK8puppi_ptJEC_5,4)-1.001e-17*pow(jetAK8puppi_ptJEC_5,5);
               recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_5+8.37e-07*pow(jetAK8puppi_ptJEC_5,2)-5.204e-10*pow(jetAK8puppi_ptJEC_5,3)+1.454e-13*pow(jetAK8puppi_ptJEC_5,4)-1.504e-17*pow(jetAK8puppi_ptJEC_5,5);
               if (fabs(jetAK8puppi_eta_5)<=1.3){jetAK8puppi_sdcorr_5=jetAK8puppi_sd_5*gencorrect*recocorrect_0eta1p3;}
               else if (fabs(jetAK8puppi_eta_5)<2.5 && fabs(jetAK8puppi_eta_5)>1.3){jetAK8puppi_sdcorr_5=jetAK8puppi_sd_5*gencorrect*recocorrect_1p3eta2p5;}
           }



           int usenumber6 = -1; double pt_larger6=0;
           for( int inum = 0; inum< numvhad; inum++){
               const pat::Jet& Vpuppi = puppijets_->at(inum);
               if(looseJetID(Vpuppi)<1) continue;
               if(jetAK8puppi_pt1[inum] > pt_larger6 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber5 && inum != usenumber4 && inum != usenumber3 && inum != usenumber2 && inum != usenumber1 && inum<8) {pt_larger6 = jetAK8puppi_pt1[inum]; usenumber6 = inum; continue;}
           }


          if(usenumber6>-1)  {
               const pat::Jet& hadronicVpuppi_6 = puppijets_->at(usenumber6);
               jetAK8puppi_ptJEC_6       = jetAK8puppi_pt1[usenumber6]; // unpruned corrected jet pt
               jetAK8puppi_eta_6     = jetAK8puppi_eta1[usenumber6]; // unpruned (w/o jec) jet eta
               jetAK8puppi_phi_6      = hadronicVpuppi_6.phi(); // unpruned (w/o jec) jet phi
               jetAK8puppi_tau1_6         = hadronicVpuppi_6.userFloat("NjettinessAK8:tau1");
               jetAK8puppi_tau2_6         = hadronicVpuppi_6.userFloat("NjettinessAK8:tau2");
               jetAK8puppi_tau3_6         = hadronicVpuppi_6.userFloat("NjettinessAK8:tau3");
               jetAK8puppi_tau21_6        = jetAK8puppi_tau2_6/jetAK8puppi_tau1_6;
               jetAK8puppi_tau4_6         = hadronicVpuppi_6.userFloat("NjettinessAK8:tau4");
               jetAK8puppi_tau42_6        = jetAK8puppi_tau4_6/jetAK8puppi_tau2_6;

               jetAK8puppi_sd_6       =  hadronicVpuppi_6.userFloat("ak8PFJetsCHSSoftDropMass"); // uncorrected pruned mass
               jetAK8puppi_sdJEC_6  =corr_AK8puppiSD[usenumber6]*jetAK8puppi_sd_6;
               Double_t gencorrect=1.0;
               Double_t recocorrect_0eta1p3=1.0;
               Double_t recocorrect_1p3eta2p5=1.0;
               gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_6*0.08,-1.2);
               recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_6+3.449e-07*pow(jetAK8puppi_ptJEC_6,2)-2.681e-10*pow(jetAK8puppi_ptJEC_6,3)+8.674e-14*pow(jetAK8puppi_ptJEC_6,4)-1.001e-17*pow(jetAK8puppi_ptJEC_6,5);
               recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_6+8.37e-07*pow(jetAK8puppi_ptJEC_6,2)-5.204e-10*pow(jetAK8puppi_ptJEC_6,3)+1.454e-13*pow(jetAK8puppi_ptJEC_6,4)-1.504e-17*pow(jetAK8puppi_ptJEC_6,5);
               if (fabs(jetAK8puppi_eta_6)<=1.3){jetAK8puppi_sdcorr_6=jetAK8puppi_sd_6*gencorrect*recocorrect_0eta1p3;}
               else if (fabs(jetAK8puppi_eta_6)<2.5 && fabs(jetAK8puppi_eta_6)>1.3){jetAK8puppi_sdcorr_6=jetAK8puppi_sd_6*gencorrect*recocorrect_1p3eta2p5;}
           }


           int usenumber7 = -1; double pt_larger7=0;
           for( int inum = 0; inum< numvhad; inum++){
               const pat::Jet& Vpuppi = puppijets_->at(inum);
               if(looseJetID(Vpuppi)<1) continue;
               if(jetAK8puppi_pt1[inum] > pt_larger7 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber6 && inum != usenumber5 && inum != usenumber4 && inum != usenumber3 && inum != usenumber2 && inum != usenumber1 && inum<8) {pt_larger7 = jetAK8puppi_pt1[inum]; usenumber7 = inum; continue;}
           }


          if(usenumber7>-1)  {
               const pat::Jet& hadronicVpuppi_7 = puppijets_->at(usenumber7);
               jetAK8puppi_ptJEC_7       = jetAK8puppi_pt1[usenumber7]; // unpruned corrected jet pt
               jetAK8puppi_eta_7     = jetAK8puppi_eta1[usenumber7]; // unpruned (w/o jec) jet eta
               jetAK8puppi_phi_7      = hadronicVpuppi_7.phi(); // unpruned (w/o jec) jet phi
               jetAK8puppi_tau1_7         = hadronicVpuppi_7.userFloat("NjettinessAK8:tau1");
               jetAK8puppi_tau2_7         = hadronicVpuppi_7.userFloat("NjettinessAK8:tau2");
               jetAK8puppi_tau3_7         = hadronicVpuppi_7.userFloat("NjettinessAK8:tau3");
               jetAK8puppi_tau21_7        = jetAK8puppi_tau2_7/jetAK8puppi_tau1_7;
               jetAK8puppi_tau4_7         = hadronicVpuppi_7.userFloat("NjettinessAK8:tau4");
               jetAK8puppi_tau42_7        = jetAK8puppi_tau4_7/jetAK8puppi_tau2_7;

               jetAK8puppi_sd_7       =  hadronicVpuppi_7.userFloat("ak8PFJetsCHSSoftDropMass"); // uncorrected pruned mass
               jetAK8puppi_sdJEC_7  =corr_AK8puppiSD[usenumber7]*jetAK8puppi_sd_7;
               Double_t gencorrect=1.0;
               Double_t recocorrect_0eta1p3=1.0;
               Double_t recocorrect_1p3eta2p5=1.0;
               gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_7*0.08,-1.2);
               recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_7+3.449e-07*pow(jetAK8puppi_ptJEC_7,2)-2.681e-10*pow(jetAK8puppi_ptJEC_7,3)+8.674e-14*pow(jetAK8puppi_ptJEC_7,4)-1.001e-17*pow(jetAK8puppi_ptJEC_7,5);
               recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_7+8.37e-07*pow(jetAK8puppi_ptJEC_7,2)-5.204e-10*pow(jetAK8puppi_ptJEC_7,3)+1.454e-13*pow(jetAK8puppi_ptJEC_7,4)-1.504e-17*pow(jetAK8puppi_ptJEC_7,5);
               if (fabs(jetAK8puppi_eta_7)<=1.3){jetAK8puppi_sdcorr_7=jetAK8puppi_sd_7*gencorrect*recocorrect_0eta1p3;}
               else if (fabs(jetAK8puppi_eta_7)<2.5 && fabs(jetAK8puppi_eta_7)>1.3){jetAK8puppi_sdcorr_7=jetAK8puppi_sd_7*gencorrect*recocorrect_1p3eta2p5;}
           }


           int usenumber8 = -1; double pt_larger8=0;
           for( int inum = 0; inum< numvhad; inum++){
               const pat::Jet& Vpuppi = puppijets_->at(inum);
               if(looseJetID(Vpuppi)<1) continue;
               if(jetAK8puppi_pt1[inum] > pt_larger8 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber7 && inum != usenumber6 && inum != usenumber5 && inum != usenumber4 && inum != usenumber3 && inum != usenumber2 && inum != usenumber1 && inum<8) {pt_larger8 = jetAK8puppi_pt1[inum]; usenumber8 = inum; continue;}
           }


          if(usenumber8>-1)  {
               const pat::Jet& hadronicVpuppi_8 = puppijets_->at(usenumber8);
               jetAK8puppi_ptJEC_8       = jetAK8puppi_pt1[usenumber8]; // unpruned corrected jet pt
               jetAK8puppi_eta_8     = jetAK8puppi_eta1[usenumber8]; // unpruned (w/o jec) jet eta
               jetAK8puppi_phi_8      = hadronicVpuppi_8.phi(); // unpruned (w/o jec) jet phi
               jetAK8puppi_tau1_8         = hadronicVpuppi_8.userFloat("NjettinessAK8:tau1");
               jetAK8puppi_tau2_8         = hadronicVpuppi_8.userFloat("NjettinessAK8:tau2");
               jetAK8puppi_tau3_8         = hadronicVpuppi_8.userFloat("NjettinessAK8:tau3");
               jetAK8puppi_tau21_8        = jetAK8puppi_tau2_8/jetAK8puppi_tau1_8;
               jetAK8puppi_tau4_8         = hadronicVpuppi_8.userFloat("NjettinessAK8:tau4");
               jetAK8puppi_tau42_8        = jetAK8puppi_tau4_8/jetAK8puppi_tau2_8;

               jetAK8puppi_sd_8       =  hadronicVpuppi_8.userFloat("ak8PFJetsCHSSoftDropMass"); // uncorrected pruned mass
               jetAK8puppi_sdJEC_8  =corr_AK8puppiSD[usenumber8]*jetAK8puppi_sd_8;
               Double_t gencorrect=1.0;
               Double_t recocorrect_0eta1p3=1.0;
               Double_t recocorrect_1p3eta2p5=1.0;
               gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_8*0.08,-1.2);
               recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_8+3.449e-07*pow(jetAK8puppi_ptJEC_8,2)-2.681e-10*pow(jetAK8puppi_ptJEC_8,3)+8.674e-14*pow(jetAK8puppi_ptJEC_8,4)-1.001e-17*pow(jetAK8puppi_ptJEC_8,5);
               recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_8+8.37e-07*pow(jetAK8puppi_ptJEC_8,2)-5.204e-10*pow(jetAK8puppi_ptJEC_8,3)+1.454e-13*pow(jetAK8puppi_ptJEC_8,4)-1.504e-17*pow(jetAK8puppi_ptJEC_8,5);
               if (fabs(jetAK8puppi_eta_8)<=1.3){jetAK8puppi_sdcorr_8=jetAK8puppi_sd_8*gencorrect*recocorrect_0eta1p3;}
               else if (fabs(jetAK8puppi_eta_8)<2.5 && fabs(jetAK8puppi_eta_8)>1.3){jetAK8puppi_sdcorr_8=jetAK8puppi_sd_8*gencorrect*recocorrect_1p3eta2p5;}
           }



	int nak4 = 0;
 
        for (size_t ik=0; ik<ak4jets->size();ik++)
         {//3
            double corr = 1;
            reco::Candidate::LorentzVector uncorrJet;
             if( doCorrOnTheFly_ ){
	    uncorrJet = (*ak4jets)[ik].correctedP4(0);
              jecAK4_->setJetEta( uncorrJet.eta() );
              jecAK4_->setJetPt ( uncorrJet.pt() );
              jecAK4_->setJetE ( uncorrJet.energy() );
              jecAK4_->setRho ( fastJetRho );
              jecAK4_->setNPV ( vertices->size() );
              jecAK4_->setJetA ( (*ak4jets)[ik].jetArea() );
              corr = jecAK4_->getCorrection();
	    } else {uncorrJet = (*ak4jets)[ik].p4();}

            if( (corr*uncorrJet.pt())>20 && (fabs((*ak4jets)[ik].eta()) < 5.0) && looseJetID((*ak4jets)[ik])>0 && nak4<8){
                ak4jet_hf[nak4]=(*ak4jets)[ik].hadronFlavour();
                ak4jet_pf[nak4]=(*ak4jets)[ik].partonFlavour();
                ak4jet_pt[nak4] =  corr*uncorrJet.pt();
                ak4jet_pt_uncorr[nak4] =  uncorrJet.pt();  
                ak4jet_eta[nak4] = (*ak4jets)[ik].eta();
                ak4jet_phi[nak4] = (*ak4jets)[ik].phi();
                ak4jet_e[nak4] =   corr*uncorrJet.energy();
                ak4jet_csv[nak4] = (*ak4jets)[ik].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
                ak4jet_icsv[nak4] = (*ak4jets)[ik].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");   
                ak4jet_IDLoose[nak4] = looseJetID((*ak4jets)[ik]);
	      ak4jet_IDTight[nak4] = tightJetID((*ak4jets)[ik]);
                nak4 = nak4 + 1;
            }
         }//3
         
        

                TLorentzVector ghadronicVpuppi, gravitonpuppiJEC,ghadronicVpuppi_2, gravitonpuppiJEC_2, ghadronicVpuppi_3, gravitonpuppiJEC_3;
                ghadronicVpuppi.SetPtEtaPhiM(jetAK8puppi_ptJEC, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sdJEC);
                ghadronicVpuppi_2.SetPtEtaPhiM(jetAK8puppi_ptJEC_2, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sdJEC_2);
                ghadronicVpuppi_3.SetPtEtaPhiM(jetAK8puppi_ptJEC_3, jetAK8puppi_eta_3, jetAK8puppi_phi_3, jetAK8puppi_sdJEC_3);
	        gravitonpuppiJEC =  ghadronicVpuppi+ ghadronicVpuppi_2+ ghadronicVpuppi_3;
                candMasspuppiJEC     = gravitonpuppiJEC.Mag();
                delPhijetmet = deltaPhi(jetAK8puppi_phi, MET_phi);
                delPhijetmet_2 = deltaPhi(jetAK8puppi_phi_2, MET_phi);
                delPhijetmet_3 = deltaPhi(jetAK8puppi_phi_3, MET_phi);
           
           TLorentzVector lvw[3];
           lvw[0] = ghadronicVpuppi;
           lvw[1] = ghadronicVpuppi_2;
           lvw[2] = ghadronicVpuppi_3;

           Double_t Wpt[3];
           Wpt[0]=jetAK8puppi_ptJEC;
           Wpt[1]=jetAK8puppi_ptJEC_2;
           Wpt[2]=jetAK8puppi_ptJEC_3;

           Int_t *indexx=new Int_t[3];
           TMath::Sort(3,Wpt,indexx,1);
           massww[0] = (lvw[0]+lvw[1]).Mag();
           massww[1] = (lvw[0]+lvw[2]).Mag();
           massww[2] = (lvw[1]+lvw[2]).Mag();
           edm::Handle<edm::View<pat::Jet> > jetsAK8;
           
            iEvent.getByToken(jetsAK8Label_, jetsAK8);
            
            edm::View<pat::Jet>::const_iterator beginAK8 = jetsAK8->begin();
            edm::View<pat::Jet>::const_iterator endAK8 = jetsAK8->end();
            edm::View<pat::Jet>::const_iterator ijetAK8 = beginAK8;
       
           edm::View<pat::Jet>::const_iterator ijetAK8_j1,ijetAK8_j2,ijetAK8_j3;
           double drak8jetmatch1=10000.,drak8jetmatch2=10000.,drak8jetmatch3=10000.;
         
            for(ijetAK8 = beginAK8; ijetAK8 != endAK8; ++ijetAK8 ) {
                    if(ijetAK8->pt()>0){
                    double tmpdrak8jet1=deltaR(ijetAK8->eta(),ijetAK8->phi(),jetAK8puppi_eta,jetAK8puppi_phi);                  
                    double tmpdrak8jet2=deltaR(ijetAK8->eta(),ijetAK8->phi(),jetAK8puppi_eta_2,jetAK8puppi_phi_2);
                    double tmpdrak8jet3=deltaR(ijetAK8->eta(),ijetAK8->phi(),jetAK8puppi_eta_3,jetAK8puppi_phi_3);
                    if (tmpdrak8jet1<drak8jetmatch1) {drak8jetmatch1=tmpdrak8jet1; ijetAK8_j1=ijetAK8;}
                    if (tmpdrak8jet2<drak8jetmatch2) {drak8jetmatch2=tmpdrak8jet2; ijetAK8_j2=ijetAK8;}
                    if (tmpdrak8jet3<drak8jetmatch3) {drak8jetmatch3=tmpdrak8jet3; ijetAK8_j3=ijetAK8;}
                                       }  
            }

           auto const & sdSubjetsPuppi_1 = ijetAK8_j1->subjets("SoftDropPuppi");
           Int_t nsj1=0;
           for ( auto const & puppiSDSJ_1 : sdSubjetsPuppi_1 ) {
                      if (nsj1==0)
                           ak8sj11.SetPtEtaPhiM(puppiSDSJ_1->correctedP4(0).pt(),puppiSDSJ_1->correctedP4(0).eta(),puppiSDSJ_1->correctedP4(0).phi(),puppiSDSJ_1->correctedP4(0).mass());
                       if (nsj1==1)
                           ak8sj12.SetPtEtaPhiM(puppiSDSJ_1->correctedP4(0).pt(),puppiSDSJ_1->correctedP4(0).eta(),puppiSDSJ_1->correctedP4(0).phi(),puppiSDSJ_1->correctedP4(0).mass());
                       if (nsj1==2)
                           ak8sj13.SetPtEtaPhiM(puppiSDSJ_1->correctedP4(0).pt(),puppiSDSJ_1->correctedP4(0).eta(),puppiSDSJ_1->correctedP4(0).phi(),puppiSDSJ_1->correctedP4(0).mass());
                       if (nsj1==3)
                           ak8sj14.SetPtEtaPhiM(puppiSDSJ_1->correctedP4(0).pt(),puppiSDSJ_1->correctedP4(0).eta(),puppiSDSJ_1->correctedP4(0).phi(),puppiSDSJ_1->correctedP4(0).mass());
                       if (nsj1==4)
                           ak8sj15.SetPtEtaPhiM(puppiSDSJ_1->correctedP4(0).pt(),puppiSDSJ_1->correctedP4(0).eta(),puppiSDSJ_1->correctedP4(0).phi(),puppiSDSJ_1->correctedP4(0).mass());
                   nsj1++;
               }

           auto const & sdSubjetsPuppi_2 = ijetAK8_j2->subjets("SoftDropPuppi");
           Int_t nsj2=0;
           for ( auto const & puppiSDSJ_2 : sdSubjetsPuppi_2 ) {
                      if (nsj2==0)
                           ak8sj21.SetPtEtaPhiM(puppiSDSJ_2->correctedP4(0).pt(),puppiSDSJ_2->correctedP4(0).eta(),puppiSDSJ_2->correctedP4(0).phi(),puppiSDSJ_2->correctedP4(0).mass());
                       if (nsj2==1)
                           ak8sj22.SetPtEtaPhiM(puppiSDSJ_2->correctedP4(0).pt(),puppiSDSJ_2->correctedP4(0).eta(),puppiSDSJ_2->correctedP4(0).phi(),puppiSDSJ_2->correctedP4(0).mass());
                       if (nsj2==2)
                           ak8sj23.SetPtEtaPhiM(puppiSDSJ_2->correctedP4(0).pt(),puppiSDSJ_2->correctedP4(0).eta(),puppiSDSJ_2->correctedP4(0).phi(),puppiSDSJ_2->correctedP4(0).mass());
                       if (nsj2==3)
                           ak8sj24.SetPtEtaPhiM(puppiSDSJ_2->correctedP4(0).pt(),puppiSDSJ_2->correctedP4(0).eta(),puppiSDSJ_2->correctedP4(0).phi(),puppiSDSJ_2->correctedP4(0).mass());
                       if (nsj2==4)
                           ak8sj25.SetPtEtaPhiM(puppiSDSJ_2->correctedP4(0).pt(),puppiSDSJ_2->correctedP4(0).eta(),puppiSDSJ_2->correctedP4(0).phi(),puppiSDSJ_2->correctedP4(0).mass());
                   nsj2++;
               }
 
           auto const & sdSubjetsPuppi_3 = ijetAK8_j3->subjets("SoftDropPuppi");
           Int_t nsj3=0;
           for ( auto const & puppiSDSJ_3 : sdSubjetsPuppi_3 ) {
                      if (nsj3==0)
                           ak8sj31.SetPtEtaPhiM(puppiSDSJ_3->correctedP4(0).pt(),puppiSDSJ_3->correctedP4(0).eta(),puppiSDSJ_3->correctedP4(0).phi(),puppiSDSJ_3->correctedP4(0).mass());
                       if (nsj3==1)
                           ak8sj32.SetPtEtaPhiM(puppiSDSJ_3->correctedP4(0).pt(),puppiSDSJ_3->correctedP4(0).eta(),puppiSDSJ_3->correctedP4(0).phi(),puppiSDSJ_3->correctedP4(0).mass());
                       if (nsj3==2)
                           ak8sj33.SetPtEtaPhiM(puppiSDSJ_3->correctedP4(0).pt(),puppiSDSJ_3->correctedP4(0).eta(),puppiSDSJ_3->correctedP4(0).phi(),puppiSDSJ_3->correctedP4(0).mass());
                       if (nsj3==3)
                           ak8sj34.SetPtEtaPhiM(puppiSDSJ_3->correctedP4(0).pt(),puppiSDSJ_3->correctedP4(0).eta(),puppiSDSJ_3->correctedP4(0).phi(),puppiSDSJ_3->correctedP4(0).mass());
                       if (nsj3==4)
                           ak8sj35.SetPtEtaPhiM(puppiSDSJ_3->correctedP4(0).pt(),puppiSDSJ_3->correctedP4(0).eta(),puppiSDSJ_3->correctedP4(0).phi(),puppiSDSJ_3->correctedP4(0).mass());
                   nsj3++;
               }
         }//2
        }//1
Int_t Nj8=(jetAK8puppi_ptJEC>0)+(jetAK8puppi_ptJEC_2>0)+(jetAK8puppi_ptJEC_3>0)+(jetAK8puppi_ptJEC_4>0)+(jetAK8puppi_ptJEC_5>0)+(jetAK8puppi_ptJEC_6>0)+(jetAK8puppi_ptJEC_7>0)+(jetAK8puppi_ptJEC_8>0);

TLorentzVector AK81,AK82,AK83,AK84;
AK81.SetPtEtaPhiM(0,-99,-99,-99);
AK82.SetPtEtaPhiM(0,-99,-99,-99);
AK83.SetPtEtaPhiM(0,-99,-99,-99);
AK84.SetPtEtaPhiM(0,-99,-99,-99);
Double_t MJJ=-99;
Double_t MJJJ=-99;
Double_t MJJJJ=-99;
if(Nj8==2){
AK81.SetPtEtaPhiM(jetAK8puppi_ptJEC,jetAK8puppi_eta,jetAK8puppi_phi,jetAK8puppi_sdcorr);
AK82.SetPtEtaPhiM(jetAK8puppi_ptJEC_2,jetAK8puppi_eta_2,jetAK8puppi_phi_2,jetAK8puppi_sdcorr_2);
MJJ=(AK81+AK82).M();
          }
if(Nj8==3){
MJJJ=candMasspuppiJEC;
          }
if(Nj8==4){
AK81.SetPtEtaPhiM(jetAK8puppi_ptJEC,jetAK8puppi_eta,jetAK8puppi_phi,jetAK8puppi_sdcorr);
AK82.SetPtEtaPhiM(jetAK8puppi_ptJEC_2,jetAK8puppi_eta_2,jetAK8puppi_phi_2,jetAK8puppi_sdcorr_2);
AK83.SetPtEtaPhiM(jetAK8puppi_ptJEC_3,jetAK8puppi_eta_3,jetAK8puppi_phi_3,jetAK8puppi_sdcorr_3);
AK84.SetPtEtaPhiM(jetAK8puppi_ptJEC_4,jetAK8puppi_eta_4,jetAK8puppi_phi_4,jetAK8puppi_sdcorr_4);
MJJJJ=(AK81+AK82+AK83+AK84).M();
          }


         nVetoEle=isHEEP*nVetoEle;//Consider HEEP7 condition

         if(nVetoEle==0&&nVetoMu==0&&((Nj8==2&&MJJ>-1000&&IDLoose>0&&IDLoose_2>0)||(Nj8==3&&MJJJ>-1000&&IDLoose>0&&IDLoose_2>0&&IDLoose_3>0&&fabs(jetAK8puppi_eta_3)<2.4)||(Nj8==4&&MJJJJ>-1000&&IDLoose>0&&IDLoose_2>0&&IDLoose_3>0&&IDLoose_4>0&&fabs(jetAK8puppi_eta_4)<2.4)) && jetAK8puppi_ptJEC>400 && jetAK8puppi_ptJEC_2>200 && jetAK8puppi_sdcorr>40 && jetAK8puppi_sdcorr_2>40&& (HLT_Mu1>0 || HLT_Mu2>0 || HLT_Mu3>0 || HLT_Mu4>0 || HLT_Mu5>0 || HLT_Mu6>0 || HLT_Mu7>0 || HLT_Mu8>0 || HLT_Mu9>0) && fabs(jetAK8puppi_eta)<2.4&& fabs(jetAK8puppi_eta_2)<2.4 && passFilter_HBHE_>0 && passFilter_GlobalHalo_>0 && passFilter_HBHEIso_>0 && passFilter_ECALDeadCell_>0 && passFilter_GoodVtx_>0 && passFilter_badMuon_>0){

              outTree_->Fill();
                     }
       outTreew_->Fill();
   }
   else {
       outTreew_->Fill();
        }
}
//-------------------------------------------------------------------------------------------------------------------------------------//


void EDBRTreeMaker::setDummyValues() {
     npT=-1.;
     npIT=-1.;
     nBX=-1;
    ak8sj11.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj21.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj31.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj12.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj22.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj32.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj13.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj23.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj33.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj14.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj24.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj34.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj15.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj25.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj35.SetPtEtaPhiM(0,-99,-99,-99);

    for(int p=0;p<211;p++){
                pweight[p]=-99; 
                          }
     nLooseEle      =-99;
     nLooseMu       =-99;
     nVetoEle      =-99;
     nVetoMu       =-99;
     nVtx           = -99;
     triggerWeight  = -99;
     pileupWeight   = -99;
     lumiWeight     = -99;

     jetAK8puppi_ptJEC         = -99;
     jetAK8puppi_eta         = -99;
     jetAK8puppi_phi         = -99;
     jetAK8puppi_tau1         = -99;
     jetAK8puppi_tau2         = -99;
     jetAK8puppi_tau3         = -99;
     jetAK8puppi_tau21         = -99;
     jetAK8puppi_tau4         = -99;
     jetAK8puppi_tau42         = -99;
     jetAK8puppi_sd         = -99;
     jetAK8puppi_sdJEC         = -99; 
     jetAK8puppi_sdcorr         = -99;

     jetAK8puppi_ptJEC_2         = -99;
     jetAK8puppi_eta_2         = -99;
     jetAK8puppi_phi_2         = -99;
     jetAK8puppi_tau1_2         = -99;
     jetAK8puppi_tau2_2         = -99;
     jetAK8puppi_tau3_2         = -99;
     jetAK8puppi_tau21_2         = -99;
     jetAK8puppi_tau4_2         = -99;
     jetAK8puppi_tau42_2         = -99;
     jetAK8puppi_sd_2         = -99;
     jetAK8puppi_sdJEC_2         = -99;
     jetAK8puppi_sdcorr_2         = -99;
 
     jetAK8puppi_ptJEC_3         = -99;
     jetAK8puppi_eta_3         = -99;
     jetAK8puppi_phi_3         = -99;
     jetAK8puppi_tau1_3         = -99;
     jetAK8puppi_tau2_3         = -99;
     jetAK8puppi_tau3_3         = -99;
     jetAK8puppi_tau21_3         = -99;
     jetAK8puppi_tau4_3         = -99;
     jetAK8puppi_tau42_3         = -99;
     jetAK8puppi_sd_3         = -99;
     jetAK8puppi_sdJEC_3         = -99;
     jetAK8puppi_sdcorr_3         = -99;

     jetAK8puppi_ptJEC_4         = -99;
     jetAK8puppi_eta_4         = -99;
     jetAK8puppi_phi_4         = -99;
     jetAK8puppi_tau1_4         = -99;
     jetAK8puppi_tau2_4         = -99;
     jetAK8puppi_tau3_4         = -99;
     jetAK8puppi_tau21_4         = -99;
     jetAK8puppi_tau4_4         = -99;
     jetAK8puppi_tau42_4         = -99;
     jetAK8puppi_sd_4         = -99;
     jetAK8puppi_sdJEC_4         = -99;
     jetAK8puppi_sdcorr_4         = -99;

     jetAK8puppi_ptJEC_5         = -99;
     jetAK8puppi_eta_5         = -99;
     jetAK8puppi_phi_5         = -99;
     jetAK8puppi_tau1_5         = -99;
     jetAK8puppi_tau2_5         = -99;
     jetAK8puppi_tau3_5         = -99;
     jetAK8puppi_tau21_5         = -99;
     jetAK8puppi_tau4_5         = -99;
     jetAK8puppi_tau42_5         = -99;
     jetAK8puppi_sd_5         = -99;
     jetAK8puppi_sdJEC_5         = -99;
     jetAK8puppi_sdcorr_5         = -99;

     jetAK8puppi_ptJEC_6         = -99;
     jetAK8puppi_eta_6         = -99;
     jetAK8puppi_phi_6         = -99;
     jetAK8puppi_tau1_6         = -99;
     jetAK8puppi_tau2_6         = -99;
     jetAK8puppi_tau3_6         = -99;
     jetAK8puppi_tau21_6         = -99;
     jetAK8puppi_tau4_6         = -99;
     jetAK8puppi_tau42_6         = -99;
     jetAK8puppi_sd_6         = -99;
     jetAK8puppi_sdJEC_6         = -99;
     jetAK8puppi_sdcorr_6         = -99;

     jetAK8puppi_ptJEC_7         = -99;
     jetAK8puppi_eta_7         = -99;
     jetAK8puppi_phi_7         = -99;
     jetAK8puppi_tau1_7         = -99;
     jetAK8puppi_tau2_7         = -99;
     jetAK8puppi_tau3_7         = -99;
     jetAK8puppi_tau21_7         = -99;
     jetAK8puppi_tau4_7         = -99;
     jetAK8puppi_tau42_7         = -99;
     jetAK8puppi_sd_7         = -99;
     jetAK8puppi_sdJEC_7         = -99;
     jetAK8puppi_sdcorr_7         = -99;

     jetAK8puppi_ptJEC_8         = -99;
     jetAK8puppi_eta_8         = -99;
     jetAK8puppi_phi_8         = -99;
     jetAK8puppi_tau1_8         = -99;
     jetAK8puppi_tau2_8         = -99;
     jetAK8puppi_tau3_8         = -99;
     jetAK8puppi_tau21_8         = -99;
     jetAK8puppi_tau4_8         = -99;
     jetAK8puppi_tau42_8         = -99;
     jetAK8puppi_sd_8         = -99;
     jetAK8puppi_sdJEC_8         = -99;
     jetAK8puppi_sdcorr_8         = -99;

     met            = -99;
     metPhi         = -99;
     delPhijetmet =  -99;
     delPhijetmet_2 =  -99;
     delPhijetmet_3 =  -99;

     gen_gra_m      = -99;
     gen_gra_pt     = -99;
     gen_gra_eta     = -99;
     gen_rad_m      = -99;
     gen_rad_pt     = -99;
     gen_rad_eta     = -99;
     gen_ele_pt     = -99;
     gen_ele_eta    = -99;
     gen_ele_phi    = -99;
     gen_ele_e      = -99;
     gen_mu_pt     = -99;
     gen_mu_eta    = -99;
     gen_mu_phi    = -99;
     gen_mu_e      = -99;
     gen_ele_pt_2     = -99;
     gen_ele_eta_2    = -99;
     gen_ele_phi_2    = -99;
     gen_ele_e_2      = -99;
     gen_mu_pt_2     = -99;
     gen_mu_eta_2    = -99;
     gen_mu_phi_2    = -99;
     gen_mu_e_2      = -99;
     gen_ele_pt_3     = -99;
     gen_ele_eta_3    = -99;
     gen_ele_phi_3    = -99;
     gen_ele_e_3      = -99;
     gen_mu_pt_3     = -99;
     gen_mu_eta_3    = -99;
     gen_mu_phi_3    = -99;
     gen_mu_e_3      = -99;

     gen_tau_pt     = -99;
     gen_tau_eta    = -99;
     gen_tau_phi    = -99;
     gen_tau_e      = -99;
     gen_tau_pt_2     = -99;
     gen_tau_eta_2    = -99;
     gen_tau_phi_2    = -99;
     gen_tau_e_2      = -99;
     gen_tau_pt_3     = -99;
     gen_tau_eta_3    = -99;
     gen_tau_phi_3    = -99;
     gen_tau_e_3      = -99;

     gentop_pt  = -99;
     gentop_eta  = -99;
     gentop_phi  = -99;
     gentop_mass  = -99;
     genantitop_pt  = -99;
     genantitop_eta  = -99;
     genantitop_phi  = -99;
     genantitop_mass  = -99;
     ptGenVlep      = -99;
     etaGenVlep      = -99;
     phiGenVlep      = -99;
     massGenVlep      = -99;
     ptGenVlep_2      = -99;
     etaGenVlep_2      = -99;
     phiGenVlep_2      = -99;
     massGenVlep_2      = -99;
     ptGenVlep_3      = -99;
     etaGenVlep_3      = -99;
     phiGenVlep_3      = -99;
     massGenVlep_3      = -99;     
     ptGenV_2      = -99;
     etaGenV_2      = -99;
     phiGenV_2      = -99;
     massGenV_2      = -99;
     ptGenV_3      = -99;
     etaGenV_3      = -99;
     phiGenV_3      = -99;
     massGenV_3      = -99;
     ptGenVhad      = -99;
     etaGenVhad      = -99;
     phiGenVhad      = -99;
     massGenVhad      = -99;
     ptq11 =-99;
     etaq11 =-99;
     phiq11 =-99;
     massq11 =-99;
     ptq12 =-99;
     etaq12 =-99;
     phiq12 =-99;
     massq12 =-99; 
     ptq21 =-99;
     etaq21 =-99;
     phiq21 =-99;
     massq21 =-99;
     ptq22 =-99;
     etaq22 =-99;
     phiq22 =-99;
     massq22 =-99;
     ptq31 =-99;
     etaq31 =-99;
     phiq31 =-99;
     massq31 =-99;
     ptq32 =-99;
     etaq32 =-99;
     phiq32 =-99;
     massq32 =-99;
     status_1       =  -1;
     status_2       =  -1;
     status_3       =  -1;

     for(Int_t ii=0;ii<8;ii++){
        ak4jet_hf[ii] = -99;
        ak4jet_pf[ii] = -99;
        ak4jet_pt[ii] = -99;
        ak4jet_pt_uncorr[ii] = -99;
        ak4jet_eta[ii] = -99;
        ak4jet_phi[ii] = -99;
        ak4jet_e[ii] = -99;
        ak4jet_dr[ii] = -99;
        ak4jet_csv[ii] = -99;
        ak4jet_icsv[ii] = -99;
        ak4jet_IDLoose[ii] = -99;
        ak4jet_IDTight[ii] = -99;
    }
 

     jetAK8_mass = -99;
     jetAK8_pt = -99;
     jetAK8_jec = -99;
     jetAK8_mass1[0] = -99;
     jetAK8_mass1[1] = -99;
     jetAK8_mass1[2] = -99;
     jetAK8_SF_mass1[0] = -99;
     jetAK8_SF_mass1[1] = -99;
     jetAK8_SF_mass1[2] = -99;
     jetAK8_SF_mass2[0] = -99;
     jetAK8_SF_mass2[1] = -99;
     jetAK8_SF_mass2[2] = -99;
     jetAK8_pt1[0] = -99;
     jetAK8_pt1[1] = -99;
     jetAK8_pt1[2] = -99;
     jetAK8_jec1[0] = -99;
     jetAK8_jec1[1] = -99;
     jetAK8_jec1[2] = -99;
     corr_AK81[0] = -99;
     corr_AK81[1] = -99;
     corr_AK81[2] = -99;
     jetAK8_eta = -99;
     jetAK8_eta1[0] = -99;
     jetAK8_eta1[1] = -99;
     jetAK8_eta1[2] = -99;
     jetAK8_phi = -99;

     METraw_et = -99;
     METraw_phi = -99;
     METraw_sumEt = -99;
     MET_et = -99;
     MET_phi = -99;
     MET_sumEt = -99;
     MET_corrPx = -99;
     MET_corrPy = -99;

     candMasspuppiJEC     =  -99;
     massww[0] = -99;
     massww[1] = -99;
     massww[2] = -99;

     HLT_Mu1=-99;
     HLT_Mu2=-99;
     HLT_Mu3=-99;
     HLT_Mu4=-99;
     HLT_Mu5=-99;
     HLT_Mu6=-99;
     HLT_Mu7=-99;
     HLT_Mu8=-99;
     HLT_Mu9=-99;
     HLT_Mu10=-99;
     HLT_Mu11=-99;
     HLT_Mu12=-99;
     HLT_Mu13=-99;
     HLT_Mu14=-99;
     HLT_Mu15=-99;
     HLT_Mu16=-99;

     theWeight = -99;
     //nump = 0;
     //numm = 0;
     passFilter_HBHE_                  = false;
     passFilter_HBHEIso_               = false;
     passFilter_GlobalHalo_            = false;
     passFilter_ECALDeadCell_          = false;
     passFilter_GoodVtx_               = false;
     passFilter_EEBadSc_               = false;
     passFilter_badMuon_               = false;
     passFilter_badChargedHadron_      = false;

}

// ------------ method called once each job just before starting event loop  ------------
void 
EDBRTreeMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void EDBRTreeMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  muPaths1.clear();
  muPaths2.clear();
  muPaths3.clear();
  muPaths4.clear();
  muPaths5.clear();
  muPaths6.clear();
  muPaths7.clear();
  muPaths8.clear();
  muPaths9.clear();
  el1.clear();
  el2.clear();
  el3.clear();
  mu1.clear();
  mu2.clear();
  mu3.clear();
  mu4.clear();

  std::cout<<"-----begin-----"<<std::endl;
   bool changed;
   if ( !hltConfig.init(iRun, iSetup, "HLT", changed) ) {
        edm::LogError("HltAnalysis") << "Initialization of HLTConfigProvider failed!!";
       return;
      }
   for (size_t i = 0; i < muPaths1_.size(); i++) {
         std::vector<std::string> foundPaths1 = hltConfig.matched( hltConfig.triggerNames(), muPaths1_[i] );
         while ( !foundPaths1.empty() ){
               muPaths1.push_back( foundPaths1.back() );
               foundPaths1.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-1 Information **************\n";
   for (size_t i=0; i < muPaths1.size(); i++) std::cout << "\n HLT paths-1:   " << i<<"  "<<muPaths1[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths2_.size(); i++) {
         std::vector<std::string> foundPaths2 = hltConfig.matched( hltConfig.triggerNames(), muPaths2_[i] );
         while ( !foundPaths2.empty() ){
               muPaths2.push_back( foundPaths2.back() );
               foundPaths2.pop_back();
                                      }
                                                }

   std::cout<<"\n************** HLT-2 Information **************\n";
   for (size_t i=0; i < muPaths2.size(); i++) std::cout << "\n Muon paths-2:   " << i<<"  "<<muPaths2[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths3_.size(); i++) {
         std::vector<std::string> foundPaths3 = hltConfig.matched( hltConfig.triggerNames(), muPaths3_[i] );
         while ( !foundPaths3.empty() ){
               muPaths3.push_back( foundPaths3.back() );
               foundPaths3.pop_back();
                                      }
                                                }

   std::cout<<"\n************** HLT-3 Information **************\n";
   for (size_t i=0; i < muPaths3.size(); i++) std::cout << "\n Muon paths-3:   " << i<<"  "<<muPaths3[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths4_.size(); i++) {
         std::vector<std::string> foundPaths4 = hltConfig.matched( hltConfig.triggerNames(), muPaths4_[i] );
         while ( !foundPaths4.empty() ){
               muPaths4.push_back( foundPaths4.back() );
               foundPaths4.pop_back();
                                      }
                                                }

   std::cout<<"\n************** HLT-4 Information **************\n";
   for (size_t i=0; i < muPaths4.size(); i++) std::cout << "\n Muon paths-4:   " << i<<"  "<<muPaths4[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths5_.size(); i++) {
         std::vector<std::string> foundPaths5 = hltConfig.matched( hltConfig.triggerNames(), muPaths5_[i] );
         while ( !foundPaths5.empty() ){
               muPaths5.push_back( foundPaths5.back() );
               foundPaths5.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-5 Information **************\n";
   for (size_t i=0; i < muPaths5.size(); i++) std::cout << "\n Muon paths-5:   " << i<<"  "<<muPaths5[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths6_.size(); i++) {
         std::vector<std::string> foundPaths6 = hltConfig.matched( hltConfig.triggerNames(), muPaths6_[i] );
         while ( !foundPaths6.empty() ){
               muPaths6.push_back( foundPaths6.back() );
               foundPaths6.pop_back();
                                      }
                                                }

   std::cout<<"\n************** HLT-6 Information **************\n";
   for (size_t i=0; i < muPaths6.size(); i++) std::cout << "\n Muon paths-6:   " << i<<"  "<<muPaths6[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths7_.size(); i++) {
         std::vector<std::string> foundPaths7 = hltConfig.matched( hltConfig.triggerNames(), muPaths7_[i] );
         while ( !foundPaths7.empty() ){
               muPaths7.push_back( foundPaths7.back() );
               foundPaths7.pop_back();
                                      }
                                                }

   std::cout<<"\n************** HLT-7 Information **************\n";
   for (size_t i=0; i < muPaths7.size(); i++) std::cout << "\n Muon paths-7:   " << i<<"  "<<muPaths7[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths8_.size(); i++) {
         std::vector<std::string> foundPaths8 = hltConfig.matched( hltConfig.triggerNames(), muPaths8_[i] );
         while ( !foundPaths8.empty() ){
               muPaths8.push_back( foundPaths8.back() );
               foundPaths8.pop_back();
                                      }
                                                }

   std::cout<<"\n************** HLT-8 Information **************\n";
   for (size_t i=0; i < muPaths8.size(); i++) std::cout << "\n Muon paths-8:   " << i<<"  "<<muPaths8[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths9_.size(); i++) {
         std::vector<std::string> foundPaths9 = hltConfig.matched( hltConfig.triggerNames(), muPaths9_[i] );
         while ( !foundPaths9.empty() ){
               muPaths9.push_back( foundPaths9.back() );
               foundPaths9.pop_back();
                                      }
                                                }

   std::cout<<"\n************** HLT-9 Information **************\n";
   for (size_t i=0; i < muPaths9.size(); i++) std::cout << "\n Muon paths-9:   " << i<<"  "<<muPaths9[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < el1_.size(); i++) {
         std::vector<std::string> foundPaths10 = hltConfig.matched( hltConfig.triggerNames(), el1_[i] );
         while ( !foundPaths10.empty() ){
               el1.push_back( foundPaths10.back() );
               foundPaths10.pop_back();
                                      }         
                                                }
   std::cout<<"\n************** HLT-10 Information **************\n";
   for (size_t i=0; i < el1.size(); i++) std::cout << "\n HLT paths-10:   " << i<<"  "<<el1[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < el2_.size(); i++) {
         std::vector<std::string> foundPaths11 = hltConfig.matched( hltConfig.triggerNames(), el2_[i] );
         while ( !foundPaths11.empty() ){
               el2.push_back( foundPaths11.back() );
               foundPaths11.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-11 Information **************\n";
   for (size_t i=0; i < el2.size(); i++) std::cout << "\n HLT paths-11:   " << i<<"  "<<el2[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < el3_.size(); i++) {
         std::vector<std::string> foundPaths12 = hltConfig.matched( hltConfig.triggerNames(), el3_[i] );
         while ( !foundPaths12.empty() ){
               el3.push_back( foundPaths12.back() );
               foundPaths12.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-12 Information **************\n";
   for (size_t i=0; i < el3.size(); i++) std::cout << "\n HLT paths-12:   " << i<<"  "<<el3[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < mu1_.size(); i++) {
         std::vector<std::string> foundPaths13 = hltConfig.matched( hltConfig.triggerNames(), mu1_[i] );
         while ( !foundPaths13.empty() ){
               mu1.push_back( foundPaths13.back() );
               foundPaths13.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-13 Information **************\n";
   for (size_t i=0; i < mu1.size(); i++) std::cout << "\n HLT paths-13:   " << i<<"  "<<mu1[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";


   for (size_t i = 0; i < mu2_.size(); i++) {
         std::vector<std::string> foundPaths14 = hltConfig.matched( hltConfig.triggerNames(), mu2_[i] );
         while ( !foundPaths14.empty() ){
               mu2.push_back( foundPaths14.back() );
               foundPaths14.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-14 Information **************\n";
   for (size_t i=0; i < mu2.size(); i++) std::cout << "\n HLT paths-14:   " << i<<"  "<<mu2[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < mu3_.size(); i++) {
         std::vector<std::string> foundPaths15 = hltConfig.matched( hltConfig.triggerNames(), mu3_[i] );
         while ( !foundPaths15.empty() ){
               mu3.push_back( foundPaths15.back() );
               foundPaths15.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-15 Information **************\n";
   for (size_t i=0; i < mu3.size(); i++) std::cout << "\n HLT paths-15:   " << i<<"  "<<mu3[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < mu4_.size(); i++) {
         std::vector<std::string> foundPaths16 = hltConfig.matched( hltConfig.triggerNames(), mu4_[i] );
         while ( !foundPaths16.empty() ){
               mu4.push_back( foundPaths16.back() );
               foundPaths16.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-16 Information **************\n";
   for (size_t i=0; i < mu4.size(); i++) std::cout << "\n HLT paths-16:   " << i<<"  "<<mu4[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";



}

void EDBRTreeMaker::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
std::cout << "EDBRTreeMaker endJob()... endRun" << std::endl;
}


void
EDBRTreeMaker::endJob() {
  std::cout << "EDBRTreeMaker endJob()..." << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(EDBRTreeMaker);
