//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan  6 02:39:46 2016 by ROOT version 5.34/32
// from TTree EDBRCandidates/EDBR Candidates
// found on file: BulkGravWW750.root
//////////////////////////////////////////////////////////

#ifndef EDBR2PKUTree_h
#define EDBR2PKUTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLorentzVector.h"
//#include "DataFormats/Math/interface/deltaR.h"

#include <iostream>
#include <fstream>
using namespace std;
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxpassFilter_HBHE = 1;
   const Int_t kMaxpassFilter_HBHEIso = 1;
   const Int_t kMaxpassFilter_GlobalHalo = 1;
   const Int_t kMaxpassFilter_ECALDeadCell = 1;
   const Int_t kMaxpassFilter_GoodVtx = 1;
   const Int_t kMaxpassFilter_EEBadSc = 1;
   const Int_t kMaxpassFilter_badMuon = 1;
   const Int_t kMaxpassFilter_badChargedHadron = 1;

class EDBR2PKUTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   Int_t           ls;
   Int_t           run;
   Int_t           nLooseEle;
   Int_t           nLooseMu;
   Int_t           nbtag;
   Int_t           num_bJet;
   Int_t           num_bJet_loose;
   Int_t           num_bJet_tight;

   Float_t pfMET;
   Float_t pfMETPhi;
   Float_t l_pt;
   Float_t mtVlepnew;
   Float_t MTVlep;
   Float_t l_eta;
   Float_t l_phi;
   Float_t jet_pt;
   Float_t PTj;
   Float_t jet_eta;
   Float_t jet_phi;
   Float_t jet_mass_pruned;
   Float_t Mj;
   Float_t Mj_un;
   Float_t Mj_corr;
   Float_t jet_tau2tau1;
   Float_t tau21;
   Float_t jet_tau1_puppi;
   Float_t jet_tau2_puppi;
   Float_t jet_tau3_puppi;
   Float_t jet_tau4_puppi;
   Float_t tau31;
   Float_t tau41;
   Float_t tau32;
   Float_t tau43;
   Float_t W_pt;
   Float_t W_eta;
   Float_t W_phi;
   Float_t MJJJ;
   Float_t fjet2_pt;
   Float_t fjet2_btag;
   Float_t fjet3_pt;
   Float_t fjet3_btag; 
   Double_t MJ_j_18;
   Double_t MJJ_j_18;
   Double_t MJ_j_10;
   Double_t MJJ_j_10;
   Double_t MJ_j_12;
   Double_t MJJ_j_12;
   Double_t MJ_j_14;
   Double_t MJJ_j_14;
   Double_t MJ_j_16;
   Double_t MJJ_j_16;
   Double_t jetAK8_pt;
   Double_t jetAK8_mass;
   Double_t jetAK8_jec;
   Double_t METraw_et;
   Double_t METraw_phi;
   Double_t METraw_sumEt;
   Double_t MET_et;
   Double_t MET_phi;
   Double_t MET_sumEt;
   Int_t CategoryID;
   Int_t vTagID;//1: tau21<0.45; 0: tau21>0.45 <0.75; -1: tau21 >0.75
   Double_t isMatch;
   Double_t        weight;
   Double_t        weight_nobtag;
   Double_t        IDweight;
   Double_t        IDweightISO;
   Double_t        IDweighttrk;
   Double_t        ToppTweight;
   Double_t        trigger_eff;
   Double_t        btagweight_center;
   Double_t        btagweight_up;
   Double_t        btagweight_down;

   Long64_t           event;
   Int_t           nVtx;
   Int_t           numCands;
   Double_t        ptVlep;
   Double_t        ptVhad;
   Double_t        yVlep;
   Double_t        yVhad;
   Double_t        yVhadJEC;
   Double_t        phiVlep;
   Double_t        phiVhad;
   Double_t        massVlep;
   Double_t        mtVlep;
   Double_t        massVhad;
   Double_t        tau1;
   Double_t        tau2;
   Double_t        tau3;
   Int_t           lep;
   Int_t           channel;
   Double_t        candMass;
   Double_t        ptlep1;
   Double_t        ptlep2;
   Double_t        etalep1;
   Double_t        etalep2;
   Double_t        philep1;
   Double_t        philep2;
   Double_t        met;
   Double_t        metPhi;
   Double_t        theWeight;
   Double_t        nump;
   Double_t        numm;
   Double_t        npT;
   Double_t        npIT;
   Int_t           nBX;
   Double_t        triggerWeight;
   Double_t        lumiWeight;
   Double_t        pileupWeight;
   Double_t        delPhilepmet;
   Double_t        deltaRlepjet;
   Double_t        delPhijetmet;
   Double_t        delPhijetlep;
   Int_t           vbftag;
   Bool_t          IDLoose;
   Bool_t          IDTight;
   Double_t        trackIso;
   Double_t        muchaiso;
   Double_t        muneuiso;
   Double_t        muphoiso;
   Double_t        muPU;
   Double_t        muisolation;
   Double_t        jetAK8_pt1[3];
   Double_t        jetAK8_eta1[3];
   Double_t        jetAK8_mass1[3];
   Double_t        jetAK8_SF_mass1[3];
   Double_t        jetAK8_SF_mass2[3];
   Double_t        jetAK8_jec1[3];
   Double_t        jetAK8_eta;
   Double_t        jetAK8_phi;
   Double_t        candMassJEC;
   Double_t        ptVlepJEC;
   Double_t        yVlepJEC;
   Double_t        phiVlepJEC;
   Double_t        massVlepJEC;
   Double_t        massVhadJEC;
   Double_t        jetAK8puppi_sdJEC;
   Double_t        jetAK8puppi_sd;
   Double_t        jetAK8puppi_tau21;
   Double_t        jetAK8puppi_tau1;
   Double_t        jetAK8puppi_tau2;
   Double_t        jetAK8puppi_tau3;
   Double_t        jetAK8puppi_tau4;
   Double_t        jetAK8puppi_ptJEC;
   Double_t        jetAK8puppi_eta;
   Float_t         Etaj;
   Double_t        jetAK8puppi_phi;
   Float_t         Phij;
   Double_t        jetAK8puppi_sdcorr;
   Double_t        candMasspuppiJEC;

   Double_t        jetAK8puppi_sdJEC_2;
   Double_t        jetAK8puppi_sd_2;
   Double_t        jetAK8puppi_tau21_2;
   Double_t        jetAK8puppi_tau42_2;
   Double_t        jetAK8puppi_tau1_2;
   Double_t        jetAK8puppi_tau2_2;
   Double_t        jetAK8puppi_tau3_2;
   Double_t        jetAK8puppi_tau4_2;
   Double_t        jetAK8puppi_tau42;
   Double_t        jetAK8puppi_ptJEC_2;
   Double_t        jetAK8puppi_eta_2;
   Float_t         Etaj_2;
   Double_t        jetAK8puppi_phi_2;
   Float_t         Phij_2;
   Double_t        jetAK8puppi_sdcorr_2;
Double_t        Mj_2;
Double_t        Mj_un_2;
Double_t        Mj_corr_2;
Double_t        PTj_2;
Float_t         PTj_23;
Float_t         ST;
Float_t         HT;
Float_t         Nj4;
Float_t         Nj8;
Double_t        tau21_2;
Float_t        ak4Ptex1;
Float_t        ak4Etaex1;
Float_t        ak4Phiex1;
Float_t        ak4Eex1;
Float_t        ak4Ptex2;
Float_t        ak4Etaex2;
Float_t        ak4Phiex2;
Float_t        ak4Eex2;
Float_t        Mj_max;
Float_t        Mj_mid; 
Float_t        Mj_min;
Float_t        PTj_max;
Float_t        PTj_mid;
Float_t        PTj_min;
Float_t        Etaj_max;
Float_t        Etaj_mid;
Float_t        Etaj_min;
Float_t        Phij_max;
Float_t        Phij_mid;
Float_t        Phij_min;
Float_t        DPhi_max_mid;
Float_t        DPhi_mid_min;
Float_t        DPhi_min_max;

Float_t        DEta_max_mid;
Float_t        DEta_mid_min;
Float_t        DEta_min_max;
Float_t        DEta_12;
Float_t        DEta_13;
Float_t        DEta_23;
Float_t        DEta_max;
Float_t        DEta_min;

Float_t        DR_max_mid;
Float_t        DR_mid_min;
Float_t        DR_min_max;
Float_t        DR_max;
Float_t        DR_min;
Float_t        MJJj;
Float_t        MJJjj;
Float_t        MJJJj;
Float_t        MJJJJ;
Float_t        t21t31t41_max;
Float_t        t21t31t41_mid;
Float_t        t21t31t41_min;
Float_t        S_t21t31t41_max;
Float_t        S_t21t31t41_mid;
Float_t        S_t21t31t41_min;
Float_t        jet_tau1_puppi_2;
Float_t        jet_tau2_puppi_2;
Float_t        jet_tau3_puppi_2;
Float_t        jet_tau4_puppi_2;
Float_t        tau31_2;
Float_t        tau41_2;
Float_t        tau32_2;
Float_t        tau43_2;
Double_t       tau42;
Double_t       tau42_2;
Float_t      t21t31t41;
Float_t      t21t31t41_2;
Float_t      t21t31t41_3;
Float_t      tau21_max;
Float_t      tau31_max;
Float_t      tau41_max;
Float_t      tau32_max;
Float_t      tau42_max;
Float_t      tau43_max;
Float_t      tau21_mid;
Float_t      tau31_mid;
Float_t      tau41_mid;
Float_t      tau32_mid;
Float_t      tau42_mid;
Float_t      tau43_mid;
Float_t      tau21_min;
Float_t      tau31_min;
Float_t      tau41_min;
Float_t      tau32_min;
Float_t      tau42_min;
Float_t      tau43_min;

Float_t      pt1pt2pt3;


Double_t massww[3];
Double_t        jetAK8puppi_sdJEC_3;
Double_t        jetAK8puppi_sdJC_3;
Double_t        jetAK8puppi_sd_3;
Double_t        jetAK8puppi_tau21_3;
Double_t        jetAK8puppi_tau1_3;
Double_t        jetAK8puppi_tau2_3;
Double_t        jetAK8puppi_tau3_3;
Double_t        jetAK8puppi_tau4_3;
Double_t        jetAK8puppi_tau42_3;
Double_t        jetAK8puppi_ptJEC_3;
Double_t        jetAK8puppi_eta_3;
Float_t         Etaj_3;
Float_t         Etaj_4;
Double_t        jetAK8puppi_phi_3;
Float_t         Phij_3;
Float_t         Phij_4;
Double_t        jetAK8puppi_sdcorr_3;
Double_t        Mj_3;
Double_t        Mj_4;
Double_t        Mj_un_3;
Double_t        Mj_corr_3;
Double_t        PTj_3;
Double_t        tau21_3;
Double_t        tau42_3;

Double_t        jetAK8puppi_eta_4;
Double_t        jetAK8puppi_phi_4;
Double_t        jetAK8puppi_ptJEC_4;
Double_t        jetAK8puppi_ptJEC_5;
Double_t        jetAK8puppi_ptJEC_6;
Double_t        jetAK8puppi_ptJEC_7;
Double_t        jetAK8puppi_ptJEC_8;
Double_t        jetAK8puppi_sdJEC_4;
Double_t        jetAK8puppi_sdJEC_5;
Double_t        jetAK8puppi_sdJEC_6;
Double_t        jetAK8puppi_sdJEC_7;
Double_t        jetAK8puppi_sdJEC_8;
Double_t        gen_gra_m;
Double_t        gen_gra_pt;
Double_t        gen_gra_eta;
Double_t        gen_rad_m;
Double_t        gen_rad_pt;
Double_t        gen_rad_eta;
Double_t        gen_rad_phi;
Double_t        gen_tau_e;
Double_t        gen_tau_pt;
Double_t        gen_tau_eta;
Double_t        gen_tau_phi;
Double_t        gen_tau_e_2;
Double_t        gen_tau_pt_2;
Double_t        gen_tau_eta_2;
Double_t        gen_tau_phi_2;
Double_t        gen_tau_e_3;
Double_t        gen_tau_pt_3;
Double_t        gen_tau_eta_3;
Double_t        gen_tau_phi_3;
Double_t        ptGenVhad;
Double_t        etaGenVhad;
Double_t        phiGenVhad;
Double_t        massGenVhad;
Double_t        ptGenV_2;
Double_t        etaGenV_2;
Double_t        phiGenV_2;
Double_t        massGenV_2;
Double_t        ptGenV_3;
Double_t        etaGenV_3;
Double_t        phiGenV_3;
Double_t        massGenV_3;
Double_t        ptGenVlep;
Double_t        etaGenVlep;
Double_t        phiGenVlep;
Double_t        massGenVlep;
Double_t        ptGenVlep_2;
Double_t        etaGenVlep_2;
Double_t        phiGenVlep_2;
Double_t        massGenVlep_2;
Double_t        ptGenVlep_3;
Double_t        etaGenVlep_3;
Double_t        phiGenVlep_3;
Double_t        massGenVlep_3;
Double_t        ptq11;
Double_t        etaq11;
Double_t        phiq11;
Double_t        massq11;
Double_t        ptq21;
Double_t        etaq21;
Double_t        phiq21;
Double_t        massq21;
Double_t        ptq31;
Double_t        etaq31;
Double_t        phiq31;
Double_t        massq31;
Double_t        ptq12;
Double_t        etaq12;
Double_t        phiq12;
Double_t        massq12;
Double_t        ptq22;
Double_t        etaq22;
Double_t        phiq22;
Double_t        massq22;
Double_t        ptq32;
Double_t        etaq32;
Double_t        phiq32;
Double_t        massq32;
Int_t        status_1;
Int_t        status_2;
Int_t        status_3;
Float_t      PTj_4;
Float_t      PTj_5;
Float_t      PTj_6;
Float_t      PTj_7;
Float_t      PTj_8;
Float_t      newgen_gra_m;
Float_t      newgen_gra_pt;
Float_t      newgen_gra_eta;
Float_t      newgen_rad_m;
Float_t      newgen_rad_pt;
Float_t      newgen_rad_eta;
Float_t      newgen_tau_e;
Float_t      newgen_tau_pt;
Float_t      newgen_tau_eta;
Float_t      newgen_tau_phi;
Float_t      newgen_tau_e_2;
Float_t      newgen_tau_pt_2;
Float_t      newgen_tau_eta_2;
Float_t      newgen_tau_phi_2;
Float_t      newgen_tau_e_3;
Float_t      newgen_tau_pt_3;
Float_t      newgen_tau_eta_3;
Float_t      newgen_tau_phi_3;
Float_t      newptGenVhad;
Float_t      newetaGenVhad;
Float_t      newphiGenVhad;
Float_t      newmassGenVhad;
Float_t      newptGenV_2;
Float_t      newetaGenV_2;
Float_t      newphiGenV_2;
Float_t      newmassGenV_2;
Float_t      newptGenV_3;
Float_t      newetaGenV_3;
Float_t      newphiGenV_3;
Float_t      newmassGenV_3;
Float_t      newptGenVlep;
Float_t      newetaGenVlep;
Float_t      newphiGenVlep;
Float_t      newmassGenVlep;
Float_t      newptGenVlep_2;
Float_t      newetaGenVlep_2;
Float_t      newphiGenVlep_2;
Float_t      newmassGenVlep_2;
Float_t      newptGenVlep_3;
Float_t      newetaGenVlep_3;
Float_t      newphiGenVlep_3;
Float_t      newmassGenVlep_3;
Float_t      newptq11;
Float_t      newetaq11;
Float_t      newphiq11;
Float_t      newmassq11;
Float_t      newptq21;
Float_t      newetaq21;
Float_t      newphiq21;
Float_t      newmassq21; 
Float_t      newptq31;
Float_t      newetaq31;
Float_t      newphiq31;
Float_t      newmassq31;
Float_t      newptq12;
Float_t      newetaq12;
Float_t      newphiq12; 
Float_t      newmassq12;
Float_t      newptq22;
Float_t      newetaq22; 
Float_t      newphiq22;
Float_t      newmassq22;
Float_t      newptq32;
Float_t      newetaq32;
Float_t      newphiq32;
Float_t      newmassq32;
Float_t      Mass2j1j2;
Float_t      Mass2j3j1;
Float_t      Mass2j2j3;
Float_t      MJJ;
Float_t      Mw1w2;
Float_t      Mw1w3;
Float_t      Mw2w3;
Float_t      Pt2dPt1;
Float_t      Pt3dPt1;

Float_t Phij_12;
Float_t Phij_13;
Float_t Phij_23;
Float_t DR_12;
Float_t DR_13;
Float_t DR_23;
Float_t DPhi_max;
Float_t DPhi_min;
Float_t DPhi_max2;
Float_t DPhi_min2;
Float_t DPhi_max3;
Float_t DPhi_min3;
   Float_t jet_tau1_puppi_3;
   Float_t jet_tau2_puppi_3;
   Float_t jet_tau3_puppi_3;
   Float_t jet_tau4_puppi_3;
   Float_t tau31_3;
   Float_t tau41_3;
   Float_t tau32_3;
   Float_t tau43_3;
   Float_t Mj_mean;
   Double_t        mtVlepJEC;
   Int_t           HLT_Ele1;
   Int_t           HLT_Ele2;
   Int_t           HLT_Ele3;
   Int_t           HLT_Ele4;
   Int_t           HLT_Ele5;
   Int_t           HLT_Ele6;
   Int_t           HLT_Ele7;
   Int_t           HLT_Ele8;
   Int_t           HLT_Mu1;
   Int_t           HLT_Mu2;
   Int_t           HLT_Mu3;
   Int_t           HLT_Mu4;
   Int_t           HLT_Mu5;
   Int_t           HLT_Mu6;
   Int_t           HLT_Mu7;
   Int_t           HLT_Mu8;
   Int_t           HLT_Mu9;
   Int_t           HLT_Mu10;
   Int_t           HLT_Mu11;
   Int_t           HLT_Mu12;
   Int_t           HLT_Mu13;
   Int_t           HLT_Mu14;
   Int_t           HLT_Mu15;
   Int_t           HLT_Mu16;
   Int_t num_ak4jetsex;
   Int_t num_ak4jetsin;
   Bool_t          passFilter_HBHE;
   Bool_t          passFilter_HBHEIso;
   Bool_t          passFilter_GlobalHalo;
   Bool_t          passFilter_ECALDeadCell;
   Bool_t          passFilter_GoodVtx;
   Bool_t          passFilter_EEBadSc;
   Bool_t          passFilter_badMuon;
   Bool_t          passFilter_badChargedHadron;
   Int_t           ak4jet_hf[8];
   Int_t           ak4jet_pf[8];
   Double_t        ak4jet_pt[8];
   Double_t        ak4jet_pt_uncorr[8];
   Double_t        ak4jet_eta[8];
   Double_t        ak4jet_phi[8];
   Double_t        ak4jet_e[8];
   Double_t        ak4jet_dr[8];
   Double_t        ak4jet_csv[8];
   Double_t        ak4jet_icsv[8];
   Double_t        deltaRAK4AK8[8];
   Double_t        deltaRAK4AK8_new[8];
   Double_t        ak4jet_IDLoose[8];
   Double_t        ak4jet_IDTight[8];
   Double_t        gentop_pt;
   Double_t        gentop_eta;
   Double_t        gentop_phi;
   Double_t        gentop_mass;
   Double_t        genantitop_pt;
   Double_t        genantitop_eta;
   Double_t        genantitop_phi;
   Double_t        genantitop_mass;
//Subjet
   TLorentzVector  *ak8sj11;
   TLorentzVector  *ak8sj21;
   TLorentzVector  *ak8sj31;
   TLorentzVector  *ak8sj12;
   TLorentzVector  *ak8sj22;
   TLorentzVector  *ak8sj32;
   TLorentzVector  *ak8sj13;
   TLorentzVector  *ak8sj23;
   TLorentzVector  *ak8sj33;
   TLorentzVector  *ak8sj14;
   TLorentzVector  *ak8sj24;
   TLorentzVector  *ak8sj34;
   TLorentzVector  *ak8sj15;
   TLorentzVector  *ak8sj25;
   TLorentzVector  *ak8sj35;
   // List of branches
   TBranch        *b_ak8sj11;   //!
   TBranch        *b_ak8sj21;   //!
   TBranch        *b_ak8sj31;   //!
   TBranch        *b_ak8sj12;   //!
   TBranch        *b_ak8sj22;   //!
   TBranch        *b_ak8sj32;   //!
   TBranch        *b_ak8sj13;   //!
   TBranch        *b_ak8sj23;   //!
   TBranch        *b_ak8sj33;   //!
   TBranch        *b_ak8sj14;   //!
   TBranch        *b_ak8sj24;   //!
   TBranch        *b_ak8sj34;   //!
   TBranch        *b_ak8sj15;   //!
   TBranch        *b_ak8sj25;   //!
   TBranch        *b_ak8sj35;   //!
   TBranch        *b_run;   //!
   TBranch        *b_ls;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_theWeight;   //!
   TBranch        *b_nump;   //!
   TBranch        *b_numm;   //!
   TBranch        *b_npT;   //!
   TBranch        *b_nLooseEle;   //!
   TBranch        *b_nLooseMu;   //!
   TBranch        *b_lep;   //!
   TBranch        *b_ptlep1;   //!
   TBranch        *b_etalep1;   //!
   TBranch        *b_philep1;   //!
   TBranch        *b_trackIso;   //!
   TBranch        *b_muisolation;   //add
   TBranch        *b_muchaiso;   //add
   TBranch        *b_muneuiso;   //add
   TBranch        *b_MET_et;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_ptVlepJEC;   //!
   TBranch        *b_yVlepJEC;   //!
   TBranch        *b_phiVlepJEC;   //!
   TBranch        *b_mtVlepJEC;   //!
   TBranch        *b_massVlepJEC;   //!
   TBranch        *b_jetAK8puppi_sdJEC;   //!
   TBranch        *b_jetAK8puppi_sd;   //!
   TBranch        *b_jetAK8puppi_tau21;   //!
TBranch        *b_jetAK8puppi_tau1;   //!
TBranch        *b_jetAK8puppi_tau2;   //!
TBranch        *b_jetAK8puppi_tau3;   //!
TBranch        *b_jetAK8puppi_tau4;   //!   
TBranch        *b_jetAK8puppi_ptJEC;   //!
   TBranch        *b_jetAK8puppi_eta;   //!
   TBranch        *b_jetAK8puppi_phi;   //!
   TBranch        *b_jetAK8puppi_eta_4;   //!
   TBranch        *b_jetAK8puppi_phi_4;   //!
   TBranch        *b_jetAK8puppi_sdcorr;   //!
   TBranch        *b_IDLoose;   //!
   TBranch        *b_delPhilepmet;   //!
   TBranch        *b_deltaRlepjet;   //!
   TBranch        *b_delPhijetmet;   //!
   TBranch        *b_delPhijetlep;   //!
   TBranch        *b_vbftag;   //add
   TBranch        *b_candMasspuppiJEC;   //!
   TBranch        *b_gentop_pt; //add old
   TBranch        *b_gentop_eta; //add old
   TBranch        *b_genantitop_pt; //add old
   TBranch        *b_genantitop_eta; //add old
   TBranch        *b_HLT_Ele1;
   TBranch        *b_HLT_Ele2;
   TBranch        *b_HLT_Ele3;   //!
   TBranch        *b_HLT_Ele4;
   TBranch        *b_HLT_Ele5;
   TBranch        *b_HLT_Ele6;   //!
   TBranch        *b_HLT_Ele7;
   TBranch        *b_HLT_Ele8;
   TBranch        *b_HLT_Mu1;
   TBranch        *b_HLT_Mu2;   //!
   TBranch        *b_HLT_Mu3;
   TBranch        *b_HLT_Mu4;   //!
   TBranch        *b_HLT_Mu5;
   TBranch        *b_HLT_Mu6;   //!
   TBranch        *b_HLT_Mu7;
   TBranch        *b_HLT_Mu8;   //!
   TBranch        *b_HLT_Mu9;
   TBranch        *b_HLT_Mu10;   //!
   TBranch        *b_HLT_Mu11;
   TBranch        *b_HLT_Mu12;   //!
   TBranch        *b_HLT_Mu13;
   TBranch        *b_HLT_Mu14;   //!
   TBranch        *b_HLT_Mu15;
   TBranch        *b_HLT_Mu16;   //!
   TBranch        *b_ak4jet_hf;
   TBranch        *b_ak4jet_pf;
   TBranch        *b_ak4jet_pt;   //!
   TBranch        *b_ak4jet_pt_uncorr;   //!
   TBranch        *b_ak4jet_eta;   //!
   TBranch        *b_ak4jet_phi;   //!
   TBranch        *b_ak4jet_e;   //!
   TBranch        *b_ak4jet_dr;   //!
   TBranch        *b_ak4jet_csv;   //!
   TBranch        *b_ak4jet_icsv;   //!
   TBranch        *b_deltaRAK4AK8;   //!
   TBranch        *b_ak4jet_IDLoose;   //!
   TBranch        *b_ak4jet_IDTight;   //!
   TBranch        *b_passFilter_HBHE_;   //!
   TBranch        *b_passFilter_HBHEIso_;   //!
   TBranch        *b_passFilter_GlobalHalo_;   //!
   TBranch        *b_passFilter_ECALDeadCell_;   //!
   TBranch        *b_passFilter_GoodVtx_;   //!
   TBranch        *b_passFilter_EEBadSc_;   //!
   TBranch        *b_passFilter_badMuon_;   //!
   TBranch        *b_passFilter_badChargedHadron_;   //!
   TBranch        *b_muphoiso;  //add
   TBranch        *b_muPU;  //add
   TBranch        *b_IDTight;   //!
   TBranch        *b_vbfeta;   //add
   TBranch        *b_vbfmjj;   //add
   TBranch        *b_nj1;   //add
   TBranch        *b_nj2;   //add
   TBranch        *b_ptVlep;   //!
   TBranch        *b_yVlep;   //!
   TBranch        *b_phiVlep;   //!
   TBranch        *b_massVlep;   //!
   TBranch        *b_mtVlep;   //! 
   TBranch        *b_ptVhad;   //!
   TBranch        *b_jetAK8_pt;   //!
   TBranch        *b_yVhad;   //!
   TBranch        *b_yVhadJEC;   //!
   TBranch        *b_phiVhad;   //!
   TBranch        *b_massVhad;   //!
   TBranch        *b_massVhadJEC;   //!
   TBranch        *b_tau1;   //!
   TBranch        *b_tau2;   //!
   TBranch        *b_tau3;   //!
   TBranch        *b_candMass;   //!
   TBranch        *b_numCands;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_ptlep2;   //!
   TBranch        *b_etalep2;   //!
   TBranch        *b_philep2;   //!
   TBranch        *b_met;   //!
   TBranch        *b_metPhi;   //!
   TBranch        *b_npIT;   //!
   TBranch        *b_nBX;   //!
   TBranch        *b_triggerWeight;   //!
   TBranch        *b_lumiWeight;   //!
   TBranch        *b_pileupWeight;   //!
   TBranch        *b_METraw_et;   //!
   TBranch        *b_METraw_phi;   //!
   TBranch        *b_METraw_sumEt;   //!
   TBranch        *b_MET_sumEt;   //!
   TBranch        *b_jetAK8_mass;   //!
   TBranch        *b_jetAK8_jec;   //!
   TBranch        *b_jetAK8_pt1;   //!
   TBranch        *b_jetAK8_eta1;   //!
   TBranch        *b_jetAK8_mass1;   //!
   TBranch        *b_jetAK8_SF_mass1;   //!
   TBranch        *b_jetAK8_SF_mass2;   //!
   TBranch        *b_jetAK8_jec1;   //!
   TBranch        *b_jetAK8_eta;   //!
   TBranch        *b_jetAK8_phi;   //!
   TBranch        *b_candMassJEC;   //!
   

   TBranch        *b_jetAK8puppi_tau42;   //!
   TBranch        *b_jetAK8puppi_sdJEC_2;   //!
   TBranch        *b_jetAK8puppi_sd_2;   //!
   TBranch        *b_jetAK8puppi_tau21_2;   //!
   TBranch        *b_jetAK8puppi_tau1_2;   //!
   TBranch        *b_jetAK8puppi_tau2_2;   //!
   TBranch        *b_jetAK8puppi_tau3_2;   //!
   TBranch        *b_jetAK8puppi_tau4_2;   //!  
   TBranch        *b_jetAK8puppi_tau42_2;   //!
   TBranch        *b_jetAK8puppi_ptJEC_2;   //!
   TBranch        *b_jetAK8puppi_eta_2;   //!
   TBranch        *b_jetAK8puppi_phi_2;   //!
   TBranch        *b_jetAK8puppi_sdcorr_2;   //!
   TBranch        *b_jetAK8puppi_sdJEC_3;   //!
   TBranch         *b_status_1;
   TBranch         *b_status_2;
   TBranch         *b_status_3;
   TBranch        *b_jetAK8puppi_sd_3;   //!
   TBranch        *b_jetAK8puppi_tau21_3;   //!
   TBranch        *b_jetAK8puppi_tau1_3;   //!
   TBranch        *b_jetAK8puppi_tau2_3;   //!
   TBranch        *b_jetAK8puppi_tau3_3;   //!
   TBranch        *b_jetAK8puppi_tau4_3;   //!  
   TBranch        *b_jetAK8puppi_tau42_3;   //!
   TBranch        *b_jetAK8puppi_ptJEC_3;   //!
   TBranch        *b_jetAK8puppi_eta_3;   //!
   TBranch        *b_jetAK8puppi_phi_3;   //!
   TBranch        *b_jetAK8puppi_sdcorr_3;   //!
   TBranch        *b_massww;   //!

   TBranch        *b_gentop_phi;   //add old
   TBranch        *b_gentop_mass;   //add old
   TBranch        *b_genantitop_phi;   //add old
   TBranch        *b_genantitop_mass;   //add old

TBranch  *b_jetAK8puppi_ptJEC_4;
TBranch  *b_jetAK8puppi_ptJEC_5;
TBranch  *b_jetAK8puppi_ptJEC_6;
TBranch  *b_jetAK8puppi_ptJEC_7;
TBranch  *b_jetAK8puppi_ptJEC_8;
TBranch  *b_jetAK8puppi_sdJEC_4;
TBranch  *b_jetAK8puppi_sdJEC_5;
TBranch  *b_jetAK8puppi_sdJEC_6;
TBranch  *b_jetAK8puppi_sdJEC_7;
TBranch  *b_jetAK8puppi_sdJEC_8;
TBranch  *b_gen_gra_m;
TBranch  *b_gen_gra_pt;
TBranch  *b_gen_gra_eta;
TBranch  *b_gen_rad_m;
TBranch  *b_gen_rad_pt;
TBranch  *b_gen_rad_eta;
TBranch  *b_gen_rad_phi;
TBranch  *b_gen_tau_e;
TBranch  *b_gen_tau_pt;
TBranch  *b_gen_tau_eta;
TBranch  *b_gen_tau_phi;
TBranch  *b_gen_tau_e_2;
TBranch  *b_gen_tau_pt_2;
TBranch  *b_gen_tau_eta_2;
TBranch  *b_gen_tau_phi_2;
TBranch  *b_gen_tau_e_3;
TBranch  *b_gen_tau_pt_3;
TBranch  *b_gen_tau_eta_3;
TBranch  *b_gen_tau_phi_3;
TBranch  *b_ptGenVhad;
TBranch  *b_etaGenVhad;
TBranch  *b_phiGenVhad;
TBranch  *b_massGenVhad;
TBranch  *b_ptGenV_2;
TBranch  *b_etaGenV_2;
TBranch  *b_phiGenV_2;
TBranch  *b_massGenV_2;
TBranch  *b_ptGenV_3;
TBranch  *b_etaGenV_3;
TBranch  *b_phiGenV_3;
TBranch  *b_massGenV_3;
TBranch  *b_ptGenVlep;
TBranch  *b_etaGenVlep;
TBranch  *b_phiGenVlep;
TBranch  *b_massGenVlep;
TBranch  *b_ptGenVlep_2;
TBranch  *b_etaGenVlep_2;
TBranch  *b_phiGenVlep_2;
TBranch  *b_massGenVlep_2;
TBranch  *b_ptGenVlep_3;
TBranch  *b_etaGenVlep_3;
TBranch  *b_phiGenVlep_3;
TBranch  *b_massGenVlep_3;
TBranch  *b_ptq11;
TBranch  *b_etaq11;
TBranch  *b_phiq11;
TBranch  *b_massq11;
TBranch  *b_ptq21;
TBranch  *b_etaq21;
TBranch  *b_phiq21;
TBranch  *b_massq21;
TBranch  *b_ptq31;
TBranch  *b_etaq31;
TBranch  *b_phiq31;
TBranch  *b_massq31;
TBranch  *b_ptq12;
TBranch  *b_etaq12;
TBranch  *b_phiq12;
TBranch  *b_massq12;
TBranch  *b_ptq22;
TBranch  *b_etaq22;
TBranch  *b_phiq22;
TBranch  *b_massq22;
TBranch  *b_ptq32;
TBranch  *b_etaq32;
TBranch  *b_phiq32;
TBranch  *b_massq32; 

   TString m_dataset;
   EDBR2PKUTree(TTree *tree=0, TString dataset="");

   virtual ~EDBR2PKUTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString channelname, Double_t XS,TTree *treew, Int_t IsData);// channelname= "mu" or "el"
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     endJob() ;

   private:
   TTree *ExTree;
   TFile *fout; 
   ofstream *file_cutflow;

};

#endif

#ifdef EDBR2PKUTree_cxx
EDBR2PKUTree::EDBR2PKUTree(TTree *tree, TString dataset) : fChain(0) 
{
//// if parameter tree is not specified (or zero), connect the file
//// used to generate this class and read the Tree.
//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("BulkGravWW750.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("BulkGravWW750.root");
//      }
//      TDirectory * dir = (TDirectory*)f->Get("BulkGravWW750.root:/treeDumper");
//      dir->GetObject("EDBRCandidates",tree);
//
//   }
   m_dataset=dataset;
   Init(tree);
}

EDBR2PKUTree::~EDBR2PKUTree()
{
   if (!fChain) return;
   //delete fChain->GetCurrentFile();
}

Int_t EDBR2PKUTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EDBR2PKUTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EDBR2PKUTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
   ak8sj11 = 0;
   ak8sj21 = 0;
   ak8sj31 = 0;
   ak8sj12 = 0;
   ak8sj22 = 0;
   ak8sj32 = 0;
   ak8sj13 = 0;
   ak8sj23 = 0;
   ak8sj33 = 0;
   ak8sj14 = 0;
   ak8sj24 = 0;
   ak8sj34 = 0;
   ak8sj15 = 0;
   ak8sj25 = 0;
   ak8sj35 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fout = new TFile(m_dataset, "RECREATE");
   ExTree = new TTree("PKUTree","PKUTree");
   file_cutflow =new ofstream(m_dataset+"_eventnum.txt");
   ExTree->Branch("num_ak4jetsex", &num_ak4jetsex, "num_ak4jetsex/I");
   ExTree->Branch("num_ak4jetsin", &num_ak4jetsin, "num_ak4jetsin/I");
   ExTree->Branch("Mj_mean", &Mj_mean, "Mj_mean/F");
   ExTree->Branch("CategoryID", &CategoryID, "CategoryID/I");
   ExTree->Branch("vTagID", &vTagID, "vTagID/I");
   ExTree->Branch("Pt2dPt1", &Pt2dPt1, "Pt2dPt1/F");
   ExTree->Branch("Pt3dPt1", &Pt3dPt1, "Pt3dPt1/F");
   ExTree->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
   ExTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/D");

   ExTree->Branch("Mass2j1j2", &Mass2j1j2, "Mass2j1j2/F");
   ExTree->Branch("Mass2j3j1", &Mass2j3j1, "Mass2j3j1/F");
   ExTree->Branch("Mass2j2j3", &Mass2j2j3, "Mass2j2j3/F");
   ExTree->Branch("MJ_j_18", &MJ_j_18, "MJ_j_18/D");
   ExTree->Branch("MJJ_j_18", &MJJ_j_18, "MJJ_j_18/D");
   ExTree->Branch("MJ_j_10", &MJ_j_10, "MJ_j_10/D");
   ExTree->Branch("MJJ_j_10", &MJJ_j_10, "MJJ_j_10/D");
   ExTree->Branch("MJ_j_12", &MJ_j_12, "MJ_j_12/D");
   ExTree->Branch("MJJ_j_12", &MJJ_j_12, "MJJ_j_12/D");
   ExTree->Branch("MJ_j_14", &MJ_j_14, "MJ_j_14/D");
   ExTree->Branch("MJJ_j_14", &MJJ_j_14, "MJJ_j_14/D");  
   ExTree->Branch("MJ_j_16", &MJ_j_16, "MJ_j_16/D");
   ExTree->Branch("MJJ_j_16", &MJJ_j_16, "MJJ_j_16/D");
   ExTree->Branch("MJJ", &MJJ, "MJJ/F");
   ExTree->Branch("Mw1w2", &Mw1w2, "Mw1w2/F");
   ExTree->Branch("Mw1w3", &Mw1w3, "Mw1w3/F");
   ExTree->Branch("Mw2w3", &Mw2w3, "Mw2w3/F");
   ExTree->Branch("event",&event,"event/L");
   ExTree->Branch("lumi",&ls,"lumi/I");
   ExTree->Branch("nPV",&nVtx,"nPV/I");
   ExTree->Branch("pfMET",&pfMET,"pfMET/F");
   ExTree->Branch("pfMETPhi",&pfMETPhi,"pfMETPhi/F");
   ExTree->Branch("weight", &weight, "weight/D");
   ExTree->Branch("weight_nobtag", &weight_nobtag, "weight_nobtag/D");
   ExTree->Branch("isMatch", &isMatch, "isMatch/D");

   ExTree->Branch("ak4Ptex1",&ak4Ptex1,"ak4Ptex1/F");
   ExTree->Branch("ak4Etaex1",&ak4Etaex1,"ak4Etaex1/F");
   ExTree->Branch("ak4Phiex1",&ak4Phiex1,"ak4Phiex1/F");
   ExTree->Branch("ak4Eex1",&ak4Eex1,"ak4Eex1/F");
   ExTree->Branch("ak4Ptex2",&ak4Ptex2,"ak4Ptex2/F");
   ExTree->Branch("ak4Etaex2",&ak4Etaex2,"ak4Etaex2/F");
   ExTree->Branch("ak4Phiex2",&ak4Phiex2,"ak4Phiex2/F");
   ExTree->Branch("ak4Eex2",&ak4Eex2,"ak4Eex2/F");
   ExTree->Branch("Mj",&Mj,"Mj/F");
  ExTree->Branch("Mj_max",&Mj_max,"Mj_max/F");  
  ExTree->Branch("Mj_mid",&Mj_mid,"Mj_mid/F");
  ExTree->Branch("Mj_min",&Mj_min,"Mj_min/F");
  ExTree->Branch("PTj_max",&PTj_max,"PTj_max/F");
  ExTree->Branch("PTj_mid",&PTj_mid,"PTj_mid/F");
  ExTree->Branch("PTj_min",&PTj_min,"PTj_min/F");
  ExTree->Branch("Etaj_max",&Etaj_max,"Etaj_max/F");
  ExTree->Branch("Etaj_mid",&Etaj_mid,"Etaj_mid/F");
  ExTree->Branch("Etaj_min",&Etaj_min,"Etaj_min/F");
  ExTree->Branch("Phij_max",&Phij_max,"Phij_max/F");
  ExTree->Branch("Phij_mid",&Phij_mid,"Phij_mid/F");
  ExTree->Branch("Phij_min",&Phij_min,"Phij_min/F");
  ExTree->Branch("DPhi_max_mid",&DPhi_max_mid,"DPhi_max_mid/F");
  ExTree->Branch("DPhi_mid_min",&DPhi_mid_min,"DPhi_mid_min/F");
  ExTree->Branch("DPhi_min_max",&DPhi_min_max,"DPhi_min_max/F");
  ExTree->Branch("DEta_max_mid",&DEta_max_mid,"DEta_max_mid/F");
  ExTree->Branch("DEta_mid_min",&DEta_mid_min,"DEta_mid_min/F");
  ExTree->Branch("DEta_min_max",&DEta_min_max,"DEta_min_max/F");
  ExTree->Branch("DEta_12",&DEta_12,"DEta_12/F");
  ExTree->Branch("DEta_13",&DEta_13,"DEta_13/F");
  ExTree->Branch("DEta_23",&DEta_23,"DEta_23/F");
  ExTree->Branch("DEta_max",&DEta_max,"DEta_max/F");
  ExTree->Branch("DEta_min",&DEta_min,"DEta_min/F");
  ExTree->Branch("DR_max_mid",&DR_max_mid,"DR_max_mid/F");
  ExTree->Branch("DR_mid_min",&DR_mid_min,"DR_mid_min/F");
  ExTree->Branch("DR_min_max",&DR_min_max,"DR_min_max/F");
  ExTree->Branch("DR_min",&DR_min,"DR_min/F");
  ExTree->Branch("DR_max",&DR_max,"DR_max/F");
  ExTree->Branch("MJJj",&MJJj,"MJJj/F");
  ExTree->Branch("MJJjj",&MJJjj,"MJJjj/F");
  ExTree->Branch("MJJJj",&MJJJj,"MJJJj/F");
  ExTree->Branch("MJJJJ",&MJJJJ,"MJJJJ/F");
  ExTree->Branch("t21t31t41_max",&t21t31t41_max,"t21t31t41_max/F");
  ExTree->Branch("t21t31t41_min",&t21t31t41_min,"t21t31t41_min/F");
  ExTree->Branch("t21t31t41_mid",&t21t31t41_mid,"t21t31t41_mid/F");
  ExTree->Branch("S_t21t31t41_max",&S_t21t31t41_max,"S_t21t31t41_max/F");
  ExTree->Branch("S_t21t31t41_min",&S_t21t31t41_min,"S_t21t31t41_min/F");
  ExTree->Branch("S_t21t31t41_mid",&S_t21t31t41_mid,"S_t21t31t41_mid/F");
  ExTree->Branch("PTj",&PTj,"PTj/F");
  ExTree->Branch("PTj_23",&PTj_23,"PTj_23/F");   
ExTree->Branch("ST",&ST,"ST/F");
ExTree->Branch("HT",&HT,"HT/F");
ExTree->Branch("Nj4",&Nj4,"Nj4/F");
ExTree->Branch("Nj8",&Nj8,"Nj8/F");
ExTree->Branch("tau21",&tau21,"tau21/F");
ExTree->Branch("tau31",&tau31,"tau31/F");
ExTree->Branch("tau41",&tau41,"tau41/F");
ExTree->Branch("tau32",&tau32,"tau32/F");
ExTree->Branch("tau43",&tau43,"tau43/F");
ExTree->Branch("pt1pt2pt3",&pt1pt2pt3,"pt1pt2pt3/F");
ExTree->Branch("t21t31t41",&t21t31t41,"t21t31t41/F");
ExTree->Branch("t21t31t41_2",&t21t31t41_2,"t21t31t41_2/F");
ExTree->Branch("t21t31t41_3",&t21t31t41_3,"t21t31t41_3/F");
   ExTree->Branch("Etaj",&Etaj,"Etaj/F");
   ExTree->Branch("Phij",&Phij,"Phij/F");
   ExTree->Branch("Mj_2",&Mj_2,"Mj_2/D");
   ExTree->Branch("PTj_2",&PTj_2,"PTj_2/D");
   ExTree->Branch("tau21_2",&tau21_2,"tau21_2/D");
ExTree->Branch("tau31_2",&tau31_2,"tau31_2/F");
ExTree->Branch("tau41_2",&tau41_2,"tau41_2/F");
ExTree->Branch("tau32_2",&tau32_2,"tau32_2/F");
ExTree->Branch("tau43_2",&tau43_2,"tau43_2/F");
   ExTree->Branch("tau42",&tau42,"tau42/D");
   ExTree->Branch("tau42_2",&tau42_2,"tau42_2/D");
   ExTree->Branch("Etaj_2",&Etaj_2,"Etaj_2/F");
   ExTree->Branch("Phij_2",&Phij_2,"Phij_2/F");
   ExTree->Branch("Etaj_4",&Etaj_4,"Etaj_4/F");
   ExTree->Branch("Phij_4",&Phij_4,"Phij_4/F");
   ExTree->Branch("Mj_3",&Mj_3,"Mj_3/D");
   ExTree->Branch("Mj_4",&Mj_4,"Mj_4/D");
   ExTree->Branch("PTj_3",&PTj_3,"PTj_3/D");
   ExTree->Branch("tau21_3",&tau21_3,"tau21_3/D");
   ExTree->Branch("tau42_3",&tau42_3,"tau42_3/D");
ExTree->Branch("tau31_3",&tau31_3,"tau31_3/F");
ExTree->Branch("tau41_3",&tau41_3,"tau41_3/F");
ExTree->Branch("tau32_3",&tau32_3,"tau32_3/F");
ExTree->Branch("tau43_3",&tau43_3,"tau43_3/F");
ExTree->Branch("tau21_max",&tau21_max,"tau21_max/F");
ExTree->Branch("tau31_max",&tau31_max,"tau31_max/F");
ExTree->Branch("tau41_max",&tau41_max,"tau41_max/F");
ExTree->Branch("tau21_mid",&tau21_mid,"tau21_mid/F");
ExTree->Branch("tau31_mid",&tau31_mid,"tau31_mid/F");
ExTree->Branch("tau41_mid",&tau41_mid,"tau41_mid/F");
ExTree->Branch("tau21_min",&tau21_min,"tau21_min/F");
ExTree->Branch("tau31_min",&tau31_min,"tau31_min/F");
ExTree->Branch("tau41_min",&tau41_min,"tau41_min/F");
ExTree->Branch("tau32_max",&tau32_max,"tau32_max/F");
ExTree->Branch("tau42_max",&tau42_max,"tau42_max/F");
ExTree->Branch("tau43_max",&tau43_max,"tau43_max/F");
ExTree->Branch("tau32_mid",&tau32_mid,"tau32_mid/F");
ExTree->Branch("tau42_mid",&tau42_mid,"tau42_mid/F");
ExTree->Branch("tau43_mid",&tau43_mid,"tau43_mid/F");
ExTree->Branch("tau32_min",&tau32_min,"tau32_min/F");
ExTree->Branch("tau42_min",&tau42_min,"tau42_min/F");
ExTree->Branch("tau43_min",&tau43_min,"tau43_min/F");
ExTree->Branch("Etaj_3",&Etaj_3,"Etaj_3/F");
ExTree->Branch("Phij_3",&Phij_3,"Phij_3/F");
ExTree->Branch("Phij_13",&Phij_13,"Phij_13/F");
ExTree->Branch("Phij_12",&Phij_12,"Phij_12/F");
ExTree->Branch("Phij_23",&Phij_23,"Phij_23/F");
ExTree->Branch("DPhi_max",&DPhi_max,"DPhi_max/F");
ExTree->Branch("DPhi_min",&DPhi_min,"DPhi_min/F");
ExTree->Branch("DPhi_max2",&DPhi_max2,"DPhi_max2/F");
ExTree->Branch("DPhi_min2",&DPhi_min2,"DPhi_min2/F");
ExTree->Branch("DPhi_max3",&DPhi_max3,"DPhi_max3/F");
ExTree->Branch("DPhi_min3",&DPhi_min3,"DPhi_min3/F");
ExTree->Branch("DR_12",&DR_12,"DR_12/F");
ExTree->Branch("DR_13",&DR_13,"DR_13/F");
ExTree->Branch("DR_23",&DR_23,"DR_23/F");
   ExTree->Branch("MJJJ",&MJJJ,"MJJJ/F");
   ExTree->Branch("nbtag",&nbtag,"nbtag/I");
   ExTree->Branch("num_bJet",&num_bJet,"num_bJet/I");
   ExTree->Branch("num_bJet_loose",&num_bJet_loose,"num_bJet_loose/I");
   ExTree->Branch("num_bJet_tight",&num_bJet_tight,"num_bJet_tight/I");
   ExTree->Branch("deltaRAK4AK8_new",deltaRAK4AK8_new,"deltaRAK4AK8_new[8]/D");
   ExTree->Branch("METraw_et",&METraw_et,"METraw_et/D");
   ExTree->Branch("METraw_phi",&METraw_phi,"METraw_phi/D");
   ExTree->Branch("METraw_sumEt",&METraw_sumEt,"METraw_sumEt/D");
   ExTree->Branch("MET_et",&MET_et,"MET_et/D");
   ExTree->Branch("MET_phi",&MET_phi,"MET_phi/D");
   ExTree->Branch("MET_sumEt",&MET_sumEt,"MET_sumEt/D");
   ExTree->Branch("muisolation",&muisolation,"muisolation/I");
ExTree->Branch("newgen_gra_m",&newgen_gra_m,"newgen_gra_m/F");
ExTree->Branch("newgen_gra_eta",&newgen_gra_eta,"newgen_gra_eta/F");
ExTree->Branch("newgen_gra_pt",&newgen_gra_pt,"newgen_gra_pt/F");
ExTree->Branch("newgen_rad_m",&newgen_rad_m,"newgen_rad_m/F");
ExTree->Branch("newgen_rad_eta",&newgen_rad_eta,"newgen_rad_eta/F");
ExTree->Branch("newgen_rad_pt",&newgen_rad_pt,"newgen_rad_pt/F");
ExTree->Branch("newgen_tau_e",&newgen_tau_e,"newgen_tau_e/F");
ExTree->Branch("newgen_tau_eta",&newgen_tau_eta,"newgen_tau_eta/F");
ExTree->Branch("newgen_tau_pt",&newgen_tau_pt,"newgen_tau_pt/F");
ExTree->Branch("newgen_tau_phi",&newgen_tau_phi,"newgen_tau_phi/F");
ExTree->Branch("newgen_tau_e_2",&newgen_tau_e_2,"newgen_tau_e_2/F");
ExTree->Branch("newgen_tau_eta_2",&newgen_tau_eta_2,"newgen_tau_eta_2/F");
ExTree->Branch("newgen_tau_pt_2",&newgen_tau_pt_2,"newgen_tau_pt_2/F");
ExTree->Branch("newgen_tau_phi_2",&newgen_tau_phi_2,"newgen_tau_phi_2/F");
ExTree->Branch("newgen_tau_e_3",&newgen_tau_e_3,"newgen_tau_e_3/F");
ExTree->Branch("newgen_tau_eta_3",&newgen_tau_eta_3,"newgen_tau_eta_3/F");
ExTree->Branch("newgen_tau_pt_3",&newgen_tau_pt_3,"newgen_tau_pt_3/F");
ExTree->Branch("newgen_tau_phi_3",&newgen_tau_phi_3,"newgen_tau_phi_3/F");
ExTree->Branch("newptGenVhad",&newptGenVhad,"newptGenVhad/F");
ExTree->Branch("newetaGenVhad",&newetaGenVhad,"newetaGenVhad/F");
ExTree->Branch("newphiGenVhad",&newphiGenVhad,"newphiGenVhad/F");
ExTree->Branch("newmassGenVhad",&newmassGenVhad,"newmassGenVhad/F");
ExTree->Branch("newptGenVlep",&newptGenVlep,"newptGenVlep/F");
ExTree->Branch("newetaGenVlep",&newetaGenVlep,"newetaGenVlep/F");
ExTree->Branch("newphiGenVlep",&newphiGenVlep,"newphiGenVlep/F");
ExTree->Branch("newmassGenVlep",&newmassGenVlep,"newmassGenVlep/F");
ExTree->Branch("newptGenV_2",&newptGenV_2,"newptGenV_2/F");
ExTree->Branch("newetaGenV_2",&newetaGenV_2,"newetaGenV_2/F");
ExTree->Branch("newphiGenV_2",&newphiGenV_2,"newphiGenV_2/F");
ExTree->Branch("newmassGenV_2",&newmassGenV_2,"newmassGenV_2/F");
ExTree->Branch("newptGenVlep_2",&newptGenVlep_2,"newptGenVlep_2/F");
ExTree->Branch("newetaGenVlep_2",&newetaGenVlep_2,"newetaGenVlep_2/F");
ExTree->Branch("newphiGenVlep_2",&newphiGenVlep_2,"newphiGenVlep_2/F");
ExTree->Branch("newmassGenVlep_2",&newmassGenVlep_2,"newmassGenVlep_2/F");
ExTree->Branch("newptGenV_3",&newptGenV_3,"newptGenV_3/F");
ExTree->Branch("newetaGenV_3",&newetaGenV_3,"newetaGenV_3/F");
ExTree->Branch("newphiGenV_3",&newphiGenV_3,"newphiGenV_3/F");
ExTree->Branch("newmassGenV_3",&newmassGenV_3,"newmassGenV_3/F");
ExTree->Branch("newptGenVlep_3",&newptGenVlep_3,"newptGenVlep_3/F");
ExTree->Branch("newetaGenVlep_3",&newetaGenVlep_3,"newetaGenVlep_3/F");
ExTree->Branch("newphiGenVlep_3",&newphiGenVlep_3,"newphiGenVlep_3/F");
ExTree->Branch("newmassGenVlep_3",&newmassGenVlep_3,"newmassGenVlep_3/F");
ExTree->Branch("newptq11",&newptq11,"newptq11/F");
ExTree->Branch("newetaq11",&newetaq11,"newetaq11/F");
ExTree->Branch("newphiq11",&newphiq11,"newphiq11/F");
ExTree->Branch("newmassq11",&newmassq11,"newmassq11/F");
ExTree->Branch("newptq12",&newptq12,"newptq12/F");
ExTree->Branch("newetaq12",&newetaq12,"newetaq12/F");
ExTree->Branch("newphiq12",&newphiq12,"newphiq12/F");
ExTree->Branch("newmassq12",&newmassq12,"newmassq12/F");
ExTree->Branch("newptq21",&newptq21,"newptq21/F");
ExTree->Branch("newetaq21",&newetaq21,"newetaq21/F");
ExTree->Branch("newphiq21",&newphiq21,"newphiq21/F");
ExTree->Branch("newmassq21",&newmassq21,"newmassq21/F");
ExTree->Branch("newptq22",&newptq22,"newptq22/F");
ExTree->Branch("newetaq22",&newetaq22,"newetaq22/F");
ExTree->Branch("newphiq22",&newphiq22,"newphiq22/F");
ExTree->Branch("newmassq22",&newmassq22,"newmassq22/F");
ExTree->Branch("newptq31",&newptq31,"newptq31/F");
ExTree->Branch("newetaq31",&newetaq31,"newetaq31/F");
ExTree->Branch("newphiq31",&newphiq31,"newphiq31/F");
ExTree->Branch("newmassq31",&newmassq31,"newmassq31/F");
ExTree->Branch("newptq32",&newptq32,"newptq32/F");
ExTree->Branch("newetaq32",&newetaq32,"newetaq32/F");
ExTree->Branch("newphiq32",&newphiq32,"newphiq32/F");
ExTree->Branch("newmassq32",&newmassq32,"newmassq32/F");
   fChain->SetBranchAddress("ak8sj11", &ak8sj11, &b_ak8sj11);
   fChain->SetBranchAddress("ak8sj21", &ak8sj21, &b_ak8sj21);
   fChain->SetBranchAddress("ak8sj31", &ak8sj31, &b_ak8sj31);
   fChain->SetBranchAddress("ak8sj12", &ak8sj12, &b_ak8sj12);
   fChain->SetBranchAddress("ak8sj22", &ak8sj22, &b_ak8sj22);
   fChain->SetBranchAddress("ak8sj32", &ak8sj32, &b_ak8sj32);
   fChain->SetBranchAddress("ak8sj13", &ak8sj13, &b_ak8sj13);
   fChain->SetBranchAddress("ak8sj23", &ak8sj23, &b_ak8sj23);
   fChain->SetBranchAddress("ak8sj33", &ak8sj33, &b_ak8sj33);
   fChain->SetBranchAddress("ak8sj14", &ak8sj14, &b_ak8sj14);
   fChain->SetBranchAddress("ak8sj24", &ak8sj24, &b_ak8sj24);
   fChain->SetBranchAddress("ak8sj34", &ak8sj34, &b_ak8sj34);
   fChain->SetBranchAddress("ak8sj15", &ak8sj15, &b_ak8sj15);
   fChain->SetBranchAddress("ak8sj25", &ak8sj25, &b_ak8sj25);
   fChain->SetBranchAddress("ak8sj35", &ak8sj35, &b_ak8sj35);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("ls", &ls, &b_ls);
   fChain->SetBranchAddress("nLooseEle", &nLooseEle, &b_nLooseEle);
   fChain->SetBranchAddress("nLooseMu", &nLooseMu, &b_nLooseMu);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("ptVhad", &ptVhad, &b_ptVhad);
   fChain->SetBranchAddress("yVhad", &yVhad, &b_yVhad);
   fChain->SetBranchAddress("yVhadJEC", &yVhadJEC, &b_yVhadJEC);
   fChain->SetBranchAddress("phiVhad", &phiVhad, &b_phiVhad);
   fChain->SetBranchAddress("massVhad", &massVhad, &b_massVhad);
   fChain->SetBranchAddress("tau1", &tau1, &b_tau1);
   fChain->SetBranchAddress("tau2", &tau2, &b_tau2);
   fChain->SetBranchAddress("tau3", &tau3, &b_tau3);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("candMass", &candMass, &b_candMass);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("theWeight", &theWeight, &b_theWeight);
   fChain->SetBranchAddress("nump", &nump, &b_nump);
   fChain->SetBranchAddress("numm", &numm, &b_numm);
   fChain->SetBranchAddress("npT", &npT, &b_npT);
   fChain->SetBranchAddress("npIT", &npIT, &b_npIT);
   fChain->SetBranchAddress("nBX", &nBX, &b_nBX);
   fChain->SetBranchAddress("triggerWeight", &triggerWeight, &b_triggerWeight);
   fChain->SetBranchAddress("lumiWeight", &lumiWeight, &b_lumiWeight);
   fChain->SetBranchAddress("pileupWeight", &pileupWeight, &b_pileupWeight);
   fChain->SetBranchAddress("delPhilepmet", &delPhilepmet, &b_delPhilepmet);
   fChain->SetBranchAddress("deltaRlepjet", &deltaRlepjet, &b_deltaRlepjet);
   fChain->SetBranchAddress("delPhijetmet", &delPhijetmet, &b_delPhijetmet);
   fChain->SetBranchAddress("delPhijetlep", &delPhijetlep, &b_delPhijetlep);
   fChain->SetBranchAddress("vbftag", &vbftag, &b_vbftag);
   fChain->SetBranchAddress("trackIso", &trackIso, &b_trackIso);
   fChain->SetBranchAddress("METraw_et", &METraw_et, &b_METraw_et);
   fChain->SetBranchAddress("METraw_phi", &METraw_phi, &b_METraw_phi);
   fChain->SetBranchAddress("METraw_sumEt", &METraw_sumEt, &b_METraw_sumEt);
   fChain->SetBranchAddress("MET_et", &MET_et, &b_MET_et);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   fChain->SetBranchAddress("jetAK8_pt", &jetAK8_pt, &b_jetAK8_pt);
   fChain->SetBranchAddress("jetAK8_mass", &jetAK8_mass, &b_jetAK8_mass);
   fChain->SetBranchAddress("jetAK8_jec", &jetAK8_jec, &b_jetAK8_jec);
   fChain->SetBranchAddress("jetAK8_pt1", jetAK8_pt1, &b_jetAK8_pt1);
   fChain->SetBranchAddress("jetAK8_eta1", jetAK8_eta1, &b_jetAK8_eta1);
   fChain->SetBranchAddress("jetAK8_mass1", jetAK8_mass1, &b_jetAK8_mass1);
   fChain->SetBranchAddress("jetAK8_SF_mass1", jetAK8_SF_mass1, &b_jetAK8_SF_mass1);
   fChain->SetBranchAddress("jetAK8_SF_mass2", jetAK8_SF_mass2, &b_jetAK8_SF_mass2);
   fChain->SetBranchAddress("jetAK8_jec1", jetAK8_jec1, &b_jetAK8_jec1);
   fChain->SetBranchAddress("jetAK8_eta", &jetAK8_eta, &b_jetAK8_eta);
   fChain->SetBranchAddress("jetAK8_phi", &jetAK8_phi, &b_jetAK8_phi);
   fChain->SetBranchAddress("candMassJEC", &candMassJEC, &b_candMassJEC);
   fChain->SetBranchAddress("massVhadJEC", &massVhadJEC, &b_massVhadJEC);
   fChain->SetBranchAddress("jetAK8puppi_sdJEC", &jetAK8puppi_sdJEC, &b_jetAK8puppi_sdJEC);
   fChain->SetBranchAddress("jetAK8puppi_sd", &jetAK8puppi_sd, &b_jetAK8puppi_sd);
   fChain->SetBranchAddress("jetAK8puppi_tau21", &jetAK8puppi_tau21, &b_jetAK8puppi_tau21);
fChain->SetBranchAddress("jetAK8puppi_tau1", &jetAK8puppi_tau1, &b_jetAK8puppi_tau1);
fChain->SetBranchAddress("jetAK8puppi_tau2", &jetAK8puppi_tau2, &b_jetAK8puppi_tau2);
fChain->SetBranchAddress("jetAK8puppi_tau3", &jetAK8puppi_tau3, &b_jetAK8puppi_tau3);
fChain->SetBranchAddress("jetAK8puppi_tau4", &jetAK8puppi_tau4, &b_jetAK8puppi_tau4);
   fChain->SetBranchAddress("jetAK8puppi_ptJEC", &jetAK8puppi_ptJEC, &b_jetAK8puppi_ptJEC);
   fChain->SetBranchAddress("jetAK8puppi_eta", &jetAK8puppi_eta, &b_jetAK8puppi_eta);
   fChain->SetBranchAddress("jetAK8puppi_phi", &jetAK8puppi_phi, &b_jetAK8puppi_phi);
   fChain->SetBranchAddress("jetAK8puppi_sdcorr", &jetAK8puppi_sdcorr, &b_jetAK8puppi_sdcorr);
   fChain->SetBranchAddress("candMasspuppiJEC", &candMasspuppiJEC, &b_candMasspuppiJEC);

   fChain->SetBranchAddress("jetAK8puppi_sdJEC_2", &jetAK8puppi_sdJEC_2, &b_jetAK8puppi_sdJEC_2);
   fChain->SetBranchAddress("jetAK8puppi_sd_2", &jetAK8puppi_sd_2, &b_jetAK8puppi_sd_2);
   fChain->SetBranchAddress("jetAK8puppi_tau21_2", &jetAK8puppi_tau21_2, &b_jetAK8puppi_tau21_2);
fChain->SetBranchAddress("jetAK8puppi_tau1_2", &jetAK8puppi_tau1_2, &b_jetAK8puppi_tau1_2);
fChain->SetBranchAddress("jetAK8puppi_tau2_2", &jetAK8puppi_tau2_2, &b_jetAK8puppi_tau2_2);
fChain->SetBranchAddress("jetAK8puppi_tau3_2", &jetAK8puppi_tau3_2, &b_jetAK8puppi_tau3_2);
fChain->SetBranchAddress("jetAK8puppi_tau4_2", &jetAK8puppi_tau4_2, &b_jetAK8puppi_tau4_2);
fChain->SetBranchAddress("jetAK8puppi_tau42_2", &jetAK8puppi_tau42_2, &b_jetAK8puppi_tau42_2);
fChain->SetBranchAddress("jetAK8puppi_tau42", &jetAK8puppi_tau42, &b_jetAK8puppi_tau42);
   fChain->SetBranchAddress("jetAK8puppi_ptJEC_2", &jetAK8puppi_ptJEC_2, &b_jetAK8puppi_ptJEC_2);
   fChain->SetBranchAddress("jetAK8puppi_eta_2", &jetAK8puppi_eta_2, &b_jetAK8puppi_eta_2);
   fChain->SetBranchAddress("jetAK8puppi_phi_2", &jetAK8puppi_phi_2, &b_jetAK8puppi_phi_2);
   fChain->SetBranchAddress("jetAK8puppi_sdcorr_2", &jetAK8puppi_sdcorr_2, &b_jetAK8puppi_sdcorr_2);
fChain->SetBranchAddress("status_1",&status_1 ,&b_status_1);
fChain->SetBranchAddress("status_2",&status_2 ,&b_status_2);
fChain->SetBranchAddress("status_3",&status_3 ,&b_status_3);
fChain->SetBranchAddress("jetAK8puppi_sdJEC_3", &jetAK8puppi_sdJEC_3, &b_jetAK8puppi_sdJEC_3);
   fChain->SetBranchAddress("jetAK8puppi_sd_3", &jetAK8puppi_sd_3, &b_jetAK8puppi_sd_3);
   fChain->SetBranchAddress("jetAK8puppi_tau21_3", &jetAK8puppi_tau21_3, &b_jetAK8puppi_tau21_3);
fChain->SetBranchAddress("jetAK8puppi_tau1_3", &jetAK8puppi_tau1_3, &b_jetAK8puppi_tau1_3);
fChain->SetBranchAddress("jetAK8puppi_tau2_3", &jetAK8puppi_tau2_3, &b_jetAK8puppi_tau2_3);
fChain->SetBranchAddress("jetAK8puppi_tau3_3", &jetAK8puppi_tau3_3, &b_jetAK8puppi_tau3_3);
fChain->SetBranchAddress("jetAK8puppi_tau4_3", &jetAK8puppi_tau4_3, &b_jetAK8puppi_tau4_3);
   fChain->SetBranchAddress("jetAK8puppi_tau42_3", &jetAK8puppi_tau42_3, &b_jetAK8puppi_tau42_3);
   fChain->SetBranchAddress("jetAK8puppi_ptJEC_3", &jetAK8puppi_ptJEC_3, &b_jetAK8puppi_ptJEC_3);
   fChain->SetBranchAddress("jetAK8puppi_eta_3", &jetAK8puppi_eta_3, &b_jetAK8puppi_eta_3);
   fChain->SetBranchAddress("jetAK8puppi_phi_3", &jetAK8puppi_phi_3, &b_jetAK8puppi_phi_3);
   fChain->SetBranchAddress("jetAK8puppi_sdcorr_3", &jetAK8puppi_sdcorr_3, &b_jetAK8puppi_sdcorr_3);
   fChain->SetBranchAddress("jetAK8puppi_eta_4", &jetAK8puppi_eta_4, &b_jetAK8puppi_eta_4);
   fChain->SetBranchAddress("jetAK8puppi_phi_4", &jetAK8puppi_phi_4, &b_jetAK8puppi_phi_4);
   fChain->SetBranchAddress("massww",&massww);
   fChain->SetBranchAddress("HLT_Mu1", &HLT_Mu1, &b_HLT_Mu1);
   fChain->SetBranchAddress("HLT_Mu2", &HLT_Mu2, &b_HLT_Mu2);
   fChain->SetBranchAddress("HLT_Mu3", &HLT_Mu3, &b_HLT_Mu3);
   fChain->SetBranchAddress("HLT_Mu4", &HLT_Mu4, &b_HLT_Mu4);
   fChain->SetBranchAddress("HLT_Mu5", &HLT_Mu5, &b_HLT_Mu5);
   fChain->SetBranchAddress("HLT_Mu6", &HLT_Mu6, &b_HLT_Mu6);
   fChain->SetBranchAddress("HLT_Mu7", &HLT_Mu7, &b_HLT_Mu7);
   fChain->SetBranchAddress("HLT_Mu8", &HLT_Mu8, &b_HLT_Mu8);
   fChain->SetBranchAddress("HLT_Mu9", &HLT_Mu9, &b_HLT_Mu9);
   fChain->SetBranchAddress("HLT_Mu10", &HLT_Mu10, &b_HLT_Mu10);
   fChain->SetBranchAddress("HLT_Mu11", &HLT_Mu11, &b_HLT_Mu11);
   fChain->SetBranchAddress("HLT_Mu12", &HLT_Mu12, &b_HLT_Mu12);
   fChain->SetBranchAddress("HLT_Mu13", &HLT_Mu13, &b_HLT_Mu13);
   fChain->SetBranchAddress("HLT_Mu14", &HLT_Mu14, &b_HLT_Mu14);
   fChain->SetBranchAddress("HLT_Mu15", &HLT_Mu15, &b_HLT_Mu15);
   fChain->SetBranchAddress("HLT_Mu16", &HLT_Mu16, &b_HLT_Mu16);
   fChain->SetBranchAddress("passFilter_HBHE", &passFilter_HBHE, &b_passFilter_HBHE_);
   fChain->SetBranchAddress("passFilter_HBHEIso", &passFilter_HBHEIso, &b_passFilter_HBHEIso_);
   fChain->SetBranchAddress("passFilter_GlobalHalo", &passFilter_GlobalHalo, &b_passFilter_GlobalHalo_);
   fChain->SetBranchAddress("passFilter_ECALDeadCell", &passFilter_ECALDeadCell, &b_passFilter_ECALDeadCell_);
   fChain->SetBranchAddress("passFilter_GoodVtx", &passFilter_GoodVtx, &b_passFilter_GoodVtx_);
   fChain->SetBranchAddress("passFilter_EEBadSc", &passFilter_EEBadSc, &b_passFilter_EEBadSc_);
   fChain->SetBranchAddress("passFilter_badChargedHadron", &passFilter_badChargedHadron, &b_passFilter_badChargedHadron_);
   fChain->SetBranchAddress("ak4jet_hf", ak4jet_hf, &b_ak4jet_hf);
   fChain->SetBranchAddress("ak4jet_pf", ak4jet_pf, &b_ak4jet_pf);
   fChain->SetBranchAddress("ak4jet_pt", ak4jet_pt, &b_ak4jet_pt);
   fChain->SetBranchAddress("ak4jet_pt_uncorr", ak4jet_pt_uncorr, &b_ak4jet_pt_uncorr);
   fChain->SetBranchAddress("ak4jet_eta", ak4jet_eta, &b_ak4jet_eta);
   fChain->SetBranchAddress("ak4jet_phi", ak4jet_phi, &b_ak4jet_phi);
   fChain->SetBranchAddress("ak4jet_e", ak4jet_e, &b_ak4jet_e);
   fChain->SetBranchAddress("ak4jet_dr", ak4jet_dr, &b_ak4jet_dr);
   fChain->SetBranchAddress("ak4jet_csv", ak4jet_csv, &b_ak4jet_csv);
   fChain->SetBranchAddress("ak4jet_icsv", ak4jet_icsv, &b_ak4jet_icsv);
   fChain->SetBranchAddress("deltaRAK4AK8", deltaRAK4AK8, &b_deltaRAK4AK8);
   fChain->SetBranchAddress("ak4jet_IDLoose", ak4jet_IDLoose, &b_ak4jet_IDLoose);
   fChain->SetBranchAddress("ak4jet_IDTight", ak4jet_IDTight, &b_ak4jet_IDTight);
   fChain->SetBranchAddress("gen_gra_m", &gen_gra_m, &b_gen_gra_m);
   fChain->SetBranchAddress("gen_gra_pt", &gen_gra_pt, &b_gen_gra_pt);
   fChain->SetBranchAddress("gen_gra_eta", &gen_gra_eta, &b_gen_gra_eta);
   fChain->SetBranchAddress("gen_rad_m", &gen_rad_m, &b_gen_rad_m);
   fChain->SetBranchAddress("gen_rad_pt", &gen_rad_pt, &b_gen_rad_pt);
   fChain->SetBranchAddress("gen_rad_eta", &gen_rad_eta, &b_gen_rad_eta);
   fChain->SetBranchAddress("gen_rad_phi", &gen_rad_phi, &b_gen_rad_phi);
   fChain->SetBranchAddress("gen_tau_eta", &gen_tau_eta, &b_gen_tau_eta);
   fChain->SetBranchAddress("gen_tau_phi", &gen_tau_phi, &b_gen_tau_phi);
   fChain->SetBranchAddress("gen_tau_pt", &gen_tau_pt, &b_gen_tau_pt);
   fChain->SetBranchAddress("gen_tau_e", &gen_tau_e, &b_gen_tau_e);
   fChain->SetBranchAddress("gen_tau_eta_2", &gen_tau_eta_2, &b_gen_tau_eta_2);
   fChain->SetBranchAddress("gen_tau_phi_2", &gen_tau_phi_2, &b_gen_tau_phi_2);
   fChain->SetBranchAddress("gen_tau_pt_2", &gen_tau_pt_2, &b_gen_tau_pt_2);
   fChain->SetBranchAddress("gen_tau_e_2", &gen_tau_e_2, &b_gen_tau_e_2);
   fChain->SetBranchAddress("gen_tau_eta_3", &gen_tau_eta_3, &b_gen_tau_eta_3);
   fChain->SetBranchAddress("gen_tau_phi_3", &gen_tau_phi_3, &b_gen_tau_phi_3);
   fChain->SetBranchAddress("gen_tau_pt_3", &gen_tau_pt_3, &b_gen_tau_pt_3);
   fChain->SetBranchAddress("gen_tau_e_3", &gen_tau_e_3, &b_gen_tau_e_3);

   fChain->SetBranchAddress("gentop_pt", &gentop_pt, &b_gentop_pt);
   fChain->SetBranchAddress("gentop_eta", &gentop_eta, &b_gentop_eta);
   fChain->SetBranchAddress("gentop_phi", &gentop_phi, &b_gentop_phi);
   fChain->SetBranchAddress("gentop_mass", &gentop_mass, &b_gentop_mass);
   fChain->SetBranchAddress("genantitop_pt", &genantitop_pt, &b_genantitop_pt);
   fChain->SetBranchAddress("genantitop_eta", &genantitop_eta, &b_genantitop_eta);
   fChain->SetBranchAddress("genantitop_phi", &genantitop_phi, &b_genantitop_phi);
   fChain->SetBranchAddress("genantitop_mass", &genantitop_mass, &b_genantitop_mass);


   fChain->SetBranchAddress("ptGenVlep", &ptGenVlep, &b_ptGenVlep);
   fChain->SetBranchAddress("etaGenVlep", &etaGenVlep, &b_etaGenVlep);
   fChain->SetBranchAddress("phiGenVlep", &phiGenVlep, &b_phiGenVlep);
   fChain->SetBranchAddress("massGenVlep", &massGenVlep, &b_massGenVlep);
   fChain->SetBranchAddress("ptGenVhad", &ptGenVhad, &b_ptGenVhad);
   fChain->SetBranchAddress("etaGenVhad", &etaGenVhad, &b_etaGenVhad);
   fChain->SetBranchAddress("phiGenVhad", &phiGenVhad, &b_phiGenVhad);
   fChain->SetBranchAddress("massGenVhad", &massGenVhad, &b_massGenVhad);
   fChain->SetBranchAddress("ptGenVlep_2", &ptGenVlep_2, &b_ptGenVlep_2);
   fChain->SetBranchAddress("etaGenVlep_2", &etaGenVlep_2, &b_etaGenVlep_2);
   fChain->SetBranchAddress("phiGenVlep_2", &phiGenVlep_2, &b_phiGenVlep_2);
   fChain->SetBranchAddress("massGenVlep_2", &massGenVlep_2, &b_massGenVlep_2);
   fChain->SetBranchAddress("ptGenV_2", &ptGenV_2, &b_ptGenV_2);
   fChain->SetBranchAddress("etaGenV_2", &etaGenV_2, &b_etaGenV_2);
   fChain->SetBranchAddress("phiGenV_2", &phiGenV_2, &b_phiGenV_2);
   fChain->SetBranchAddress("massGenV_2", &massGenV_2, &b_massGenV_2);
   fChain->SetBranchAddress("ptGenVlep_3", &ptGenVlep_3, &b_ptGenVlep_3);
   fChain->SetBranchAddress("etaGenVlep_3", &etaGenVlep_3, &b_etaGenVlep_3);
   fChain->SetBranchAddress("phiGenVlep_3", &phiGenVlep_3, &b_phiGenVlep_3);
   fChain->SetBranchAddress("massGenVlep_3", &massGenVlep_3, &b_massGenVlep_3);
   fChain->SetBranchAddress("ptGenV_3", &ptGenV_3, &b_ptGenV_3);
   fChain->SetBranchAddress("etaGenV_3", &etaGenV_3, &b_etaGenV_3);
   fChain->SetBranchAddress("phiGenV_3", &phiGenV_3, &b_phiGenV_3);
   fChain->SetBranchAddress("massGenV_3", &massGenV_3, &b_massGenV_3);
fChain->SetBranchAddress("jetAK8puppi_ptJEC_4", &jetAK8puppi_ptJEC_4, &b_jetAK8puppi_ptJEC_4);
fChain->SetBranchAddress("jetAK8puppi_ptJEC_5", &jetAK8puppi_ptJEC_5, &b_jetAK8puppi_ptJEC_5);
fChain->SetBranchAddress("jetAK8puppi_ptJEC_6", &jetAK8puppi_ptJEC_6, &b_jetAK8puppi_ptJEC_6);
fChain->SetBranchAddress("jetAK8puppi_ptJEC_7", &jetAK8puppi_ptJEC_7, &b_jetAK8puppi_ptJEC_7);
fChain->SetBranchAddress("jetAK8puppi_ptJEC_8", &jetAK8puppi_ptJEC_8, &b_jetAK8puppi_ptJEC_8);
fChain->SetBranchAddress("jetAK8puppi_sdJEC_4", &jetAK8puppi_sdJEC_4, &b_jetAK8puppi_sdJEC_4);
fChain->SetBranchAddress("jetAK8puppi_sdJEC_5", &jetAK8puppi_sdJEC_5, &b_jetAK8puppi_sdJEC_5);
fChain->SetBranchAddress("jetAK8puppi_sdJEC_6", &jetAK8puppi_sdJEC_6, &b_jetAK8puppi_sdJEC_6);
fChain->SetBranchAddress("jetAK8puppi_sdJEC_7", &jetAK8puppi_sdJEC_7, &b_jetAK8puppi_sdJEC_7);
fChain->SetBranchAddress("jetAK8puppi_sdJEC_8", &jetAK8puppi_sdJEC_8, &b_jetAK8puppi_sdJEC_8);
   fChain->SetBranchAddress("ptq11", &ptq11, &b_ptq11);
   fChain->SetBranchAddress("etaq11", &etaq11, &b_etaq11);
   fChain->SetBranchAddress("phiq11", &phiq11, &b_phiq11);
   fChain->SetBranchAddress("massq11", &massq11, &b_massq11);
   fChain->SetBranchAddress("ptq12", &ptq12, &b_ptq12);
   fChain->SetBranchAddress("etaq12", &etaq12, &b_etaq12);
   fChain->SetBranchAddress("phiq12", &phiq12, &b_phiq12);
   fChain->SetBranchAddress("massq12", &massq12, &b_massq12);
   fChain->SetBranchAddress("ptq21", &ptq21, &b_ptq21);
   fChain->SetBranchAddress("etaq21", &etaq21, &b_etaq21);
   fChain->SetBranchAddress("phiq21", &phiq21, &b_phiq21);
   fChain->SetBranchAddress("massq21", &massq21, &b_massq21);
   fChain->SetBranchAddress("ptq22", &ptq22, &b_ptq22);
   fChain->SetBranchAddress("etaq22", &etaq22, &b_etaq22);
   fChain->SetBranchAddress("phiq22", &phiq22, &b_phiq22);
   fChain->SetBranchAddress("massq22", &massq22, &b_massq22);   
   fChain->SetBranchAddress("ptq31", &ptq31, &b_ptq31);
   fChain->SetBranchAddress("etaq31", &etaq31, &b_etaq31);
   fChain->SetBranchAddress("phiq31", &phiq31, &b_phiq31);
   fChain->SetBranchAddress("massq31", &massq31, &b_massq31);
   fChain->SetBranchAddress("ptq32", &ptq32, &b_ptq32);
   fChain->SetBranchAddress("etaq32", &etaq32, &b_etaq32);
   fChain->SetBranchAddress("phiq32", &phiq32, &b_phiq32);
   fChain->SetBranchAddress("massq32", &massq32, &b_massq32);
Notify();
}

Bool_t EDBR2PKUTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EDBR2PKUTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

void EDBR2PKUTree::endJob() {
   fout->cd();
   ExTree->Write();
   fout->Write();
   fout->Close();
}

Int_t EDBR2PKUTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EDBR2PKUTree_cxx
