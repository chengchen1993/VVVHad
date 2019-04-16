#define EDBR2PKUTree_cxx
#include "EDBR2PKUTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#define Pi 3.141593
#include "TLorentzVector.h"
//#include "BTagCalibrationStandalone.h"
//
vector<Double_t> generate_weights(TH1* data_npu_estimated, Int_t isForSynch){
	// see SimGeneral/MixingModule/python/mix_2015_25ns_Startup_PoissonOOTPU_cfi.pyy; copy and paste from there:
	const Double_t npu_probs[75] = {
                1.78653e-05 ,
                2.56602e-05 ,
                5.27857e-05 ,
                8.88954e-05 ,
                0.000109362 ,
                0.000140973 ,
                0.000240998 ,
                0.00071209 ,
                0.00130121 ,
                0.00245255 , 
                0.00502589 ,
                0.00919534 ,
                0.0146697 ,
                0.0204126 , 
                0.0267586 ,
                0.0337697 ,
                0.0401478 ,
                0.0450159 ,
                0.0490577 ,
                0.0524855 ,
                0.0548159 ,
                0.0559937 ,
                0.0554468 ,
                0.0537687 ,
                0.0512055 ,
                0.0476713 ,
                0.0435312 ,
                0.0393107 ,
                0.0349812 ,
                0.0307413 ,
                0.0272425 ,
                0.0237115 ,
                0.0208329 ,
                0.0182459 ,
                0.0160712 ,
                0.0142498 ,
                0.012804 ,
                0.011571 ,
                0.010547 ,
                0.00959489 ,
                0.00891718 ,
                0.00829292 ,
                0.0076195 ,
                0.0069806 ,
                0.0062025 ,
                0.00546581 ,
                0.00484127 ,
                0.00407168 ,
                0.00337681 ,
                0.00269893 ,
                0.00212473 ,
                0.00160208 , 
                0.00117884 ,
                0.000859662 ,
                0.000569085 ,
                0.000365431 ,
                0.000243565 ,
                0.00015688 ,
                9.88128e-05 ,
                6.53783e-05 ,
                3.73924e-05 ,
                2.61382e-05 ,
                2.0307e-05 ,
                1.73032e-05 ,
                1.435e-05 ,
                1.36486e-05 ,
                1.35555e-05 ,
                1.37491e-05 ,
                1.34255e-05 ,
                1.33987e-05 ,
                1.34061e-05 ,
                1.34211e-05 ,
                1.34177e-05 ,
                1.32959e-05 ,
                1.33287e-05 };
	if (isForSynch==0) { //OFFICIAL RECIPE
		vector<Double_t> result(75);
		Double_t s = 0.0;
		for(Int_t npu=0; npu<75; ++npu){
			Double_t npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));              
			result[npu] = npu_estimated / npu_probs[npu];
			s += npu_estimated;
		}
		// normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
		for(Int_t npu=0; npu<75; ++npu){
			result[npu] /= s;
		}
		return result;
	}
	else { //THIS IS FOR THE SYNCH ONLY. THIS IS NOT THE OFFICIAL RECIPE!
		vector<Double_t> result(60);
		for(Int_t npu=0; npu<60; ++npu){
			if (data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu))==NULL)
			  result[npu] = 0.;
			else {
				Double_t npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));            
				result[npu] = npu_estimated;
			}
		}
		return result;
	}

}


Double_t bsv (Int_t cud, Double_t x ) // cud=1,2,3 for central,up,down; x for pt
{
  double result=1.0;

  if (cud==1) {  //central
   result=0.498094*((1.+(0.422991*x))/(1.+(0.210944*x)));
  }
  else if (cud==2) { //up
   if(x<20) {result=1;}
   else if(x<30)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.039442747831344604;}
   else if(x<50)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.012669667601585388;}
   else if(x<70) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.011243613436818123;}
   else if(x<100) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.010686126537621021;}
   else if(x<140) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.010994619689881802;}
   else if(x<200) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.011888998560607433;}
   else if(x<300) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.015397069044411182;}
   else if(x<600) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.020918292924761772;}
   else if(x<1000) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.031905386596918106;}
   else {result=1;}
  }
 else if (cud==3) {//down
   if(x<20) {result=1;}
   else if(x<30)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.039442747831344604;}
   else if(x<50)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.012669667601585388;}
   else if(x<70) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.011243613436818123;}
   else if(x<100) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.010686126537621021;}
   else if(x<140) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.010994619689881802;}
   else if(x<200) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.011888998560607433;}
   else if(x<300) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.015397069044411182;}
   else if(x<600) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.020918292924761772;}
   else if(x<1000) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.031905386596918106;}
   else {result=1;}
  }
  return result;
}


Double_t csv (Int_t cud, Double_t x ) // c-jet; cud=1,2,3 for central,up,down; x for pt
{
  double result=1.0;

  if (cud==1) {  //central
   result=0.498094*((1.+(0.422991*x))/(1.+(0.210944*x)));
  }
  else if (cud==2) { //up 
   if(x<20) {result=1;}
   else if(x<30)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.078885495662689209;}
   else if(x<50)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.025339335203170776;}
   else if(x<70) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.022487226873636246;}
   else if(x<100) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.021372253075242043;}
   else if(x<140) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.021989239379763603;}
   else if(x<200) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.023777997121214867;}
   else if(x<300) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.030794138088822365;}
   else if(x<600) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.041836585849523544;}
   else if(x<1000) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.063810773193836212;}
   else {result=1;}
  } 
  else if (cud==3) {//down
   if(x<20) {result=1;}
   else if(x<30)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.078885495662689209;}
   else if(x<50)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.025339335203170776;}
   else if(x<70) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.022487226873636246;}
   else if(x<100) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.021372253075242043;}
   else if(x<140) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.021989239379763603;}
   else if(x<200) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.023777997121214867;}
   else if(x<300) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.030794138088822365;}
   else if(x<600) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.041836585849523544;}
   else if(x<1000) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.063810773193836212;}
   else {result=1;}
  }
  return result;
}
  


Double_t lsv (Int_t cud, Double_t x ) // light flavor; cud=1,2,3 for central,up,down; x for pt
{
  double result=1.0;

  if (cud==1) {  //central
   result=1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x;
  }
  else if (cud==2) { //up 
   result=(1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x)*(1+(0.100485+3.95509e-05*x+-4.90326e-08*x*x));
  }
  else if (cud==3) {//down
   result=(1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x)*(1-(0.100485+3.95509e-05*x+-4.90326e-08*x*x));
  }
  return result;
}

void EDBR2PKUTree::Loop(TString channelname, Double_t XS, TTree *treew, Int_t IsData) {

	std::vector<Double_t> weights_pu1; //these are made with our recipe
	std::vector<Double_t> weights_pu2; //these are made with the official recipe
	//TFile* pileupFile1 = TFile::Open("rereco_2016full_69p2mb_Summer16.root");//pileupDataRun2016BH_63mb_80X.root");  
	TFile* pileupFile1 = TFile::Open("MyDataPileupHistogram.root");
	TH1F* pileupHisto1 = (TH1F*)pileupFile1->Get("pileup");  
	weights_pu1 = generate_weights(pileupHisto1,0);
	pileupFile1->Close();

	//  TFile* pileupFile2 = TFile::Open("puweights.root");  
	TFile* pileupFile2 = TFile::Open("PUxSynch.root");  
	TH1F *pileupHisto2 = (TH1F*)pileupFile2->Get("puweights");
	weights_pu2 = generate_weights(pileupHisto2,1);
	pileupFile2->Close();

	//TFile * input1 = new TFile ("puweights.root");
	//TH1F* hR1= (TH1F*)input1->Get("puweights");
	//zixu
	TFile * input1 = new TFile ("puweight.root");	
	TH1F* hR1= (TH1F*)input1->Get("h2");
	//TFile * input1 = new TFile ("test_mu.root");
	//TH1F* hR1= (TH1F*)input1->Get("hRatio"); //"pileup");//hRatio");


	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();

	Double_t n_deltaRlepjet = 0; 
	Double_t n_delPhijetlep = 0; 
	Double_t ntau = 0;
	Double_t number_qq = 0; 
	Double_t nmassVhad = 0; 
	Double_t nptVlepJEC = 0;
	Double_t nID_e = 0;
	Double_t npt_e = 0;
	Double_t nmet_e = 0; 
	Double_t nnum_bJet_e = 0; 
	Double_t n_delPhijetmet = 0; 

	Double_t nID_mu = 0;
	Double_t npt_mu = 0;
	Double_t nmet_mu = 0; 
	Double_t nnum_bJet_mu = 0; 
	//Double_t nbtb_mu = 0; 

	Double_t nptVhad = 0;
	Double_t yields = 0;
	//TLorentzVector jetV, genjetV;
	//some constants inside this analysis
	Double_t pi_2=1.57079632679;
	Long64_t npp = treew->GetEntries("theWeight>0.");
	Long64_t nmm = treew->GetEntries("theWeight<0.");
//	cout<<"npp="<<npp<<" nmm="<<nmm<<" totaleventnumber="<<totaleventnumber<<endl;

	Double_t nn;
	Double_t eff_and_pu_Weight;
	Double_t eff_and_pu_Weight1;
	Float_t Identical_lumiWeight = XS;//All the events inside a sample are same lumiweight
	//Float_t Identical_lumiWeight = XS/totaleventnumber;//All the events inside a sample are same lumiweight

	Long64_t nbytes = 0, nb = 0;
	//for (Long64_t jentry=0; jentry<10;jentry++)
	for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
                   if (jentry%4000000==0)
                      {std::cout<<jentry<<std::endl;}
		nb = fChain->GetEntry(jentry); 
		nbytes += nb;

                if(event<0){event=event+pow(2,32);}
               Double_t gencorrect=1.0;
                Double_t recocorrect_0eta1p3=1.0;
                Double_t recocorrect_1p3eta2p5=1.0;
                gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC*0.08,-1.2);
                recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC+3.449e-07*pow(jetAK8puppi_ptJEC,2)-2.681e-10*pow(jetAK8puppi_ptJEC,3)+8.674e-14*pow(jetAK8puppi_ptJEC,4)-1.001e-17*pow(jetAK8puppi_ptJEC,5);
                recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC+8.37e-07*pow(jetAK8puppi_ptJEC,2)-5.204e-10*pow(jetAK8puppi_ptJEC,3)+1.454e-13*pow(jetAK8puppi_ptJEC,4)-1.504e-17*pow(jetAK8puppi_ptJEC,5);
                if (fabs(jetAK8puppi_eta)<=1.3){jetAK8puppi_sdJEC=jetAK8puppi_sd*gencorrect*recocorrect_0eta1p3;}
                else if (fabs(jetAK8puppi_eta)<2.5 && fabs(jetAK8puppi_eta)>1.3){jetAK8puppi_sdJEC=jetAK8puppi_sd*gencorrect*recocorrect_1p3eta2p5;}


Double_t gencorrect_2=1.0;
                Double_t recocorrect_0eta1p3_2=1.0;
                Double_t recocorrect_1p3eta2p5_2=1.0;
gencorrect_2=1.006-1.062*pow(jetAK8puppi_ptJEC_2*0.08,-1.2);
                recocorrect_0eta1p3_2=1.093-1.501e-04*jetAK8puppi_ptJEC_2+3.449e-07*pow(jetAK8puppi_ptJEC_2,2)-2.681e-10*pow(jetAK8puppi_ptJEC_2,3)+8.674e-14*pow(jetAK8puppi_ptJEC_2,4)-1.001e-17*pow(jetAK8puppi_ptJEC_2,5);
                recocorrect_1p3eta2p5_2=1.272-5.72e-04*jetAK8puppi_ptJEC_2+8.37e-07*pow(jetAK8puppi_ptJEC_2,2)-5.204e-10*pow(jetAK8puppi_ptJEC_2,3)+1.454e-13*pow(jetAK8puppi_ptJEC_2,4)-1.504e-17*pow(jetAK8puppi_ptJEC_2,5);
                if (fabs(jetAK8puppi_eta_2)<=1.3){jetAK8puppi_sdJEC_2=jetAK8puppi_sd_2*gencorrect_2*recocorrect_0eta1p3_2;}
                else if (fabs(jetAK8puppi_eta_2)<2.5 && fabs(jetAK8puppi_eta_2)>1.3){jetAK8puppi_sdJEC_2=jetAK8puppi_sd_2*gencorrect_2*recocorrect_1p3eta2p5_2;}
      

Double_t gencorrect_3=1.0; 
                Double_t recocorrect_0eta1p3_3=1.0;
                Double_t recocorrect_1p3eta2p5_3=1.0; 
gencorrect_3=1.006-1.062*pow(jetAK8puppi_ptJEC_3*0.08,-1.2);
                recocorrect_0eta1p3_3=1.093-1.501e-04*jetAK8puppi_ptJEC_3+3.449e-07*pow(jetAK8puppi_ptJEC_3,2)-2.681e-10*pow(jetAK8puppi_ptJEC_3,3)+8.674e-14*pow(jetAK8puppi_ptJEC_3,4)-1.001e-17*pow(jetAK8puppi_ptJEC_3,5);
                recocorrect_1p3eta2p5_3=1.272-5.72e-04*jetAK8puppi_ptJEC_3+8.37e-07*pow(jetAK8puppi_ptJEC_3,2)-5.204e-10*pow(jetAK8puppi_ptJEC_3,3)+1.454e-13*pow(jetAK8puppi_ptJEC_3,4)-1.504e-17*pow(jetAK8puppi_ptJEC_3,5);
                if (fabs(jetAK8puppi_eta_3)<=1.3){jetAK8puppi_sdJEC_3=jetAK8puppi_sd_3*gencorrect_3*recocorrect_0eta1p3_3;}
                else if (fabs(jetAK8puppi_eta_3)<2.5 && fabs(jetAK8puppi_eta_3)>1.3){jetAK8puppi_sdJEC_3=jetAK8puppi_sd_3*gencorrect_3*recocorrect_1p3eta2p5_3;}
//jetAK8puppi_sdJEC=corred by fitting.

		pfMET             = Float_t(met);
		pfMETPhi          = Float_t(metPhi);
		l_pt              = Float_t(ptlep1);
		l_eta             = Float_t(etalep1);
		l_phi             = Float_t(philep1);
		ptVhad            = Float_t(ptVhad);
		jet_eta           = Float_t(yVhad);
		jet_phi           = Float_t(phiVhad);
		jet_mass_pruned   = Float_t(massVhadJEC);
                Mj    = Float_t(jetAK8puppi_sdJEC);
                Mj_un = Float_t(jetAK8puppi_sd);
                tau21    = Float_t(jetAK8puppi_tau21);
tau31    = Float_t((jetAK8puppi_tau3>0)*(jetAK8puppi_tau1>0)*jetAK8puppi_tau3/jetAK8puppi_tau1-99*(jetAK8puppi_tau3<0||jetAK8puppi_tau1<0));
tau41    = Float_t((jetAK8puppi_tau4>0)*(jetAK8puppi_tau1>0)*jetAK8puppi_tau4/jetAK8puppi_tau1-99*(jetAK8puppi_tau4<0||jetAK8puppi_tau1<0));
tau32    = Float_t((jetAK8puppi_tau3>0)*(jetAK8puppi_tau2>0)*jetAK8puppi_tau3/jetAK8puppi_tau2-99*(jetAK8puppi_tau3<0||jetAK8puppi_tau2<0));
tau43    = Float_t((jetAK8puppi_tau4>0)*(jetAK8puppi_tau3>0)*jetAK8puppi_tau4/jetAK8puppi_tau3-99*(jetAK8puppi_tau4<0||jetAK8puppi_tau3<0));
jet_tau1_puppi    = Float_t(jetAK8puppi_tau1);
jet_tau2_puppi    = Float_t(jetAK8puppi_tau2);
jet_tau3_puppi    = Float_t(jetAK8puppi_tau3);
jet_tau4_puppi    = Float_t(jetAK8puppi_tau4);
                PTj      = Float_t(jetAK8puppi_ptJEC);
		jetAK8_mass       = Float_t(jetAK8_mass);
//		jet_mass_softdrop = Float_t(sdropJEC);
//		jet_tau2tau1      = Float_t(tau21);
		W_pt              = Float_t(ptVlepJEC);
		W_eta             = Float_t(yVlep);
		W_phi             = Float_t(phiVlep);
                MJJJ              = Float_t(candMasspuppiJEC);

Mj_2    = Float_t(jetAK8puppi_sdJEC_2);
                Mj_un_2 = Float_t(jetAK8puppi_sd_2);
                tau21_2    = Float_t(jetAK8puppi_tau21_2);
PTj_3      = Float_t(jetAK8puppi_ptJEC_3);
                tau42    = Float_t(jetAK8puppi_tau42);
                tau42_3    = Float_t(jetAK8puppi_tau42_3);

tau31_2    = Float_t((jetAK8puppi_tau3_2>0)*(jetAK8puppi_tau1_2>0)*jetAK8puppi_tau3_2/jetAK8puppi_tau1_2-99*(jetAK8puppi_tau3_2<0||jetAK8puppi_tau1_2<0));
tau41_2    = Float_t((jetAK8puppi_tau4_2>0)*(jetAK8puppi_tau1_2>0)*jetAK8puppi_tau4_2/jetAK8puppi_tau1_2-99*(jetAK8puppi_tau4_2<0||jetAK8puppi_tau1_2<0));
tau32_2    = Float_t((jetAK8puppi_tau3_2>0)*(jetAK8puppi_tau2_2>0)*jetAK8puppi_tau3_2/jetAK8puppi_tau2_2-99*(jetAK8puppi_tau3_2<0||jetAK8puppi_tau2_2<0));
tau43_2    = Float_t((jetAK8puppi_tau4_2>0)*(jetAK8puppi_tau3_2>0)*jetAK8puppi_tau4_2/jetAK8puppi_tau3_2-99*(jetAK8puppi_tau4_2<0||jetAK8puppi_tau3_2<0));

jet_tau1_puppi_2    = Float_t(jetAK8puppi_tau1_2);
jet_tau2_puppi_2    = Float_t(jetAK8puppi_tau2_2);
jet_tau3_puppi_2    = Float_t(jetAK8puppi_tau3_2);
jet_tau4_puppi_2    = Float_t(jetAK8puppi_tau4_2);
                PTj_2      = Float_t(jetAK8puppi_ptJEC_2);
                tau42    = Float_t(jetAK8puppi_tau42);
                tau42_2    = Float_t(jetAK8puppi_tau42_2);
Mj_3    = Float_t(jetAK8puppi_sdJEC_3);
Mj_4    = Float_t(jetAK8puppi_sdJEC_4);
Etaj_4    = Float_t(jetAK8puppi_eta_4);
Phij_4    = Float_t(jetAK8puppi_phi_4);

if(jetAK8puppi_ptJEC_2<0)
{
Mj_min=-99;
Mj_mid=-99;
Mj_max=-99;
tau21_max=-99;
tau31_max=-99;
tau41_max=-99;
tau32_max=-99;
tau42_max=-99;
tau43_max=-99;
tau21_mid=-99;
tau31_mid=-99;
tau41_mid=-99;
tau32_mid=-99;
tau42_mid=-99;
tau43_mid=-99;
tau21_min=-99; 
tau31_min=-99;
tau41_min=-99;
tau32_min=-99;
tau42_min=-99;
tau43_min=-99;
PTj_max=-99;
PTj_mid=-99;
PTj_min=-99;
Etaj_max=-99;
Etaj_mid=-99;
Etaj_min=-99;
Phij_max=-99;
Phij_mid=-99;
Phij_min=-99;
DPhi_max_mid=-99;
DPhi_mid_min=-99;
DPhi_min_max=-99;
DPhi_max=-99;
DPhi_min=-99;
DEta_max_mid=-99;
DEta_mid_min=-99;
DEta_min_max=-99;
DEta_max=-99;
DEta_min=-99;
DEta_12=-99;
DEta_13=-99;
DEta_23=-99;
DR_max_mid=-99;
DR_mid_min=-99;
DR_min_max=-99;
DR_max=-99;
DR_min=-99;
}

if(jetAK8puppi_ptJEC_3<0&&jetAK8puppi_ptJEC_2>0)
{
 if(jetAK8puppi_sdJEC_2>jetAK8puppi_sdJEC)
 {
Mj_min=Float_t(jetAK8puppi_sdJEC);
Mj_mid=-99;
Mj_max=Float_t(jetAK8puppi_sdJEC_2);
tau21_max=Float_t(jetAK8puppi_tau2_2/jetAK8puppi_tau1_2);
tau31_max=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau1_2);
tau41_max=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau1_2);
tau32_max=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau2_2);
tau42_max=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau2_2);
tau43_max=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau3_2);
tau21_mid=-99;
tau31_mid=-99;
tau41_mid=-99;
tau32_mid=-99;
tau42_mid=-99;
tau43_mid=-99;
tau21_min=Float_t(jetAK8puppi_tau2/jetAK8puppi_tau1);
tau31_min=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau1);
tau41_min=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau1);
tau32_min=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau2);
tau42_min=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau2);
tau43_min=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau3);
PTj_min=Float_t(jetAK8puppi_ptJEC);
PTj_mid=-99;
PTj_max=Float_t(jetAK8puppi_ptJEC_2);
Etaj_min=Float_t(jetAK8puppi_eta);
Etaj_mid=-99;
Etaj_max=Float_t(jetAK8puppi_eta_2);
Phij_min=Float_t(jetAK8puppi_phi);
Phij_mid=-99;
Phij_max=Float_t(jetAK8puppi_phi_2);
DPhi_max_mid=-99;
DPhi_mid_min=-99;
DPhi_min_max=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_2)-3.14159));
DEta_min_max=fabs(jetAK8puppi_eta-jetAK8puppi_eta_2);
DR_max_mid=-99;
DR_mid_min=-99;
DR_min_max=Float_t(sqrt(DPhi_min_max*DPhi_min_max+(jetAK8puppi_eta-jetAK8puppi_eta_2)*(jetAK8puppi_eta-jetAK8puppi_eta_2)));

 }
 if(jetAK8puppi_sdJEC_2<jetAK8puppi_sdJEC)
 {
Mj_min=Float_t(jetAK8puppi_sdJEC_2);
Mj_mid=-99;
Mj_max=Float_t(jetAK8puppi_sdJEC);
tau21_max=Float_t(jetAK8puppi_tau2/jetAK8puppi_tau1);
tau31_max=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau1);
tau41_max=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau1);
tau32_max=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau2);
tau42_max=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau2);
tau43_max=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau3);
tau21_mid=-99;
tau31_mid=-99;
tau41_mid=-99;
tau32_mid=-99;
tau42_mid=-99;
tau43_mid=-99;
tau21_min=Float_t(jetAK8puppi_tau2_2/jetAK8puppi_tau1_2);
tau31_min=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau1_2);
tau41_min=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau1_2);
tau32_min=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau2_2);
tau42_min=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau2_2);
tau43_min=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau3_2);
PTj_min=Float_t(jetAK8puppi_ptJEC_2);
PTj_mid=-99;
PTj_max=Float_t(jetAK8puppi_ptJEC);
Etaj_min=Float_t(jetAK8puppi_eta_2);
Etaj_mid=-99;
Etaj_max=Float_t(jetAK8puppi_eta);
Phij_min=Float_t(jetAK8puppi_phi_2);
Phij_mid=-99;
Phij_max=Float_t(jetAK8puppi_phi);
DPhi_max_mid=-99;
DPhi_mid_min=-99;
DPhi_min_max=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_2)-3.14159));
DEta_min_max=fabs(jetAK8puppi_eta-jetAK8puppi_eta_2);
DR_max_mid=-99;
DR_mid_min=-99;
DR_min_max=Float_t(sqrt(DPhi_min_max*DPhi_min_max+(jetAK8puppi_eta-jetAK8puppi_eta_2)*(jetAK8puppi_eta-jetAK8puppi_eta_2)));
 } 
}

if(jetAK8puppi_ptJEC_3>0)
{
 if(jetAK8puppi_sdJEC_3>jetAK8puppi_sdJEC_2&&jetAK8puppi_sdJEC_2>jetAK8puppi_sdJEC)
 {
Mj_min=Float_t(jetAK8puppi_sdJEC);
Mj_mid=Float_t(jetAK8puppi_sdJEC_2);
Mj_max=Float_t(jetAK8puppi_sdJEC_3);
tau21_max=Float_t(jetAK8puppi_tau2_3/jetAK8puppi_tau1_3);
tau31_max=Float_t(jetAK8puppi_tau3_3/jetAK8puppi_tau1_3);
tau41_max=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau1_3);
tau32_max=Float_t(jetAK8puppi_tau3_3/jetAK8puppi_tau2_3);
tau42_max=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau2_3);
tau43_max=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau3_3);
tau21_mid=Float_t(jetAK8puppi_tau2_2/jetAK8puppi_tau1_2);
tau31_mid=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau1_2);
tau41_mid=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau1_2);
tau32_mid=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau2_2);
tau42_mid=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau2_2);
tau43_mid=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau3_2);
tau21_min=Float_t(jetAK8puppi_tau2/jetAK8puppi_tau1);
tau31_min=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau1);
tau41_min=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau1);
tau32_min=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau2);
tau42_min=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau2);
tau43_min=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau3);
PTj_min=Float_t(jetAK8puppi_ptJEC);
PTj_mid=Float_t(jetAK8puppi_ptJEC_2);
PTj_max=Float_t(jetAK8puppi_ptJEC_3);
Etaj_min=Float_t(jetAK8puppi_eta);
Etaj_mid=Float_t(jetAK8puppi_eta_2);
Etaj_max=Float_t(jetAK8puppi_eta_3);
Phij_min=Float_t(jetAK8puppi_phi);
Phij_mid=Float_t(jetAK8puppi_phi_2);
Phij_max=Float_t(jetAK8puppi_phi_3);
DPhi_max_mid=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi_3-jetAK8puppi_phi_2)-3.14159));
DPhi_mid_min=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_2)-3.14159));
DPhi_min_max=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_3)-3.14159));
DEta_max_mid=fabs(jetAK8puppi_eta_3-jetAK8puppi_eta_2);
DEta_mid_min=fabs(jetAK8puppi_eta-jetAK8puppi_eta_2);
DEta_min_max=fabs(jetAK8puppi_eta-jetAK8puppi_eta_3);
DR_max_mid=Float_t(sqrt(DPhi_max_mid*DPhi_max_mid+(jetAK8puppi_eta_3-jetAK8puppi_eta_2)*(jetAK8puppi_eta_3-jetAK8puppi_eta_2)));
DR_mid_min=Float_t(sqrt(DPhi_mid_min*DPhi_mid_min+(jetAK8puppi_eta-jetAK8puppi_eta_2)*(jetAK8puppi_eta-jetAK8puppi_eta_2)));
DR_min_max=Float_t(sqrt(DPhi_min_max*DPhi_min_max+(jetAK8puppi_eta-jetAK8puppi_eta_3)*(jetAK8puppi_eta-jetAK8puppi_eta_3)));
 }           
 if(jetAK8puppi_sdJEC_3>jetAK8puppi_sdJEC&&jetAK8puppi_sdJEC>jetAK8puppi_sdJEC_2)
 {
Mj_min=Float_t(jetAK8puppi_sdJEC_2);
Mj_mid=Float_t(jetAK8puppi_sdJEC);
Mj_max=Float_t(jetAK8puppi_sdJEC_3);
tau21_max=Float_t(jetAK8puppi_tau2_3/jetAK8puppi_tau1_3);
tau31_max=Float_t(jetAK8puppi_tau3_3/jetAK8puppi_tau1_3);
tau41_max=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau1_3);
tau32_max=Float_t(jetAK8puppi_tau3_3/jetAK8puppi_tau2_3);
tau42_max=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau2_3);
tau43_max=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau3_3);
tau21_mid=Float_t(jetAK8puppi_tau2/jetAK8puppi_tau1);
tau31_mid=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau1);
tau41_mid=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau1);
tau32_mid=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau2);
tau42_mid=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau2);
tau43_mid=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau3);
tau21_min=Float_t(jetAK8puppi_tau2_2/jetAK8puppi_tau1_2);
tau31_min=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau1_2);
tau41_min=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau1_2);
tau32_min=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau2_2);
tau42_min=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau2_2);
tau43_min=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau3_2);
PTj_min=Float_t(jetAK8puppi_ptJEC_2);
PTj_mid=Float_t(jetAK8puppi_ptJEC);
PTj_max=Float_t(jetAK8puppi_ptJEC_3);
Etaj_min=Float_t(jetAK8puppi_eta_2);
Etaj_mid=Float_t(jetAK8puppi_eta);
Etaj_max=Float_t(jetAK8puppi_eta_3);
Phij_min=Float_t(jetAK8puppi_phi_2);
Phij_mid=Float_t(jetAK8puppi_phi);
Phij_max=Float_t(jetAK8puppi_phi_3);
DPhi_max_mid=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_3)-3.14159));
DPhi_mid_min=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_2)-3.14159));
DPhi_min_max=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi_3-jetAK8puppi_phi_2)-3.14159));
DEta_max_mid=fabs(jetAK8puppi_eta-jetAK8puppi_eta_3);
DEta_mid_min=fabs(jetAK8puppi_eta-jetAK8puppi_eta_2);
DEta_min_max=fabs(jetAK8puppi_eta_3-jetAK8puppi_eta_2);
DR_max_mid=Float_t(sqrt(DPhi_max_mid*DPhi_max_mid+(jetAK8puppi_eta_3-jetAK8puppi_eta)*(jetAK8puppi_eta_3-jetAK8puppi_eta)));
DR_mid_min=Float_t(sqrt(DPhi_mid_min*DPhi_mid_min+(jetAK8puppi_eta-jetAK8puppi_eta_2)*(jetAK8puppi_eta-jetAK8puppi_eta_2)));
DR_min_max=Float_t(sqrt(DPhi_min_max*DPhi_min_max+(jetAK8puppi_eta_2-jetAK8puppi_eta_3)*(jetAK8puppi_eta_2-jetAK8puppi_eta_3)));
 }
 if(jetAK8puppi_sdJEC_2>jetAK8puppi_sdJEC&&jetAK8puppi_sdJEC>jetAK8puppi_sdJEC_3)
 {
Mj_min=Float_t(jetAK8puppi_sdJEC_3);
Mj_mid=Float_t(jetAK8puppi_sdJEC);
Mj_max=Float_t(jetAK8puppi_sdJEC_2);
tau21_max=Float_t(jetAK8puppi_tau2_2/jetAK8puppi_tau1_2);
tau31_max=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau1_2);
tau41_max=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau1_2);
tau32_max=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau2_2);
tau42_max=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau2_2);
tau43_max=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau3_2);
tau21_mid=Float_t(jetAK8puppi_tau2/jetAK8puppi_tau1);
tau31_mid=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau1);
tau41_mid=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau1);
tau32_mid=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau2);
tau42_mid=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau2);
tau43_mid=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau3);
tau21_min=Float_t(jetAK8puppi_tau2_3/jetAK8puppi_tau1_3);
tau31_min=Float_t(jetAK8puppi_tau3_3/jetAK8puppi_tau1_3);
tau41_min=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau1_3);
tau32_min=Float_t(jetAK8puppi_tau3_3/jetAK8puppi_tau2_3);
tau42_min=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau2_3);
tau43_min=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau3_3);
PTj_min=Float_t(jetAK8puppi_ptJEC_3);
PTj_mid=Float_t(jetAK8puppi_ptJEC);
PTj_max=Float_t(jetAK8puppi_ptJEC_2);
Etaj_min=Float_t(jetAK8puppi_eta_3);
Etaj_mid=Float_t(jetAK8puppi_eta);
Etaj_max=Float_t(jetAK8puppi_eta_2);
Phij_min=Float_t(jetAK8puppi_phi_3);
Phij_mid=Float_t(jetAK8puppi_phi);
Phij_max=Float_t(jetAK8puppi_phi_2);
DPhi_max_mid=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_2)-3.14159));
DPhi_mid_min=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_3)-3.14159));
DPhi_min_max=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi_3-jetAK8puppi_phi_2)-3.14159));
DEta_max_mid=fabs(jetAK8puppi_eta-jetAK8puppi_eta_2);
DEta_mid_min=fabs(jetAK8puppi_eta-jetAK8puppi_eta_3);
DEta_min_max=fabs(jetAK8puppi_eta_3-jetAK8puppi_eta_2);
DR_max_mid=Float_t(sqrt(DPhi_max_mid*DPhi_max_mid+(jetAK8puppi_eta-jetAK8puppi_eta_2)*(jetAK8puppi_eta-jetAK8puppi_eta_2)));
DR_mid_min=Float_t(sqrt(DPhi_mid_min*DPhi_mid_min+(jetAK8puppi_eta-jetAK8puppi_eta_2)*(jetAK8puppi_eta-jetAK8puppi_eta_2)));
DR_min_max=Float_t(sqrt(DPhi_min_max*DPhi_min_max+(jetAK8puppi_eta_2-jetAK8puppi_eta_3)*(jetAK8puppi_eta_2-jetAK8puppi_eta_3)));
 }
 if(jetAK8puppi_sdJEC_2>jetAK8puppi_sdJEC_3&&jetAK8puppi_sdJEC_3>jetAK8puppi_sdJEC)
 {
Mj_min=Float_t(jetAK8puppi_sdJEC);
Mj_mid=Float_t(jetAK8puppi_sdJEC_3);
Mj_max=Float_t(jetAK8puppi_sdJEC_2);
tau21_max=Float_t(jetAK8puppi_tau2_2/jetAK8puppi_tau1_2);
tau31_max=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau1_2);
tau41_max=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau1_2);
tau32_max=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau2_2);
tau42_max=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau2_2);
tau43_max=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau3_2);
tau21_mid=Float_t(jetAK8puppi_tau2_3/jetAK8puppi_tau1_3);
tau31_mid=Float_t(jetAK8puppi_tau3_3/jetAK8puppi_tau1_3);
tau41_mid=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau1_3);
tau32_mid=Float_t(jetAK8puppi_tau3_3/jetAK8puppi_tau2_3);
tau42_mid=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau2_3);
tau43_mid=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau3_3);
tau21_min=Float_t(jetAK8puppi_tau2/jetAK8puppi_tau1);
tau31_min=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau1);
tau41_min=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau1);
tau32_min=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau2);
tau42_min=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau2);
tau43_min=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau3);
PTj_min=Float_t(jetAK8puppi_ptJEC);
PTj_mid=Float_t(jetAK8puppi_ptJEC_3);
PTj_max=Float_t(jetAK8puppi_ptJEC_2);
Etaj_min=Float_t(jetAK8puppi_eta);
Etaj_mid=Float_t(jetAK8puppi_eta_3);
Etaj_max=Float_t(jetAK8puppi_eta_2);
Phij_min=Float_t(jetAK8puppi_phi);
Phij_mid=Float_t(jetAK8puppi_phi_3);
Phij_max=Float_t(jetAK8puppi_phi_2);
DPhi_max_mid=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi_3-jetAK8puppi_phi_2)-3.14159));
DPhi_mid_min=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_3)-3.14159));
DPhi_min_max=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_2)-3.14159));
DEta_max_mid=fabs(jetAK8puppi_eta_3-jetAK8puppi_eta_2);
DEta_mid_min=fabs(jetAK8puppi_eta-jetAK8puppi_eta_3);
DEta_min_max=fabs(jetAK8puppi_eta-jetAK8puppi_eta_2);
DR_max_mid=Float_t(sqrt(DPhi_max_mid*DPhi_max_mid+(jetAK8puppi_eta_3-jetAK8puppi_eta_2)*(jetAK8puppi_eta_3-jetAK8puppi_eta_2)));
DR_mid_min=Float_t(sqrt(DPhi_mid_min*DPhi_mid_min+(jetAK8puppi_eta-jetAK8puppi_eta_3)*(jetAK8puppi_eta-jetAK8puppi_eta_3)));
DR_min_max=Float_t(sqrt(DPhi_min_max*DPhi_min_max+(jetAK8puppi_eta-jetAK8puppi_eta_2)*(jetAK8puppi_eta-jetAK8puppi_eta_2)));
 }
 if(jetAK8puppi_sdJEC>jetAK8puppi_sdJEC_2&&jetAK8puppi_sdJEC_2>jetAK8puppi_sdJEC_3)
 {
Mj_min=Float_t(jetAK8puppi_sdJEC_3);
Mj_mid=Float_t(jetAK8puppi_sdJEC_2);
Mj_max=Float_t(jetAK8puppi_sdJEC);
tau21_max=Float_t(jetAK8puppi_tau2/jetAK8puppi_tau1);
tau31_max=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau1);
tau41_max=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau1);
tau32_max=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau2);
tau42_max=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau2);
tau43_max=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau3);
tau21_mid=Float_t(jetAK8puppi_tau2_2/jetAK8puppi_tau1_2);
tau31_mid=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau1_2);
tau41_mid=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau1_2);
tau32_mid=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau2_2);
tau42_mid=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau2_2);
tau43_mid=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau3_2);
tau21_min=Float_t(jetAK8puppi_tau2_3/jetAK8puppi_tau1_3);
tau31_min=Float_t(jetAK8puppi_tau3_3/jetAK8puppi_tau1_3);
tau41_min=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau1_3);
tau32_min=Float_t(jetAK8puppi_tau3_3/jetAK8puppi_tau2_3);
tau42_min=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau2_3);
tau43_min=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau3_3);
PTj_min=Float_t(jetAK8puppi_ptJEC_3);
PTj_mid=Float_t(jetAK8puppi_ptJEC_2);
PTj_max=Float_t(jetAK8puppi_ptJEC);
Etaj_min=Float_t(jetAK8puppi_eta_3);
Etaj_mid=Float_t(jetAK8puppi_eta_2);
Etaj_max=Float_t(jetAK8puppi_eta);
Phij_min=Float_t(jetAK8puppi_phi_3);
Phij_mid=Float_t(jetAK8puppi_phi_2);
Phij_max=Float_t(jetAK8puppi_phi);
DPhi_max_mid=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_2)-3.14159));
DPhi_mid_min=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi_3-jetAK8puppi_phi_2)-3.14159));
DPhi_min_max=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_3)-3.14159));
DEta_max_mid=fabs(jetAK8puppi_eta-jetAK8puppi_eta_2);
DEta_mid_min=fabs(jetAK8puppi_eta_3-jetAK8puppi_eta_2);
DEta_min_max=fabs(jetAK8puppi_eta-jetAK8puppi_eta_3);
DR_max_mid=Float_t(sqrt(DPhi_max_mid*DPhi_max_mid+(jetAK8puppi_eta-jetAK8puppi_eta_2)*(jetAK8puppi_eta-jetAK8puppi_eta_2)));
DR_mid_min=Float_t(sqrt(DPhi_mid_min*DPhi_mid_min+(jetAK8puppi_eta_3-jetAK8puppi_eta_2)*(jetAK8puppi_eta_3-jetAK8puppi_eta_2)));
DR_min_max=Float_t(sqrt(DPhi_min_max*DPhi_min_max+(jetAK8puppi_eta-jetAK8puppi_eta_3)*(jetAK8puppi_eta-jetAK8puppi_eta_3)));
 }
 if(jetAK8puppi_sdJEC>jetAK8puppi_sdJEC_3&&jetAK8puppi_sdJEC_3>jetAK8puppi_sdJEC_2)
 {
Mj_min=Float_t(jetAK8puppi_sdJEC_2);
Mj_mid=Float_t(jetAK8puppi_sdJEC_3);
Mj_max=Float_t(jetAK8puppi_sdJEC);
tau21_max=Float_t(jetAK8puppi_tau2/jetAK8puppi_tau1);
tau31_max=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau1);
tau41_max=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau1);
tau32_max=Float_t(jetAK8puppi_tau3/jetAK8puppi_tau2);
tau42_max=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau2);
tau43_max=Float_t(jetAK8puppi_tau4/jetAK8puppi_tau3);
tau21_mid=Float_t(jetAK8puppi_tau2_3/jetAK8puppi_tau1_3);
tau31_mid=Float_t(jetAK8puppi_tau3_3/jetAK8puppi_tau1_3);
tau41_mid=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau1_3);
tau32_mid=Float_t(jetAK8puppi_tau3_3/jetAK8puppi_tau2_3);
tau42_mid=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau2_3);
tau43_mid=Float_t(jetAK8puppi_tau4_3/jetAK8puppi_tau3_3);
tau21_min=Float_t(jetAK8puppi_tau2_2/jetAK8puppi_tau1_2);
tau31_min=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau1_2);
tau41_min=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau1_2);
tau32_min=Float_t(jetAK8puppi_tau3_2/jetAK8puppi_tau2_2);
tau42_min=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau2_2);
tau43_min=Float_t(jetAK8puppi_tau4_2/jetAK8puppi_tau3_2);
PTj_min=Float_t(jetAK8puppi_ptJEC_2);
PTj_mid=Float_t(jetAK8puppi_ptJEC_3);
PTj_max=Float_t(jetAK8puppi_ptJEC);
Etaj_min=Float_t(jetAK8puppi_eta_2);
Etaj_mid=Float_t(jetAK8puppi_eta_3);
Etaj_max=Float_t(jetAK8puppi_eta);
Phij_min=Float_t(jetAK8puppi_phi_2);
Phij_mid=Float_t(jetAK8puppi_phi_3);
Phij_max=Float_t(jetAK8puppi_phi);
DPhi_max_mid=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_3)-3.14159));
DPhi_mid_min=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi_3-jetAK8puppi_phi_2)-3.14159));
DPhi_min_max=Float_t(3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_2)-3.14159));
DEta_max_mid=fabs(jetAK8puppi_eta-jetAK8puppi_eta_3);
DEta_mid_min=fabs(jetAK8puppi_eta_3-jetAK8puppi_eta_2);
DEta_min_max=fabs(jetAK8puppi_eta-jetAK8puppi_eta_2);
DR_max_mid=Float_t(sqrt(DPhi_max_mid*DPhi_max_mid+(jetAK8puppi_eta_3-jetAK8puppi_eta)*(jetAK8puppi_eta_3-jetAK8puppi_eta)));
DR_mid_min=Float_t(sqrt(DPhi_mid_min*DPhi_mid_min+(jetAK8puppi_eta_3-jetAK8puppi_eta_2)*(jetAK8puppi_eta_3-jetAK8puppi_eta_2)));
DR_min_max=Float_t(sqrt(DPhi_min_max*DPhi_min_max+(jetAK8puppi_eta-jetAK8puppi_eta_2)*(jetAK8puppi_eta-jetAK8puppi_eta_2)));
 }
}

t21t31t41_max  = Float_t(pow(tau21_max*tau31_max*tau41_max,1.0/3.0));
if(tau21_mid>0)
{
t21t31t41_mid  = Float_t(pow(tau21_mid*tau31_mid*tau41_mid,1.0/3.0));
}
else
{
t21t31t41_mid  = -99;
}
t21t31t41_min  = Float_t(pow(tau21_min*tau31_min*tau41_min,1.0/3.0));
S_t21t31t41_max  = Float_t((tau21_max+tau31_max+tau41_max)/3.0);
S_t21t31t41_mid  = Float_t((tau21_mid+tau31_mid+tau41_mid)/3.0);
S_t21t31t41_min  = Float_t((tau21_min+tau31_min+tau41_min)/3.0);

if(jetAK8puppi_ptJEC>0)
{
Mj_mean=Float_t((jetAK8puppi_sdJEC*(jetAK8puppi_sdJEC>0)+jetAK8puppi_sdJEC_2*(jetAK8puppi_sdJEC_2>0)+jetAK8puppi_sdJEC_3*(jetAK8puppi_sdJEC_3>0))/((jetAK8puppi_sdJEC>0)+(jetAK8puppi_sdJEC_2>0)+(jetAK8puppi_sdJEC_3>0)));
}

if(jetAK8puppi_ptJEC<0)
{
Mj_mean=-99;
}

if(jetAK8puppi_ptJEC_2>0)
{
Pt2dPt1=Float_t(PTj_2/PTj);
}
if(jetAK8puppi_ptJEC_2<0)
{
Pt2dPt1=-99;
}

if(jetAK8puppi_ptJEC_3>0)
{
Pt3dPt1=Float_t(PTj_3/PTj);
}
if(jetAK8puppi_ptJEC_3<0)
{
Pt3dPt1=-99;
}
if(jetAK8puppi_ptJEC>0&&jetAK8puppi_ptJEC_2>0)
{
Phij_12=3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_2)-3.14159);
DEta_12=fabs(jetAK8puppi_eta-jetAK8puppi_eta_2);
DR_12=sqrt(Phij_12*Phij_12+(jetAK8puppi_eta-jetAK8puppi_eta_2)*(jetAK8puppi_eta-jetAK8puppi_eta_2));
}
if(jetAK8puppi_ptJEC<0||jetAK8puppi_ptJEC_2<0)
{
Phij_12=-99;
DR_12=-99;
}

if(jetAK8puppi_ptJEC>0&&jetAK8puppi_ptJEC_3>0)
{
Phij_13=3.14159-fabs(fabs(jetAK8puppi_phi-jetAK8puppi_phi_3)-3.14159);
DEta_13=fabs(jetAK8puppi_eta-jetAK8puppi_eta_3);
DR_13=sqrt(Phij_13*Phij_13+(jetAK8puppi_eta-jetAK8puppi_eta_3)*(jetAK8puppi_eta-jetAK8puppi_eta_3));
}
if(jetAK8puppi_ptJEC<0||jetAK8puppi_ptJEC_3<0)
{
Phij_13=-99;
DR_13=-99;
}

if(jetAK8puppi_ptJEC_3>0&&jetAK8puppi_ptJEC_2>0)
{
Phij_23=3.14159-fabs(fabs(jetAK8puppi_phi_3-jetAK8puppi_phi_2)-3.14159);
DEta_23=fabs(jetAK8puppi_eta_3-jetAK8puppi_eta_2);
DR_23=sqrt(Phij_23*Phij_23+(jetAK8puppi_eta_3-jetAK8puppi_eta_2)*(jetAK8puppi_eta_3-jetAK8puppi_eta_2));
}
if(jetAK8puppi_ptJEC_3<0||jetAK8puppi_ptJEC_2<0)
{
Phij_23=-99;
DR_23=-99;
}

if(jetAK8puppi_ptJEC_2<0)
{
DPhi_max=-99;
DPhi_min=-99;
DR_max=-99;
DR_min=-99;
}

if(jetAK8puppi_ptJEC_2>0&&jetAK8puppi_ptJEC_3<0)
{
DPhi_max=Phij_12;
DPhi_min=Phij_12;
DEta_max=DEta_12;
DEta_min=DEta_12;
DR_max=DR_12;
DR_min=DR_12;
}

if(jetAK8puppi_ptJEC_2>0&&jetAK8puppi_ptJEC_3>0)
{
   if(Phij_12>=Phij_13&&Phij_12>=Phij_23)
   {
     DPhi_max=Phij_12;
   }
   if(Phij_13>=Phij_12&&Phij_13>=Phij_23)
   {
     DPhi_max=Phij_13;
   }
   if(Phij_23>=Phij_12&&Phij_23>=Phij_13)
   {
     DPhi_max=Phij_23;
   }

   if(Phij_12<=Phij_13&&Phij_12<=Phij_23)
   {
     DPhi_min=Phij_12;
   }
   if(Phij_13<=Phij_12&&Phij_13<=Phij_23)
   {
     DPhi_min=Phij_13;
   }
   if(Phij_23<=Phij_12&&Phij_23<=Phij_13)
   {
     DPhi_min=Phij_23;
   }
   

   if(DEta_12>=DEta_13&&DEta_12>=DEta_23)
   {
     DEta_max=DEta_12;
   }
   if(DEta_13>=DEta_12&&DEta_13>=DEta_23)
   {
     DEta_max=DEta_13;
   }
   if(DEta_23>=DEta_12&&DEta_23>=DEta_13)
   {
     DEta_max=DEta_23;
   }

   if(DEta_12<=DEta_13&&DEta_12<=DEta_23)
   {
     DEta_min=DEta_12;
   }
   if(DEta_13<=DEta_12&&DEta_13<=DEta_23)
   {
     DEta_min=DEta_13;
   }
   if(DEta_23<=DEta_12&&DEta_23<=DEta_13)
   {
     DEta_min=DEta_23;
   }


   if(DR_12>=DR_13&&DR_12>=DR_23)
   {
     DR_max=DR_12;
   }
   if(DR_13>=DR_12&&DR_13>=DR_23)
   {
     DR_max=DR_13;
   }
   if(DR_23>=DR_12&&DR_23>=DR_13)
   {
     DR_max=DR_23;
   }

   if(DR_12<=DR_13&&DR_12<=DR_23)
   {
     DR_min=DR_12;
   }
   if(DR_13<=DR_12&&DR_13<=DR_23)
   {
     DR_min=DR_13;
   }
   if(DR_23<=DR_12&&DR_23<=DR_13)
   {
     DR_min=DR_23;
   }
}

DPhi_min2=DPhi_min;
DPhi_min3=DPhi_min;
DPhi_max2=DPhi_max;
DPhi_max3=DPhi_max;

if(jetAK8puppi_ptJEC_2>0&&jetAK8puppi_ptJEC_3>0)
{
PTj_23 =sqrt(PTj_2*PTj_2+PTj_3*PTj_3+2.0*PTj_2*PTj_3*cos(Phij_23));
}
if(jetAK8puppi_ptJEC_2<0||jetAK8puppi_ptJEC_3<0)
{
PTj_23 = -99;
}

HT=Float_t((jetAK8puppi_ptJEC>0)*jetAK8puppi_ptJEC+(jetAK8puppi_ptJEC_2>0)*jetAK8puppi_ptJEC_2+(jetAK8puppi_ptJEC_3>0)*jetAK8puppi_ptJEC_3+(jetAK8puppi_ptJEC_4>0)*jetAK8puppi_ptJEC_4+(jetAK8puppi_ptJEC_5>0)*jetAK8puppi_ptJEC_5+(jetAK8puppi_ptJEC_6>0)*jetAK8puppi_ptJEC_6+(jetAK8puppi_ptJEC_7>0)*jetAK8puppi_ptJEC_7+(jetAK8puppi_ptJEC_8>0)*jetAK8puppi_ptJEC_8);

ST=Float_t(HT+MET_et);

Nj4=(ak4jet_pt[0]>0)+(ak4jet_pt[1]>0)+(ak4jet_pt[2]>0)+(ak4jet_pt[3]>0)+(ak4jet_pt[4]>0)+(ak4jet_pt[5]>0)+(ak4jet_pt[6]>0)+(ak4jet_pt[7]>0);

Nj8=(jetAK8puppi_ptJEC>0)+(jetAK8puppi_ptJEC_2>0)+(jetAK8puppi_ptJEC_3>0)+(jetAK8puppi_ptJEC_4>0)+(jetAK8puppi_ptJEC_5>0)+(jetAK8puppi_ptJEC_6>0)+(jetAK8puppi_ptJEC_7>0)+(jetAK8puppi_ptJEC_8>0);

Mj_un_3 = Float_t(jetAK8puppi_sd_3);
tau21_3    = Float_t(jetAK8puppi_tau21_3);
tau31_3    = Float_t((jetAK8puppi_tau3_3>0)*(jetAK8puppi_tau1_3>0)*jetAK8puppi_tau3_3/jetAK8puppi_tau1_3-99*(jetAK8puppi_tau3_3<0||jetAK8puppi_tau1_3<0));
tau41_3    = Float_t((jetAK8puppi_tau4_3>0)*(jetAK8puppi_tau1_3>0)*jetAK8puppi_tau4_3/jetAK8puppi_tau1_3-99*(jetAK8puppi_tau4_3<0||jetAK8puppi_tau1_3<0));
tau32_3    = Float_t((jetAK8puppi_tau3_3>0)*(jetAK8puppi_tau2_3>0)*jetAK8puppi_tau3_3/jetAK8puppi_tau2_3-99*(jetAK8puppi_tau3_3<0||jetAK8puppi_tau2_3<0));
tau43_3    = Float_t((jetAK8puppi_tau4_3>0)*(jetAK8puppi_tau3_3>0)*jetAK8puppi_tau4_3/jetAK8puppi_tau3_3-99*(jetAK8puppi_tau4_3<0||jetAK8puppi_tau3_3<0));
jet_tau1_puppi_3    = Float_t(jetAK8puppi_tau1_3);
jet_tau2_puppi_3    = Float_t(jetAK8puppi_tau2_3);
jet_tau3_puppi_3    = Float_t(jetAK8puppi_tau3_3);
jet_tau4_puppi_3    = Float_t(jetAK8puppi_tau4_3);
t21t31t41  = Float_t(pow(tau41*tau31*tau21,1.0/3.0));
t21t31t41_2  = Float_t(pow(tau41_2*tau31_2*tau21_2,1.0/3.0));
t21t31t41_3  = Float_t(pow(tau41_3*tau31_3*tau21_3,1.0/3.0));

if(jetAK8puppi_ptJEC>0&&jetAK8puppi_ptJEC_2>0&&jetAK8puppi_ptJEC_3>0)
{
pt1pt2pt3  = Float_t(pow(PTj*PTj_2*PTj_3,1.0/3.0));
}
if(jetAK8puppi_ptJEC<0||jetAK8puppi_ptJEC_2<0||jetAK8puppi_ptJEC_3<0)
{
pt1pt2pt3  = -99.;
}

Etaj=Float_t(jetAK8puppi_eta);
Etaj_2=Float_t(jetAK8puppi_eta_2);
Etaj_3=Float_t(jetAK8puppi_eta_3);

Phij=Float_t(jetAK8puppi_phi);
Phij_2=Float_t(jetAK8puppi_phi_2);
Phij_3=Float_t(jetAK8puppi_phi_3);

PTj_4=Float_t(jetAK8puppi_ptJEC_4);
PTj_5=Float_t(jetAK8puppi_ptJEC_5);
PTj_6=Float_t(jetAK8puppi_ptJEC_6);
PTj_7=Float_t(jetAK8puppi_ptJEC_7);
PTj_8=Float_t(jetAK8puppi_ptJEC_8);
newgen_gra_m=Float_t(gen_gra_m);
newgen_gra_pt=Float_t(gen_gra_pt);
newgen_gra_eta=Float_t(gen_gra_eta);
newgen_rad_m=Float_t(gen_rad_m);
newgen_rad_pt=Float_t(gen_rad_pt);
newgen_rad_eta=Float_t(gen_rad_eta);
newgen_tau_e=Float_t(gen_tau_e);
newgen_tau_pt=Float_t(gen_tau_pt);
newgen_tau_eta=Float_t(gen_tau_eta);
newgen_tau_phi=Float_t(gen_tau_phi);
newgen_tau_e_2=Float_t(gen_tau_e_2);
newgen_tau_pt_2=Float_t(gen_tau_pt_2);
newgen_tau_eta_2=Float_t(gen_tau_eta_2);
newgen_tau_phi_2=Float_t(gen_tau_phi_2);
newgen_tau_e_3=Float_t(gen_tau_e_3);
newgen_tau_pt_3=Float_t(gen_tau_pt_3);
newgen_tau_eta_3=Float_t(gen_tau_eta_3);
newgen_tau_phi_3=Float_t(gen_tau_phi_3);
newptGenVhad=Float_t(ptGenVhad);
newetaGenVhad=Float_t(etaGenVhad);
newphiGenVhad=Float_t(phiGenVhad);
newmassGenVhad=Float_t(massGenVhad);
newptGenV_2=Float_t(ptGenV_2);
newetaGenV_2=Float_t(etaGenV_2);
newphiGenV_2=Float_t(phiGenV_2);
newmassGenV_2=Float_t(massGenV_2);
newptGenV_3=Float_t(ptGenV_3);
newetaGenV_3=Float_t(etaGenV_3);
newphiGenV_3=Float_t(phiGenV_3);
newmassGenV_3=Float_t(massGenV_3);
newptGenVlep=Float_t(ptGenVlep);
newetaGenVlep=Float_t(etaGenVlep);
newphiGenVlep=Float_t(phiGenVlep);
newmassGenVlep=Float_t(massGenVlep);
newptGenVlep_2=Float_t(ptGenVlep_2);
newetaGenVlep_2=Float_t(etaGenVlep_2);
newphiGenVlep_2=Float_t(phiGenVlep_2);
newmassGenVlep_2=Float_t(massGenVlep_2);
newptGenVlep_3=Float_t(ptGenVlep_3);
newetaGenVlep_3=Float_t(etaGenVlep_3);
newphiGenVlep_3=Float_t(phiGenVlep_3);
newmassGenVlep_3=Float_t(massGenVlep_3);
newptq11=Float_t(ptq11);
newetaq11=Float_t(etaq11);
newphiq11=Float_t(phiq11);
newmassq11=Float_t(massq11);
newptq21=Float_t(ptq21);
newetaq21=Float_t(etaq21);
newphiq21=Float_t(phiq21);
newmassq21=Float_t(massq21);
newptq31=Float_t(ptq31);
newetaq31=Float_t(etaq31);
newphiq31=Float_t(phiq31);
newmassq31=Float_t(massq31);
newptq12=Float_t(ptq12);
newetaq12=Float_t(etaq12);
newphiq12=Float_t(phiq12);
newmassq12=Float_t(massq12);
newptq22=Float_t(ptq22);
newetaq22=Float_t(etaq22);
newphiq22=Float_t(phiq22);
newmassq22=Float_t(massq22);
newptq32=Float_t(ptq32);
newetaq32=Float_t(etaq32);
newphiq32=Float_t(phiq32);
newmassq32=Float_t(massq32);

MJ_j_18=-99;
MJJ_j_18=-99;
MJ_j_10=-99;
MJJ_j_10=-99;
MJ_j_12=-99;
MJJ_j_12=-99;
MJ_j_14=-99;
MJJ_j_14=-99;
MJ_j_16=-99;
MJJ_j_16=-99;

TLorentzVector AK41,AK42,AK43,AK44,AK45,AK46,LMJ_j,LMJJ_j;
AK41.SetPtEtaPhiM(0,-99,-99,-99);
AK42.SetPtEtaPhiM(0,-99,-99,-99);
AK43.SetPtEtaPhiM(0,-99,-99,-99);
AK44.SetPtEtaPhiM(0,-99,-99,-99);
AK45.SetPtEtaPhiM(0,-99,-99,-99);
AK46.SetPtEtaPhiM(0,-99,-99,-99);
LMJ_j.SetPtEtaPhiM(0,-99,-99,-99);
LMJJ_j.SetPtEtaPhiM(0,-99,-99,-99);

ak4Ptex1=-99;
ak4Etaex1=-99;
ak4Phiex1=-99;
ak4Eex1=-99;
ak4Ptex2=-99;
ak4Etaex2=-99;
ak4Phiex2=-99;
ak4Eex2=-99;

double deltaRak4sj[8]={0};
double deltaRak4sj2[8]={0};
double deltaRak4sj3[8]={0};
double deltaRak4sj4[8]={0};

for(Int_t i=0;i<8;i++){
deltaRak4sj[i]=sqrt(pow(fabs(ak4jet_eta[i]-Etaj),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-Phij),2*Pi-fabs(ak4jet_phi[i]-Phij)),2));
deltaRak4sj2[i]=sqrt(pow(fabs(ak4jet_eta[i]-Etaj_2),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-Phij_2),2*Pi-fabs(ak4jet_phi[i]-Phij_2)),2));
       if(Nj8==3){
deltaRak4sj3[i]=sqrt(pow(fabs(ak4jet_eta[i]-Etaj_3),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-Phij_3),2*Pi-fabs(ak4jet_phi[i]-Phij_3)),2));
           }
       else{
deltaRak4sj3[i]=10000;
           }     
       if(Nj8==4){
deltaRak4sj4[i]=sqrt(pow(fabs(ak4jet_eta[i]-Etaj_4),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-Phij_4),2*Pi-fabs(ak4jet_phi[i]-Phij_4)),2));
           }
       else{
deltaRak4sj4[i]=10000;
           }
}

num_ak4jetsex=0;
num_ak4jetsin=0;

Double_t DR=0.8;

for(Int_t ii=0; ii<8; ii++) {
bool cutsj=0;
bool cutsj2=0;
bool cutsj3=0;
bool cutsj4=0;
bool cutsjin=0; 
bool cutsjin2=0;
bool cutsjin3=0;
bool cutsjin4=0;
bool sjout=1,sjin=0;

cutsj=(ak4jet_pt[ii]>0)&&(deltaRak4sj[ii]>DR);
cutsjin=(ak4jet_pt[ii]>0)&&(deltaRak4sj[ii]<DR);

cutsj2=(ak4jet_pt[ii]>0)&&(deltaRak4sj2[ii]>DR);
cutsjin2=(ak4jet_pt[ii]>0)&&(deltaRak4sj2[ii]<DR);

cutsj3=(ak4jet_pt[ii]>0)&&(deltaRak4sj3[ii]>DR);
cutsjin3=(ak4jet_pt[ii]>0)&&(deltaRak4sj3[ii]<DR);

cutsj4=(ak4jet_pt[ii]>0)&&(deltaRak4sj4[ii]>DR);
cutsjin4=(ak4jet_pt[ii]>0)&&(deltaRak4sj4[ii]<DR);

sjout*=cutsj*cutsj2*cutsj3*cutsj4;
sjin+=cutsjin+cutsjin2+cutsjin3+cutsjin4;

if(sjout==1) 
              {num_ak4jetsex++;
                   if(num_ak4jetsex==1) {
         AK41.SetPtEtaPhiE(ak4jet_pt[ii],ak4jet_eta[ii],ak4jet_phi[ii],ak4jet_e[ii]);
                               ak4Ptex1=ak4jet_pt[ii];
                               ak4Etaex1=ak4jet_eta[ii];
                               ak4Phiex1=ak4jet_phi[ii];
                               ak4Eex1=ak4jet_e[ii];
                                        }
                   if(num_ak4jetsex==2) {
         AK42.SetPtEtaPhiE(ak4jet_pt[ii],ak4jet_eta[ii],ak4jet_phi[ii],ak4jet_e[ii]);
                               ak4Ptex2=ak4jet_pt[ii];
                               ak4Etaex2=ak4jet_eta[ii];
                               ak4Phiex2=ak4jet_phi[ii];
                               ak4Eex2=ak4jet_e[ii];
                                        }
                   if(num_ak4jetsex==3) {
         AK43.SetPtEtaPhiE(ak4jet_pt[ii],ak4jet_eta[ii],ak4jet_phi[ii],ak4jet_e[ii]);
                                        }
                   if(num_ak4jetsex==4) {
         AK44.SetPtEtaPhiE(ak4jet_pt[ii],ak4jet_eta[ii],ak4jet_phi[ii],ak4jet_e[ii]);
                                        }
                   if(num_ak4jetsex==5) {
         AK45.SetPtEtaPhiE(ak4jet_pt[ii],ak4jet_eta[ii],ak4jet_phi[ii],ak4jet_e[ii]);
                                        }
                   if(num_ak4jetsex==6) {
         AK46.SetPtEtaPhiE(ak4jet_pt[ii],ak4jet_eta[ii],ak4jet_phi[ii],ak4jet_e[ii]);
                                        }
              }
if(sjin>=1) num_ak4jetsin++;
}

TLorentzVector AK81,AK82,AK83,AK84;
AK81.SetPtEtaPhiM(0,-99,-99,-99);
AK82.SetPtEtaPhiM(0,-99,-99,-99);
AK83.SetPtEtaPhiM(0,-99,-99,-99);
AK84.SetPtEtaPhiM(0,-99,-99,-99);
MJJj=-99;
MJJjj=-99;
MJJJj=-99;
MJJJJ=-99;
double AK8MMIN_n_phi,AK8MMIN_n_eta;
double DR1,DR2,DR3,DR4,DR5,DR6;
DR1=99;
DR2=99;
DR3=99;
DR4=99;
DR5=99;
DR6=99;
double DR01,DR02,DR03,DR04,DR05;
DR01=1.8;
DR02=1.0;
DR03=1.2;
DR04=1.4;
DR05=1.6;

AK8MMIN_n_phi=((Phij_min+Pi)>2*Pi)*(Phij_min-Pi)+((Phij_min+Pi)<2*Pi)*(Phij_min+Pi);
AK8MMIN_n_eta=-Etaj_min;
if(Nj8==2)
{

if(num_ak4jetsex==6){
DR1=sqrt(pow(fabs(AK41.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR2=sqrt(pow(fabs(AK42.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR3=sqrt(pow(fabs(AK43.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR4=sqrt(pow(fabs(AK44.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR5=sqrt(pow(fabs(AK45.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR6=sqrt(pow(fabs(AK46.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
                    }
if(num_ak4jetsex==5){
DR1=sqrt(pow(fabs(AK41.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR2=sqrt(pow(fabs(AK42.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR3=sqrt(pow(fabs(AK43.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR4=sqrt(pow(fabs(AK44.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR5=sqrt(pow(fabs(AK45.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
                    }
if(num_ak4jetsex==4){
DR1=sqrt(pow(fabs(AK41.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR2=sqrt(pow(fabs(AK42.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR3=sqrt(pow(fabs(AK43.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR4=sqrt(pow(fabs(AK44.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
                    }
if(num_ak4jetsex==3){
DR1=sqrt(pow(fabs(AK41.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR2=sqrt(pow(fabs(AK42.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR3=sqrt(pow(fabs(AK43.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
                    }
if(num_ak4jetsex==2){
DR1=sqrt(pow(fabs(AK41.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
DR2=sqrt(pow(fabs(AK42.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
                    }
if(num_ak4jetsex==1){
DR1=sqrt(pow(fabs(AK41.Eta()-AK8MMIN_n_eta),2)+pow(TMath::Min(fabs(AK41.Phi()-AK8MMIN_n_phi),2*Pi-fabs(AK41.Phi()-AK8MMIN_n_phi)),2));
                    }

AK82.SetPtEtaPhiM(PTj_min,Etaj_min,Phij_min,Mj_min);
AK81.SetPtEtaPhiM(PTj_max,Etaj_max,Phij_max,Mj_max);
//order by mass                   
LMJ_j=AK81;
LMJJ_j=AK81+AK82;
if(DR1<DR01){
LMJ_j=LMJ_j+AK41;
LMJJ_j=LMJJ_j+AK41;
           }
if(DR2<DR01){
LMJ_j=LMJ_j+AK42;
LMJJ_j=LMJJ_j+AK42;
           }
if(DR3<DR01){
LMJ_j=LMJ_j+AK43;
LMJJ_j=LMJJ_j+AK43;
           }
if(DR4<DR01){
LMJ_j=LMJ_j+AK44;
LMJJ_j=LMJJ_j+AK44;
           }
if(DR5<DR01){
LMJ_j=LMJ_j+AK45;
LMJJ_j=LMJJ_j+AK45;
           }
if(DR6<DR01){
LMJ_j=LMJ_j+AK46;
LMJJ_j=LMJJ_j+AK46;
           }
MJ_j_18=LMJ_j.M();
MJJ_j_18=LMJJ_j.M();


LMJ_j=AK81;
LMJJ_j=AK81+AK82;
if(DR1<DR02){
LMJ_j=LMJ_j+AK41;
LMJJ_j=LMJJ_j+AK41;
           }
if(DR2<DR02){
LMJ_j=LMJ_j+AK42;
LMJJ_j=LMJJ_j+AK42;
           }
if(DR3<DR02){
LMJ_j=LMJ_j+AK43;
LMJJ_j=LMJJ_j+AK43;
           }
if(DR4<DR02){
LMJ_j=LMJ_j+AK44;
LMJJ_j=LMJJ_j+AK44;
           }
if(DR5<DR02){
LMJ_j=LMJ_j+AK45;
LMJJ_j=LMJJ_j+AK45;
           }
if(DR6<DR02){
LMJ_j=LMJ_j+AK46;
LMJJ_j=LMJJ_j+AK46;
           }
MJ_j_10=LMJ_j.M();
MJJ_j_10=LMJJ_j.M();


LMJ_j=AK81;
LMJJ_j=AK81+AK82;
if(DR1<DR03){
LMJ_j=LMJ_j+AK41;
LMJJ_j=LMJJ_j+AK41;
           }
if(DR2<DR03){
LMJ_j=LMJ_j+AK42;
LMJJ_j=LMJJ_j+AK42;
           }
if(DR3<DR03){
LMJ_j=LMJ_j+AK43;
LMJJ_j=LMJJ_j+AK43;
           }
if(DR4<DR03){
LMJ_j=LMJ_j+AK44;
LMJJ_j=LMJJ_j+AK44;
           }
if(DR5<DR03){
LMJ_j=LMJ_j+AK45;
LMJJ_j=LMJJ_j+AK45;
           }
if(DR6<DR03){
LMJ_j=LMJ_j+AK46;
LMJJ_j=LMJJ_j+AK46;
           }
MJ_j_12=LMJ_j.M();
MJJ_j_12=LMJJ_j.M();


LMJ_j=AK81;
LMJJ_j=AK81+AK82;
if(DR1<DR04){
LMJ_j=LMJ_j+AK41;
LMJJ_j=LMJJ_j+AK41;
           }
if(DR2<DR04){
LMJ_j=LMJ_j+AK42;
LMJJ_j=LMJJ_j+AK42;
           }
if(DR3<DR04){
LMJ_j=LMJ_j+AK43;
LMJJ_j=LMJJ_j+AK43;
           }
if(DR4<DR04){
LMJ_j=LMJ_j+AK44;
LMJJ_j=LMJJ_j+AK44;
           }
if(DR5<DR04){
LMJ_j=LMJ_j+AK45;
LMJJ_j=LMJJ_j+AK45;
           }
if(DR6<DR04){
LMJ_j=LMJ_j+AK46;
LMJJ_j=LMJJ_j+AK46;
           }
MJ_j_14=LMJ_j.M();
MJJ_j_14=LMJJ_j.M();


LMJ_j=AK81;
LMJJ_j=AK81+AK82;
if(DR1<DR05){
LMJ_j=LMJ_j+AK41;
LMJJ_j=LMJJ_j+AK41;
           }
if(DR2<DR05){
LMJ_j=LMJ_j+AK42;
LMJJ_j=LMJJ_j+AK42;
           }
if(DR3<DR05){
LMJ_j=LMJ_j+AK43;
LMJJ_j=LMJJ_j+AK43;
           }
if(DR4<DR05){
LMJ_j=LMJ_j+AK44;
LMJJ_j=LMJJ_j+AK44;
           }
if(DR5<DR05){
LMJ_j=LMJ_j+AK45;
LMJJ_j=LMJJ_j+AK45;
           }
if(DR6<DR05){
LMJ_j=LMJ_j+AK46;
LMJJ_j=LMJJ_j+AK46;
           }
MJ_j_16=LMJ_j.M();
MJJ_j_16=LMJJ_j.M();


AK81.SetPtEtaPhiM(PTj_max,Etaj_max,Phij_max,Mj_max);
AK82.SetPtEtaPhiM(PTj_min,Etaj_min,Phij_min,Mj_min);
//order by mass
      MJJJj=-99;
    if(AK41.Pt()>0)
          {
      MJJj=(AK81+AK82+AK41).M();
          }
    else
          {
      MJJj=(AK81+AK82).M();
          }
    
    if(AK41.Pt()>0&&AK42.Pt()>0)
          {
      MJJjj=(AK81+AK82+AK41+AK42).M();
          }
    else if(AK41.Pt()>0)
          {
      MJJjj=(AK81+AK82+AK41).M();
          }
    else
      MJJjj=(AK81+AK82).M();
}

if(Nj8==3)
{
AK81.SetPtEtaPhiM(PTj_max,Etaj_max,Phij_max,Mj_max);
AK82.SetPtEtaPhiM(PTj_mid,Etaj_mid,Phij_mid,Mj_mid);
AK83.SetPtEtaPhiM(PTj_min,Etaj_min,Phij_min,Mj_min);
    if(AK41.Pt()>0)
          {
      MJJJj=(AK81+AK82+AK83+AK41).M();
          }
    else
          {
      MJJJj=(AK81+AK82+AK83).M();
          }

    if(AK41.Pt()>0)
          {
      MJJj=(AK81+AK82+AK41).M();
          }
    else
          {
      MJJj=(AK81+AK82).M();
          }

    if(AK41.Pt()>0&&AK42.Pt()>0)
          {
      MJJjj=(AK81+AK82+AK41+AK42).M();
          }
    else if(AK41.Pt()>0)
          {
      MJJjj=(AK81+AK82+AK41).M();
          }
    else
      MJJjj=(AK81+AK82).M();
}

if(Nj8==4)
{
AK81.SetPtEtaPhiM(PTj_max,Etaj_max,Phij_max,Mj_max);
AK82.SetPtEtaPhiM(PTj_mid,Etaj_mid,Phij_mid,Mj_mid);
AK83.SetPtEtaPhiM(PTj_min,Etaj_min,Phij_min,Mj_min);
AK84.SetPtEtaPhiM(PTj_4,Etaj_4,Phij_4,Mj_4);
    if(AK41.Pt()>0)
          {
      MJJJj=(AK81+AK82+AK83+AK41).M();
          }
    else  
          {
      MJJJj=(AK81+AK82+AK83).M();
          }
    
    if(AK41.Pt()>0)
          {
      MJJj=(AK81+AK82+AK41).M();
          }
    else  
          {
      MJJj=(AK81+AK82).M();
          }
    
    if(AK41.Pt()>0&&AK42.Pt()>0)
          {
      MJJjj=(AK81+AK82+AK41+AK42).M();
          }
    else if(AK41.Pt()>0)
          {
      MJJjj=(AK81+AK82+AK41).M();
          }
    else
      MJJjj=(AK81+AK82).M();
    MJJJJ=(AK81+AK82+AK83+AK84).M();
}

MJJ=(AK81+AK82).M();
Mass2j1j2 = Float_t(massww[0]*massww[0]);
Mass2j3j1 = Float_t(massww[1]*massww[1]);
Mass2j2j3 = Float_t(massww[2]*massww[2]);
Mw1w2=Float_t(massww[0]);
Mw1w3=Float_t(massww[1]);
Mw2w3=Float_t(massww[2]);

        mtVlepnew         = Float_t(sqrt(2*ptlep1*met*(1.0-cos(philep1-metPhi))));
        MTVlep            = Float_t(sqrt(2*ptlep1*MET_et*(1.0-cos(philep1-MET_phi))));
		Double_t deltaRWhadGen = sqrt(pow(etaGenVhad-jetAK8puppi_eta,2) + pow(phiGenVhad-jetAK8puppi_phi,2));

		//Weight Calculation
		Int_t bin = hR1->FindBin(npT);
		pileupWeight = hR1->GetBinContent(bin);		

		eff_and_pu_Weight = 0;
		eff_and_pu_Weight1 = 0;
		if(IsData>0) {
			if(npT < weights_pu1.size()){
				eff_and_pu_Weight = weights_pu1[npT];
			}
			if(npT < weights_pu2.size()){
				eff_and_pu_Weight1 = weights_pu2[npT];
			}
		}

                trigger_eff=1.0;
		if(theWeight>0) nn=1;
		else nn= -1;
		if(npp>0) lumiWeight=Identical_lumiWeight/(npp-nmm);
		else lumiWeight=Identical_lumiWeight/nentries;
		weight_nobtag=lumiWeight*triggerWeight*eff_and_pu_Weight*nn*trigger_eff;
		//weight=lumiWeight*triggerWeight*pileupWeight*nn;
		if (IsData>1 ) weight_nobtag = weight_nobtag*1.21;

		//lumiWeight=Identical_lumiWeight;
		//if(npp>0) weight=lumiWeight*triggerWeight*pileupWeight/(npp-nmm)*nn*0.04024;//0.00559;
		//else weight=lumiWeight*triggerWeight*pileupWeight/nentries*0.04024;//0.00559;
		if ( IsData==0 ) weight_nobtag=1;
		//Weight Calculation Done


      IDweight=1.0;
      IDweightISO=1.0;
      IDweighttrk=1.0;

        ToppTweight=1.0;

//--BSF----------------------------
      btagweight_center=1.0;
      btagweight_up=1.0;
      btagweight_down=1.0;
      
                if(theWeight>0) nn=1;
                else nn= -1;
                if(npp>0) lumiWeight=Identical_lumiWeight/(npp-nmm);
                else lumiWeight=Identical_lumiWeight/nentries;
                weight=lumiWeight*triggerWeight*eff_and_pu_Weight*nn*trigger_eff*btagweight_center*IDweight*IDweightISO*IDweighttrk*ToppTweight;
                if (IsData==2 ) weight = weight*1.21*1.06684;
                if (IsData==3 ) weight = weight*1.21*1.046;
                if (IsData==4 ) weight = weight*1.21*0.99246;
                if (IsData==5 ) weight = weight*1.21*0.9155;
                if (IsData==6 ) weight = weight*1.21*0.8093;
                if (IsData==7 ) weight = weight*1.21*0.6498;
                if (IsData==8 ) weight = weight*1.21*0.4843;
                if ( IsData==0 ) weight=1;

		//number of bjet calculation
		num_bJet=0.;
		num_bJet_loose=0.;
		num_bJet_tight=0.;
		for(Int_t i=0; i<8; i++)  {
                        deltaRAK4AK8_new[i]=0.;
                        deltaRAK4AK8_new[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.8484 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 ) {num_bJet=num_bJet+1;}
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.5426 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 ) {num_bJet_loose=num_bJet_loose+1;}
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.9535 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 ) {num_bJet_tight=num_bJet_tight+1;}
		}
		nbtag=num_bJet;
		//number of bjet calculation Done

		Int_t nLooseLep=nLooseEle+nLooseMu;//the tight Lep included

		Double_t isAnaHP=1.;
		Double_t isAnaLP=1.;
		Double_t isAnaNP=1.;
		Double_t isTTBarControl=1.;
		Int_t tmp_categoryID_channel=0;
		if( channelname=="el" ){
		tmp_categoryID_channel=-1;// -1 for el; 1 for mu

			//HP: 0<jetAK8puppi_tau21<=0.5;
			if (isAnaHP>0 && lep==11 && nLooseLep==1){ nID_e = nID_e+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptlep1>55 && fabs(etalep1)<2.5){ npt_e = npt_e+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && MET_et>80) { nmet_e = nmet_e+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptVlepJEC > 200.) { nptVlepJEC = nptVlepJEC +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && PTj>200 && fabs(jetAK8puppi_eta)<2.4 && IDLoose>0 ){ nptVhad = nptVhad+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && num_bJet<1){ nnum_bJet_e = nnum_bJet_e +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && deltaRlepjet>pi_2) {n_deltaRlepjet = n_deltaRlepjet+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && fabs(delPhijetmet) >2.0)  {n_delPhijetmet = n_delPhijetmet+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && fabs(delPhijetlep)>2.0) {n_delPhijetlep = n_delPhijetlep +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && jetAK8puppi_tau21>0. && jetAK8puppi_tau21<=0.45) {ntau = ntau+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0)
			{
				nmassVhad = nmassVhad +1;
				yields = yields + weight;
				(*file_cutflow)<<"event:"<<event<<endl;
			} else{ isAnaHP=-1; }

			//LP: 0.5<jetAK8puppi_tau21<=0.75;
			if ( lep==11 && nLooseLep==1 && ptlep1>55 && fabs(etalep1)<2.5 && MET_et>80 && ptVlepJEC > 200.  && PTj>200 && fabs(jetAK8puppi_eta)<2.4 && IDLoose>0 && num_bJet<1 && deltaRlepjet>pi_2 && fabs(delPhijetmet) >2.0 && fabs(delPhijetlep)>2.0 && jetAK8puppi_tau21>0.45 && jetAK8puppi_tau21<=0.75)// && (( jetAK8puppi_sdcorr >0 &&  jetAK8puppi_sdcorr< 150 )) )
			{ isAnaLP=1.; } 
			else{ isAnaLP=-1.; }
			//NP: 0.75<jetAK8puppi_tau21
			if ( lep==11 && nLooseLep==1 && ptlep1>55 && fabs(etalep1)<2.5 && MET_et>80 && ptVlepJEC > 200.  && PTj>200 && fabs(jetAK8puppi_eta)<2.4 && IDLoose>0 && num_bJet<1 && deltaRlepjet>pi_2 && fabs(delPhijetmet) >2.0 && fabs(delPhijetlep)>2.0 && jetAK8puppi_tau21>0.75)// && (( jetAK8puppi_sdcorr >0 &&  jetAK8puppi_sdcorr< 150 )) )
			{ isAnaNP=1.; } 
			else{ isAnaNP=-1.; }


			//TTbar control
			if ( lep==11 && nLooseLep==1 && ptlep1>55 && fabs(etalep1)<2.5 && MET_et>80 && ptVlepJEC > 200. && PTj>200 && fabs(jetAK8puppi_eta)<2.4 && IDLoose>0  && num_bJet>0)// &&  jetAK8puppi_sdcorr>0 && jetAK8puppi_sdcorr <150)
			{ isTTBarControl=1.; } 
			else{ isTTBarControl=-1.; }
		}
		else if( channelname=="mu" ){
		tmp_categoryID_channel=1;// -1 for el; 1 for mu
			//HP: 0<jetAK8puppi_tau21<=0.5;
               	        if (isAnaHP>0 &&(HLT_Mu1>0 || HLT_Mu2>0 || HLT_Mu3>0 || HLT_Mu4>0 || HLT_Mu5>0 || HLT_Mu6>0 || HLT_Mu7>0 || HLT_Mu8>0 || HLT_Mu9>0)) { isAnaHP=1; } else{ isAnaHP=-1; }

			if ( ((Nj8==2&&MJJ>1000)||(Nj8==3&&MJJJ>1000&&fabs(jetAK8puppi_eta_3)<2.4)||(Nj8==4&&MJJJ>1000&&fabs(jetAK8puppi_eta_4)<2.4))  && isAnaHP>0 && PTj>400 && PTj_2>200 && Mj>40 && Mj_2>40 && fabs(jetAK8puppi_eta)<2.4 && fabs(jetAK8puppi_eta_2)<2.4 && passFilter_HBHE>0 && passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_EEBadSc>0){ isAnaHP=1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && jetAK8puppi_tau21>0. && jetAK8puppi_tau21<=0.45) {ntau = ntau+1;} else{ isAnaHP=-1; }
			if (isAnaHP>0)// && (( jetAK8puppi_sdcorr >0&& jetAK8puppi_sdcorr< 150 )))// && m_lvj>100 && m_lvj<5000 )
			{ 
				nmassVhad = nmassVhad +1; 
				(*file_cutflow)<<"event:"<<event<<endl;
			} else{ isAnaHP=-1; }

			//LP: 0.5<jetAK8puppi_tau21<=0.75;
			 if ( ((Nj8==2&&MJJ>1000)||(Nj8==3&&MJJJ>1000&&fabs(jetAK8puppi_eta_3)<2.4)||(Nj8==4&&MJJJ>1000&&fabs(jetAK8puppi_eta_4)<2.4)) && PTj>400 && PTj_2>200 && Mj>40 && Mj_2>40 && (HLT_Mu1>0 || HLT_Mu2>0 || HLT_Mu3>0 || HLT_Mu4>0 || HLT_Mu5>0 || HLT_Mu6>0 || HLT_Mu7>0 || HLT_Mu8>0 || HLT_Mu9>0) && fabs(jetAK8puppi_eta)<2.4 && jetAK8puppi_tau21>0.45 && jetAK8puppi_tau21<=0.75&& fabs(jetAK8puppi_eta_2)<2.4 && passFilter_HBHE>0 && passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_EEBadSc>0)// && (( jetAK8puppi_sdcorr >0&& jetAK8puppi_sdcorr< 150 )))

			{ isAnaLP=1.; } 
			else{ isAnaLP=-1.; }

			//NP: 0.75<jetAK8puppi_tau21;
		
			if ( ((Nj8==2&&MJJ>1000)||(Nj8==3&&MJJJ>1000&&fabs(jetAK8puppi_eta_3)<2.4)||(Nj8==4&&MJJJ>1000&&fabs(jetAK8puppi_eta_4)<2.4))   && PTj>400 && PTj_2>200 && Mj>40 && Mj_2>40 && (HLT_Mu1>0 || HLT_Mu2>0 || HLT_Mu3>0 || HLT_Mu4>0 || HLT_Mu5>0 || HLT_Mu6>0 || HLT_Mu7>0 || HLT_Mu8>0 || HLT_Mu9>0) && fabs(jetAK8puppi_eta)<2.4 && jetAK8puppi_tau21>0.75&& fabs(jetAK8puppi_eta_2)<2.4 && passFilter_HBHE>0 && passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_EEBadSc>0)

			{ isAnaNP=1.; } 
			else{ isAnaNP=-1.; }

			//TTbar control
			if (((Nj8==2&&MJJ>1000)||(Nj8==3&&MJJJ>1000&&fabs(jetAK8puppi_eta_3)<2.4)||(Nj8==4&&MJJJ>1000&&fabs(jetAK8puppi_eta_4)<2.4))  && PTj>400 && PTj_2>200 && Mj>40 && Mj_2>40 && (HLT_Mu1>0 || HLT_Mu2>0 || HLT_Mu3>0 || HLT_Mu4>0 || HLT_Mu5>0 || HLT_Mu6>0 || HLT_Mu7>0 || HLT_Mu8>0 || HLT_Mu9>0) && fabs(jetAK8puppi_eta)<2.4&& fabs(jetAK8puppi_eta_2)<2.4 && passFilter_HBHE>0 && passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_EEBadSc>0)// && jetAK8puppi_sdcorr>0 && jetAK8puppi_sdcorr <150)
			{ isTTBarControl=1.; } 
			else{ isTTBarControl=-1.; }

		}else{
			cout<<"We don't know channelname:"<<channelname<<endl;
		}

		Int_t tmp_categoryID_eventselection=0;
		if(isAnaHP>0)tmp_categoryID_eventselection=1;
		else if(isAnaLP>0)tmp_categoryID_eventselection=2;
		else if(isAnaNP>0)tmp_categoryID_eventselection=4;
		else if(isTTBarControl>0)tmp_categoryID_eventselection=3;
		else tmp_categoryID_eventselection=100;

		CategoryID=tmp_categoryID_channel* tmp_categoryID_eventselection;

		isMatch=1.;
		if(deltaRWhadGen >= 0.3) isMatch=-1;
		if(jetAK8puppi_tau21<=0){vTagID=2;}
		else if(jetAK8puppi_tau21>0.45 && jetAK8puppi_tau21<=0.60){vTagID=1;}
		else if(jetAK8puppi_tau21>0.60 && jetAK8puppi_tau21<=0.75){vTagID=0;}
		else if(jetAK8puppi_tau21>0.75 && jetAK8puppi_tau21<=1){vTagID=-1;}
		else {vTagID=-2;}

		if(TMath::Abs(CategoryID)<10) ExTree->Fill();
	}

	if(channelname=="el"){ 
		std::cout << "nID_e" << nID_e << "; npt_e" << npt_e << "; nmet_e" << nmet_e << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<<"; nnum_bJet_e" << nnum_bJet_e <<"; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad << "; number_qq" << number_qq << "; yields " << yields << std::endl;
		(*file_cutflow) << "nID_e" << nID_e << "; npt_e" << npt_e << "; nmet_e" << nmet_e << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<<"; nnum_bJet_e" << nnum_bJet_e <<"; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad << "; number_qq" << number_qq << std::endl;
	}
	if(channelname=="mu"){
		std::cout << "nID_mu" << nID_mu << "; npt_mu" << npt_mu << "; nmet_mu" << nmet_mu << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<< "; nnum_bJet_mu" << nnum_bJet_mu << "; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad<< "; number_qq" << number_qq << std::endl;
		(*file_cutflow) << "nID_mu" << nID_mu << "; npt_mu" << npt_mu << "; nmet_mu" << nmet_mu << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<< "; nnum_bJet_mu" << nnum_bJet_mu << "; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad<< "; number_qq" << number_qq << std::endl;
	}

}
