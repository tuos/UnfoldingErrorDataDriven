#include <iostream>
#include <stdio.h>

#include <TRandom2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include "TRandom3.h"

#include "utilities.h"
#include "bayesianUnfold.h"
#include "prior.h"
#include "generateSmoothMatrix.h"

using namespace std;


//==============================================================================
// Unfolding Ying Lu 08 07 11
// Update Yen-Jie Lee 06.22.12
//==============================================================================

void pAUnfoldSmearedMatrix(int algo= 3,bool useSpectraFromFile=0, bool useMatrixFromFile=0, int doToy = 0, int isMC = 0,char *spectraFileName = (char*)"ppb_spectra_akPu3PF.root",double recoJetPtCut = 30.,double trackMaxPtCut = 0, int nBayesianIter = 4, int doBjets=0, int doTrigEffCorr=0) // algo 2 =akpu2 ; 3 =akpu3 ; 4 =akpu4 ;1 = icpu5
{
  
  gStyle->SetErrorX(0.);
  gStyle->SetPaintTextFormat("3.2f");
  gStyle->SetOptLogz(1);
  gStyle->SetPadRightMargin(0.13);	
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);


  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  //TString coll ;
      coll= getenv("coll") ; 
  const float ppsigma=70.;  //units in mb
  const float pPbSigma=2.061;  //units in b
//  const float Lumi = 15.78 ; //units in nb
   float Lumi ;
    if(coll=="PPb") Lumi = 20.7 ; //units in nb
    else Lumi = 14. ;
  const float frac = 0.85 ;
  const float nMB = 2.60259184640000000e+10 ; 
  const float Ncoll = 6.9 ; // for inclusive pPb 
//  const float Ncoll = 7.5 ; // for pPb in 0-90% 
   float etamin ; //= -1.0 ;
   float etamax ; //= 1.0 ;
   etamin=atof(getenv("MIN"));
   etamax=atof(getenv("MAX"));
   int ieta ;
  if(coll=="PPb" ) {
   if(etamin==-1 && etamax == 1.) ieta = 0 ;
   if(etamin==-1.5 && etamax == -1.0) ieta = 1 ;
   if(etamin==-1.0 && etamax == -0.5) ieta = 2 ;
   if(etamin==-0.5 && etamax == 0.5) ieta = 3 ;
   if(etamin==0.5 && etamax == 1.0) ieta = 4 ;
   if(etamin==1.0 && etamax == 1.5) ieta = 5 ;
   if(etamin==1.5 && etamax == 2.0) ieta = 6 ;
   if(etamin==2.0 && etamax == 2.5) ieta = 7 ;
  }  
  else {
   if(etamin==-1 && etamax == 1.) ieta = 0 ;
   if(etamin==-2.5 && etamax == -1.0) ieta = 1 ;
   if(etamin==-2.0 && etamax == -1.5) ieta = 2 ;
   if(etamin==-1.5 && etamax == -1.0) ieta = 3 ;
   if(etamin==-1.0 && etamax == -0.5) ieta = 4 ;
   if(etamin==-0.5 && etamax == 0.5) ieta = 5 ;
   if(etamin==0.5 && etamax == 1.0) ieta = 6 ;
   if(etamin==1.0 && etamax == 1.5) ieta = 7 ;
  }
  cout << "eta min = " << etamin << " -  max =" <<etamax <<endl ;
  const bool SavePlot=kFALSE;
  bool doParameterizedMatrix = 0;
  const bool DoSmear=kFALSE ;
   TF1 * fsmear = new TF1("fsmear","[0]/pow(x,[1])",30,700);  
   if(algo==2) fsmear->SetParameters(0.4183, 0.4244);
   if(algo==3) fsmear->SetParameters(1.052, 0.5261);
   TF1 * fgaus = new TF1("fgaus", "[0]*exp(-0.5*((x-[1])/[2])**2)", -20, 20);
   fgaus->SetParameters(1, 0, 1);
  
double fpara[8][5];
double fparaerr[8][5];
  double fppb[8][5]={{4.00741e+12 , -78.6815 ,  -6.38764 , 4.24115,   0}, //! -1 to 1
                     {1.30459e+11 , -51.4578 ,  -5.77698 , 7.11119,   -1.25}, //! -1.5 to -1
                     {1.16972e+12 , -81.4762 ,  -6.18544 , 4.05456,   -0.75},
                     {7.11276e+11 , -72.4556 ,  -6.08591 , 4.7646,   0},
                     {7.70754e+11 , -75.3034 ,  -6.10602 , 4.62161,   0.75},
                     {1.14373e+11 , -52.6737 ,  -5.73268 , 7.42508,   1.25},
                     {4.91178e+10 , -48.3787 ,  -5.56167 , 8.35963,   1.75},
                     {1.98136e+09 , -26.2293 ,  -4.82503 , 11.5059,   2.25} };
 double fppberr[8][5]={{4.3985e+10 , 0.331003 ,  0.331003 , 0.331003,   0},
                        {4.25335e+09 , 1.01527 ,  1.01527 , 1.01527,   0},
                        {2.97339e+10 , 0.77854 ,  0.77854 , 0.77854,   0},
                        {1.15465e+10 , 0.530106 ,  0.530106 , 0.530106,   0},
                        {1.95303e+10 , 0.796188 ,  0.796188 , 0.796188,   0},
                        {3.63768e+09 , 0.992593 ,  0.992593 , 0.992593,   0},
                        {2.03389e+09 , 1.29328 ,  1.29328 , 1.29328,   0},
                        {1.15011e+08 , 1.69676 ,  1.69676 , 1.69676,   0} };
  double fpbp[8][5]={{3.45139e+12 , -76.868 ,  -6.35762 , 4.58932,   0},
                     {6.13561e+08 , -8.32984 ,  -4.60887 , 11.7227,   -2.25},
                     {6.35639e+10 , -51.5085 ,  -5.61183 , 8.2267,   -1.75},
                     {6.18049e+10 , -44.8095 ,  -5.61109 , 8.09173,   -1.25},
                     {6.47069e+11 , -73.1118 ,  -6.07436 , 4.7806,   -0.75},
                     {8.25068e+11 , -75.8003 ,  -6.11028 , 4.74121,   0},
                     {7.92363e+11 , -76.2709 ,  -6.10601 , 4.79662,   0.75},
                     {9.08053e+10 , -47.0361 ,  -5.70471 , 7.51127,   1.25} };
    double fpbperr[8][5]={{3.76538e+10 , 0.331963 ,  0.331963 , 0.331963,   0},
                          {3.60767e+07 , 1.71314 ,  1.71314 , 1.71314,   0},
                          {2.61571e+09 , 1.18678 ,  1.18678 , 1.18678,   0},
                          {1.96405e+09 , 0.990336 ,  0.990336 , 0.990336,   0},
                          {1.70688e+11 , 3.91135 ,  3.91135 , 3.91135,   0},
                          {1.31958e+10 , 0.513648 ,  0.513648 , 0.513648,   0},
                          {2.00762e+10 , 0.817164 ,  0.817164 , 0.817164,   0},
                          {3.02649e+09 , 1.01814 ,  1.01814 , 1.01814,   0} };
   for(int i = 0 ; i < 8 ; i++){
     for(int j = 0 ; j < 5 ; j++){
       if(coll=="PPb") {
          fpara[i][j]=fppb[i][j];
          fparaerr[i][j]=fppberr[i][j];
        }
       else {
          fpara[i][j]=fpbp[i][j];
          fparaerr[i][j]=fpbperr[i][j];
       }
     }
  } 
 // fit on the generator level ture distributions
  TF1 * ftrue = new TF1("ftrue", "[0]*exp([1]/x)*pow(x,[2])*pow(1-x*cosh([4])/4000.,[3])", 30, 800) ;
       for(int ip=0;ip<5;ip++){
        //   ftrue->SetParameter(ip,fpara[ieta][ip]);
           ftrue->SetParameter(ip,fpara[ieta][ip]+fparaerr[ieta][ip]);
           cout << "para i = " << ip <<  " para = " << ftrue->GetParameter(ip)<<endl;
        //   ftrue->SetParameter(ip,(fpara[ieta][ip]-fparaerr[ieta][ip]));
        } 
//   ftrue->SetParameters(1.107e+6, -23.95, -5.507, 3.561, -1.691);
 //  ftrue->SetParameters(1.107e+6+3.257e+5, -23.95+3.07, -5.507+0.062, 3.561+0.895, -1.691+0.161);
//   ftrue->SetParameters(1.107e+6-3.257e+5, -23.95-3.07, -5.507-0.062, 3.561-0.895, -1.691-0.161);

  // input files
  char *fileNamePP_mc = NULL;
  char *fileNamePbPb_mc = NULL;
  char *fileNamePP_data = NULL;
  char *fileNamePbPb_data = NULL;
  

  // pp file needs replacing
  if(doBjets)fileNamePbPb_data = (char*)"/net/hidsk0001/d00/scratch/maoyx/Btag/Unfold/bJetRAA/bFractionMCTemplate_ppPbPb1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2_binomErrors_jet55_wideBin_v2.root";
   else {
    if(coll=="PPb")
    //  fileNamePbPb_data = (char*)Form("/home/maoy/working/pA/RpA/histos/DATA%s%sJetSpectraKurtCombineJetTriggerEtaWeight7EtabinJetIDCutRecoPt.root",coll.Data(), algoName[algo]) ;
    //  fileNamePbPb_data = (char*)Form("/home/maoy/working/pA/RpA/histos/DATA%s%sJetSpectraKurtCombineJetTriggerDogaEtaWeight7EtabinJetIDCutRecoPt.root",coll.Data(), algoName[algo]) ;
      fileNamePbPb_data = (char*)Form("/home/maoy/working/pA/RpA/histos/DATA%s%sJetSpectraKurtCombineJetTriggerEtaWeight7EtabinJetIDCutRecoPt.root",coll.Data(), algoName[algo]) ;
  //    fileNamePbPb_data = (char*)Form("/home/maoy/working/pA/RpA/histos/DATA%s%sJetSpectraKurtCombineJetTriggerEtaWeight7EtabinJetIDCutRawPt.root",coll.Data(), algoName[algo]) ;
 //     fileNamePbPb_data = (char*)Form("/home/maoy/working/pA/RpA/histos/oldForest/DATA%s%sJetSpectraKurtCombineJetTriggerEtaWeight7EtabinNoJetIDCutRecoPt.root",coll.Data(), algoName[algo]) ;
     else 
      fileNamePbPb_data = (char*)Form("/home/maoy/working/pA/RpA/histos/DATA%s%sJetSpectraKurtCombineJetTriggerEtaWeight7EtabinJetIDCutRecoPt.root",coll.Data(), algoName[algo]) ;
  //    fileNamePbPb_data = (char*)Form("/home/maoy/working/pA/RpA/histos/DATA%s%sJetSpectraKurtCombineJetTriggerEtaWeight7EtabinJetIDCutRawPt.root",coll.Data(), algoName[algo]) ;
    //  fileNamePbPb_data = (char*)Form("/home/maoy/working/pA/RpA/histos/oldForest/DATA%s%sJetSpectraKurtCombineJetTriggerEtaWeight7EtabinNoJetIDCutRecoPt.root",coll.Data(), algoName[algo]) ;
    }
  if(doBjets)fileNamePbPb_mc = (char*) "/net/hisrv0001/home/mnguyen/scratch/bTaggingOutput/ntuples/PbPbBMC_pt30by3_ipHICalibCentWeight_noTrig.root";
  else {
  //  if(coll=="PPb")fileNamePbPb_mc = (char*) Form("/cms/store/user/ymao/pA5TEV/Mixing/STARTHI53V27/merged/PPbMCOfficialForestNewVzWeight_ppReco_%s_QCDjetTrigJECv8_JetPt0noIPupperCut.root", algoName[algo]) ;
    if(coll=="PPb")fileNamePbPb_mc = (char*) Form("/cms/store/user/ymao/pA5TEV/Mixing/STARTHI53V27/merged/PPbMCOfficialForestNewVzWeightAddHLT_ppReco_%s_QCDjetTrigJECv8_JetPt0pthatLowerCut.root", algoName[algo]) ;
  else 
 //   fileNamePbPb_mc = (char*) Form("/cms/store/user/ymao/pA5TEV/Mixing/STARTHI53V27/merged/PbPMCOfficialForestNewVzWeight_ppReco_%s_QCDjetTrigJECv19_JetPt0noIPupperCut.root", algoName[algo]) ;
     fileNamePbPb_mc = (char*) Form("/cms/store/user/ymao/pA5TEV/Mixing/STARTHI53V27/merged/PbPMCOfficialForestNewVzWeightAddHLT_ppReco_%s_QCDjetTrigJECv19_JetPt0pthatLowerCut.root", algoName[algo]) ;
    }


  // grab ntuples
//  TFile *infPbPb_mc = new TFile(fileNamePbPb_mc);
//  TFile *infPP_mc = new TFile(fileNamePP_mc);
  

  string bJetString = "Inc";
  if(doBjets) bJetString = "bJets";
  string level = "";
  if(doParameterizedMatrix || DoSmear) level = "SmearedMatrix" ;

  // Output file
  TFile *pbpb_Unfo;
  if (isMC) {
    if(coll=="PPb")
        //pbpb_Unfo = new TFile(Form("/home/maoy/working/pA/RpA/Unfold/outputs/%s_UnfoPriorGen_%sOfficialMC%sNoIDCut_MCJECv8_jtpt%.0f_RecoPtScaleUpEtaBin%.f_%.f_%s_v6.root",coll.Data(), algoName[algo],level.c_str(), recoJetPtCut,-etamax*10, -etamin*10, bJetString.c_str()),"RECREATE");
        pbpb_Unfo = new TFile(Form("/home/tuos/workingpA/v1/Unfold/outputs/%s_UnfoPriorGen_%sOfficialMC%sNoIDCut_MCJECv8_jtpt%.0f_RecoPtScaleUpEtaBin%.f_%.f_%s_MultiUnfolding_v6.root",coll.Data(), algoName[algo],level.c_str(), recoJetPtCut,-etamax*10, -etamin*10, bJetString.c_str()),"RECREATE");
    //else     pbpb_Unfo = new TFile(Form("/home/maoy/working/pA/RpA/Unfold/outputs/%s_UnfoPriorGen_%sOfficialMC%sNoIDCut_MCJECv19_jtpt%.0f_RecoPtScaleUpEtaBin%.f_%.f_%s_v6.root",coll.Data(), algoName[algo],level.c_str(), recoJetPtCut,etamin*10, etamax*10, bJetString.c_str()),"RECREATE");
    else     pbpb_Unfo = new TFile(Form("/home/tuos/workingpA/v1/Unfold/outputs/%s_UnfoPriorGen_%sOfficialMC%sNoIDCut_MCJECv19_jtpt%.0f_RecoPtScaleUpEtaBin%.f_%.f_%s_v6.root",coll.Data(), algoName[algo],level.c_str(), recoJetPtCut,etamin*10, etamax*10, bJetString.c_str()),"RECREATE");
    }
  else {
   if(coll=="PPb")
      pbpb_Unfo  = new TFile(Form("/home/tuos/workingpA/v1/Unfold/outputs/%s_UnfoPriorGen_%sOfficialMC%sJetEtaWeightNoIDCutJECv8_jtpt%.0f_GenPtFitUpEtaBin%.f_%.f_%s_v6.root",coll.Data(), algoName[algo],level.c_str(), recoJetPtCut,-etamax*10, -etamin*10, bJetString.c_str()),"RECREATE");
   else   pbpb_Unfo  = new TFile(Form("/home/tuos/workingpA/v1/Unfold/outputs/%s_UnfoPriorGen_%sOfficialMC%sJetEtaWeightNoIDCutJECv19_jtpt%.0f_GenPtFitUpEtaBin%.f_%.f_%s_v6.root",coll.Data(), algoName[algo],level.c_str(), recoJetPtCut,etamin*10, etamax*10, bJetString.c_str()),"RECREATE");
  } 
// Histograms used by RooUnfold
  UnfoldingHistos *uhist[nbins_cent+1];
//  UnfoldingHistos *uhist[nbins_cent];
		
  // Initialize Histograms   
     for (int i=0;i<=nbins_cent;i++) {
        uhist[i] = new UnfoldingHistos(i);
     }	
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
	
   cout <<"file: " << fileNamePbPb_mc <<endl ;
  // Initialize reweighting functions
  /*
  TCut dataSelection;
  TCut dataSelectionPP;
  TCut TriggerSelectionPP;
  TCut TriggerSelectionPbPb80;

  if(doBjets)dataSelection = "weight*(abs(refparton_flavorForB)==5&&abs(jteta)<2)";
  else dataSelection = "weight*(abs(jteta)<2)";
  */

  if (isMC) cout<<"This is a MC closure test"<<endl;
  else cout<< "This is a data analysis"<<endl;    		     

	     	
  // Setup jet data branches, basically the jet tree branches are assigned to this object when we loop over the events
	
  JetDataPbPb *dataPbPb   = new JetDataPbPb(fileNamePbPb_mc,(char*)"nt"); // PbPb data	
	
  TFile *fSpectra(0);		
  if (useSpectraFromFile||useMatrixFromFile){
    fSpectra = new TFile(spectraFileName,"read");
  }
  
  // Come back to the output file dir
  pbpb_Unfo->cd();

  cout <<"MC = " << isMC <<endl;	
  
//  TH1D *hCent = new TH1D("hCent","",nbins_cent,boundaries_cent);
 

  // if you change the binning you have to change these, too
  // inclusive trig eff
  /*
    float trigEffInc[6]={0.777265,
			 0.95765,
			 0.998357,
			 0.999941,
			 1.,
			 1.};
  */

    // b-jet trig eff
    /*
    float trigEffbJet[6]={0.660276,
		      0.908997,
		      0.980793,
		      0.998767,
		      0.999442,
		      1.};
    */



 
		
  // Fill pPb MC
   double etashift = 0. ;
   if(coll=="PbP") etashift = -0.465 ;
   if(coll=="PPb") etashift = 0.465 ;   
  if (!useMatrixFromFile) {
    for (Long64_t jentry2=0; jentry2<dataPbPb->tJet->GetEntries();jentry2++) {
      dataPbPb->tJet->GetEntry(jentry2);
      
      // change when we switch to centrality binning
      int cBin = 0;
  
      if(! (dataPbPb->HLT_PAJet20_noJetID_v1) && ! (dataPbPb->HLT_PAJet40_noJetID_v1) && ! (dataPbPb->HLT_PAJet60_noJetID_v1) && ! (dataPbPb->HLT_PAJet80_noJetID_v1) && ! (dataPbPb->HLT_PAJet100_noJetID_v1)  ) continue ;   
      if ( dataPbPb->refpt  < 0. ) continue;
    //  if ( dataPbPb->refpt  < 20. ) continue;
   //   if ( dataPbPb->jtpt  < 20. ) continue;
      if(dataPbPb->subid!=0) continue ;
      if(fabs(dataPbPb->vz)>15.) continue ;
     if(dataPbPb->rawpt<22) continue ;
      //if(dataPbPb->chargedMax/dataPbPb->jtpt < 0.05) continue ;
  //   if(dataPbPb->jtpt>4*(dataPbPb->pthat)) continue;
   //  if(dataPbPb->rawpt<25) continue ;
 //    if((dataPbPb->chargedSum+dataPbPb->photonSum+dataPbPb->neutralSum+dataPbPb->muSum+dataPbPb->eSum)/dataPbPb->jtpt>1.01) continue ;
 //    if ( (dataPbPb->jteta+etashift)  > 2. || (dataPbPb->jteta+etashift) < -2.5 ) continue;
    if ( (dataPbPb->jteta+etashift)  > etamax || (dataPbPb->jteta+etashift) < etamin ) continue;
  //  if ( (dataPbPb->refeta+etashift)  > etamax || (dataPbPb->refeta+etashift) < etamin ) continue;
   //  if ( (dataPbPb->jteta)  > etamax || (dataPbPb->jteta) < etamin ) continue;
 //     if ( (dataPbPb->refeta+etashift)  > etamax || (dataPbPb->refeta+etashift) < etamin ) continue;
    //  if ( fabs(dataPbPb->jteta+0.465)  > 1. ) continue;
      if ( dataPbPb->refpt<0) dataPbPb->refpt=0;
      if (doBjets && fabs(dataPbPb->refparton_flavorForB)!=5) continue;
      if (doBjets && dataPbPb->jtptB < recoJetPtCut) continue;
      if (!doBjets && dataPbPb->jtpt < recoJetPtCut) continue;
  //    if(TMath::Abs(dataPbPb->refparton_flavor)>6 ) continue ;  //! quark jet martix
    //  if((dataPbPb->refparton_flavor)!=21 ) continue ; //gluon jet matrix
      //if (!doTrigEffCorr && dataPbPb->isTrig <1) continue;
     // if ( dataPbPb->isTrig <1) continue;
      
      //if(!doBjets)if(dataPbPb->refpt < 50 && dataPbPb->jtptA>120) continue;
      //if(doBjets)if(dataPbPb->refpt < 50 && dataPbPb->jtptB>120) continue;

      //for centrality selection using HF sum energy with |eta|>4
   //   if((dataPbPb->hiHFplusEta4+dataPbPb->hiHFminusEta4)<2.87) continue ;

      float weight = dataPbPb->weight;

      if(doTrigEffCorr){
	for(int iBin = 0; iBin<nbins_rec; iBin++){
	  float myJtPt = 0.;
	  if(doBjets) myJtPt = dataPbPb->jtptB;
	  else myJtPt = dataPbPb->jtpt;
	  if(myJtPt > boundaries_rec[iBin] && myJtPt < boundaries_rec[iBin+1]){
	    if(doBjets) weight/= trigEffbJet[iBin];
	    else weight/= trigEffInc[iBin];
	  }							  
	}
      }
     if(isMC) {
      if (jentry2 % 2 == 1) {
	if(doBjets){
                uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->jtptB,weight);
                uhist[cBin]-> hMeasJECSys->Fill(dataPbPb->jtptB*(1.+0.02/nbins_cent*(nbins_cent-0)),weight);
          } 
	else { 
         uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->jtpt,weight);
    //    uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->rawpt,weight); 
//	 uhist[cBin]-> hGen->Fill(dataPbPb->refpt,weight);   
     //   if((dataPbPb->refparton_flavor)==21) uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->rawpt,weight);
     //    if(TMath::Abs(dataPbPb->refparton_flavor)<=6 && TMath::Abs(dataPbPb->refparton_flavor)>0) uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->rawpt,weight);
           }
      }	  
      if (jentry2 % 2 == 0) {
  //    uhist[cBin]-> hGen->Fill(dataPbPb->refpt*(1-0.02),weight);
	uhist[cBin]-> hGen->Fill(dataPbPb->refpt,weight);   
     //    if((dataPbPb->jteta+etashift)  > etamin && (dataPbPb->jteta+etashift) < etamax) uhist[cBin]-> hMeas->Fill(dataPbPb->jtpt,weight);  	 
         uhist[cBin]-> hMeas->Fill(dataPbPb->jtpt,weight);  	 
     //     uhist[cBin]-> hMeas->Fill(dataPbPb->rawpt,weight);  	 
	//uhist[cBin]-> hMeasJECSys->Fill(dataPbPb->jtpt*(1.+0.02/nbins_cent*(nbins_cent-i)),weight); 
	// FIXME!!!!!!  i is supposed to be a loop over centrality !!!!
          uhist[cBin]-> hMeasJECSys->Fill(dataPbPb->jtpt*(1.+0.02/nbins_cent*(nbins_cent-0)),weight); 
      }
   }
  else {
   //   uhist[cBin]-> hGen->Fill(dataPbPb->refpt*(1-0.02),weight);
      uhist[cBin]-> hGen->Fill(dataPbPb->refpt,weight);
      uhist[cBin]-> hMeas->Fill(dataPbPb->jtpt,weight);
   //   uhist[cBin]-> hMeas->Fill(dataPbPb->rawpt,weight);
     uhist[cBin]-> hMeasJECSys->Fill(dataPbPb->jtpt*(1.+0.02/nbins_cent*(nbins_cent-0)),weight);
      if(DoSmear) 
     //   uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->jtpt*(1+(fsmear->Eval(dataPbPb->jtpt))*fgaus->GetRandom()),weight);
        uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,gRandom->Gaus(dataPbPb->jtpt, fsmear->Eval(dataPbPb->jtpt)),weight);
     else   uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->jtpt,weight);
   //  else   uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->jtpt*(1-0.02),weight);
  //   else   uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->rawpt,weight);
      } 

   } //PbPb entries loop 

    //pp will just fill the last index of the centrality array

  } //Fill matrix for pp and PbPb

  if (isMC==0) {
    // Use measured histogram from Matt & Kurt's file
	   
    // PbPb file:

    TFile *infMatt = new TFile(fileNamePbPb_data);
    TH1D *hMattPbPb= NULL;
    TH1D *hRebinPbPb ; // = NULL;
    TH1D *hTagEffPbPb = NULL;
    if(doBjets){
      hMattPbPb = (TH1D*) infMatt->Get("hRawBData");
      hTagEffPbPb = (TH1D*) infMatt->Get("hBEfficiencyMC");
    }
      else {
        if(TMath::Abs(etamin)==1. && TMath::Abs(etamax)==1. )
          hMattPbPb= (TH1D*) infMatt->Get("jetptEtaBin-10_10");
        else
          hMattPbPb= (TH1D*) infMatt->Get(Form("jetptEtaBin%.f_%.f", etamin*10, etamax*10));
 
    }
  
   //   hMattPbPb->Scale(pPbSigma*1.e3);
   //    hMattPbPb->Scale(1./(frac*Lumi*1.e6));
  //     hMattPbPb->Scale(1./(frac*Lumi*1.e9*pPbSigma));
//       hMattPbPb->Scale(1./nMB);
   //  hMattPbPb->Scale(1./(etamax-etamin));    
   //  hMattPbPb->Scale(1./2.);    
  //  else hMattPbPb = (TH1D*) infMatt->Get("hjtpt");
  //  divideBinWidth(hMattPbPb);
           
    // Need to match the binning carefully, please double check whenever you change the binning
     hRebinPbPb = (TH1D*)hMattPbPb->Rebin(nbins_rec, Form("DataJetPtEtaBin%.f_%.f",etamin*10, etamax*10),boundaries_rec);

       //fixed by Yaxian
/*       for(int ibin=1;ibin<=hRebinPbPb->GetNbinsX();ibin++){
          cout <<"data bins =" << hRebinPbPb->GetBinCenter(ibin) << " bin content =" << hRebinPbPb->GetBinContent(ibin)<<endl;
      }
*/       for(int ibin=1;ibin<=uhist[0]->hMeas->GetNbinsX();ibin++){
        uhist[0]->hMeas->SetBinContent(ibin,0);
        uhist[0]->hMeas->SetBinError(ibin,0); 
        float binCenter = uhist[0]->hMeas->GetBinCenter(ibin);
        float testBinWidth = uhist[0]->hMeas->GetBinWidth(ibin);
        float testLowEdge = uhist[0]->hMeas->GetBinLowEdge(ibin);
        float binContent =  hRebinPbPb->GetBinContent(hRebinPbPb->FindBin(binCenter));
         float binError = hRebinPbPb->GetBinError(hRebinPbPb->FindBin(binCenter));   

      //  cout <<"get bins " << testLowEdge<<" content =" <<binContent<<endl ;
        uhist[0]->hMeas->SetBinContent(ibin,binContent);  
        uhist[0]->hMeas->SetBinError(ibin,binError);  
       }
  } // real data histograms
  
    


   //=========================================Response Matrix========================================================= 

  cout <<"Response Matrix..."<<endl;
	
 // TCanvas * cMatrix = new TCanvas("cMatrix","cMatrix",800,400);
  TCanvas * cMatrix = new TCanvas("cMatrix","cMatrix",600,600);
//  cMatrix->Divide(nbins_cent,1);
 cMatrix->cd(1);
  TH2D * hResponse[nbins_cent+1];
  for (int i=0;i<nbins_cent;i++){
 //   cMatrix->cd(i+1);
    if(doParameterizedMatrix){
                  uhist[i]->hMatrix = (TH2D*)generateSmoothMatrix(i,uhist[i]->hMatrix,20.,recoJetPtCut,etamin,etamax, coll.Data(), 1, 1)->Clone(Form("hMatrix_%d",i));
/*    if(isMC){
       TH1D *hFunMeas =  (TH1D*) uhist[i]->hMatrix->ProjectionY(Form("hFunMeas_cent%d",i));
         for(int ibin=1;ibin<=uhist[i]->hMeas->GetNbinsX();ibin++){
             uhist[i]->hMeas->SetBinContent(ibin,0);
             uhist[i]->hMeas->SetBinError(ibin,0);
             float binCenter = uhist[i]->hMeas->GetBinCenter(ibin);
             float testBinWidth = uhist[i]->hMeas->GetBinWidth(ibin);
             float testLowEdge = uhist[i]->hMeas->GetBinLowEdge(ibin);
             float binContent =  hFunMeas->GetBinContent(hFunMeas->FindBin(binCenter));
             float binError = hFunMeas->GetBinError(hFunMeas->FindBin(binCenter));

                uhist[i]->hMeas->SetBinContent(ibin,binContent);
                 uhist[i]->hMeas->SetBinError(ibin,binError);
          } 
    }
*/ 
  }
    if (!useMatrixFromFile) {
/*       if(doParameterizedMatrix){
//      uhist[i]->hMatrix = (TH2D*)FillSmoothMatrix(i,uhist[i]->hMatrix,20.,recoJetPtCut,etamin,etamax, coll.Data(), 1, 1)->Clone(Form("hMatrix_%d",i));
      uhist[i]->hMatrix = (TH2D*)generateSmoothMatrix(i,uhist[i]->hMatrix,20.,recoJetPtCut,etamin,etamax, coll.Data(), 1, 1)->Clone(Form("hMatrix_%d",i));
 //     uhist[i]-> hGen = (TH1D*) uhist[i]->hMatrix->ProjectionX(Form("hGen_cent%d",i));
 //     uhist[i]-> hMeas = (TH1D*) uhist[i]->hMatrix->ProjectionY(Form("hMeas_cent%d",i));
    }
*/
      if(uhist[i]->hMatrix->GetEntries()<=0) cout << "Wrong!!!the matrix is empty!!!"<<endl ;
      uhist[i]->hMatrix->SetMaximum(1.e-3);
      uhist[i]->hMatrix->SetMinimum(1e-12);
      TF1 *f = new TF1("f","[0]*pow(x+[2],[1])");
      f->SetParameters(1e10,-8.8,40);
      for (int y=1;y<=uhist[i]->hMatrix->GetNbinsY();y++) {
	double sum=0;
	for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {
	  if (uhist[i]->hMatrix->GetBinContent(x,y)<=1*uhist[i]->hMatrix->GetBinError(x,y)) {
	    uhist[i]->hMatrix->SetBinContent(x,y,0);
	    uhist[i]->hMatrix->SetBinError(x,y,0);
	  }
	  sum+=uhist[i]->hMatrix->GetBinContent(x,y);
	}
				
	for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {	   
	  double ratio = 1;
	  uhist[i]->hMatrix->SetBinContent(x,y,uhist[i]->hMatrix->GetBinContent(x,y)*ratio);
	  uhist[i]->hMatrix->SetBinError(x,y,uhist[i]->hMatrix->GetBinError(x,y)*ratio);
	}
      }

    }

    uhist[i]->hResponse = (TH2D*)uhist[i]->hMatrix->Clone(Form("hResponse_cent%d",i));
    TH1D *hProj = (TH1D*)uhist[i]->hResponse->ProjectionY(Form("hProj_cent%d",i));


    for (int y=1;y<=uhist[i]->hResponse->GetNbinsY();y++) {
      double sum=0;
      for (int x=1;x<=uhist[i]->hResponse->GetNbinsX();x++) {
	if (uhist[i]->hResponse->GetBinContent(x,y)<=1*uhist[i]->hResponse->GetBinError(x,y)) {
	  uhist[i]->hResponse->SetBinContent(x,y,0);
	  uhist[i]->hResponse->SetBinError(x,y,0);
	}
	sum+=uhist[i]->hResponse->GetBinContent(x,y);
      }
			
     for (int x=1;x<=uhist[i]->hResponse->GetNbinsX();x++) {  	
	if (sum==0) continue;
	double ratio =1.;
       ratio = 1./sum ;
       if (hProj->GetBinContent(y)==0) ratio = 1e-100/sum;
         else ratio = hProj->GetBinContent(y)/sum;
//          double val = uhist[i]->hResponse->GetBinContent(y,x)/sum;
//         val=val*100;
//         val=ceil(val);
//         double nearest=val/100;
//          uhist[i]->hResponse->SetBinContent(y,x,nearest);
       uhist[i]->hResponse->SetBinContent(x,y,uhist[i]->hResponse->GetBinContent(x,y)*ratio);
        uhist[i]->hResponse->SetBinError(x,y,uhist[i]->hResponse->GetBinError(x,y)*ratio);
      }
    }
		
    
    uhist[i]->hResponseNorm = (TH2D*)uhist[i]->hMatrix->Clone(Form("hResponseNorm_cent%d",i));
    for (int x=1;x<=uhist[i]->hResponseNorm->GetNbinsX();x++) {
      double sum=0;
      for (int y=1;y<=uhist[i]->hResponseNorm->GetNbinsY();y++) {
	if (uhist[i]->hResponseNorm->GetBinContent(x,y)<=1*uhist[i]->hResponseNorm->GetBinError(x,y)) {
	  uhist[i]->hResponseNorm->SetBinContent(x,y,0);
	  uhist[i]->hResponseNorm->SetBinError(x,y,0);
	}
	sum+=uhist[i]->hResponseNorm->GetBinContent(x,y);
      }
			
      for (int y=1;y<=uhist[i]->hResponseNorm->GetNbinsY();y++) {  	
	if (sum==0) continue;
	double ratio = 1./sum;
	uhist[i]->hResponseNorm->SetBinContent(x,y,uhist[i]->hResponseNorm->GetBinContent(x,y)*ratio);
	uhist[i]->hResponseNorm->SetBinError(x,y,uhist[i]->hResponseNorm->GetBinError(x,y)*ratio);
      }
    }

    if (!useMatrixFromFile) uhist[i]->hMatrixFit = (TH2D*)uhist[i]->hMatrix->Clone(Form("hMatrixFit_cent%d",i));
    uhist[i]->hMatrixFit->SetName(Form("hMatrixFit_cent%d",i));

     uhist[i]->hResponse->SetTitleOffset(1.4, "Y");
    uhist[i]->hResponse->SetTitleOffset(1.2, "X");
  //  if(isMC)uhist[i]->hResponse->SetMinimum(1.e-8);
  //  else
      uhist[i]->hResponse->SetMinimum(1.e-8);
    uhist[i]->hResponse->DrawCopy("colz");
/*    uhist[i]->hMatrix->SetTitleOffset(1.4, "Y");
    uhist[i]->hMatrix->SetTitleOffset(1.2, "X");
    uhist[i]->hMatrix->Draw("colz");
		
   cMatrix->cd(i+1)->Modified();
   cMatrix->cd(i+1);
*/ 
 }

  // cMatrix->SetSelected(cMatrix);
//  cMatrix->Update();
	
  cout << "==================================== UNFOLD ===================================" << endl;
	
  //char chmet[100]; 
	
  // ======================= Reconstructed pp and PbPb spectra =========================================================
	
  TH1D *hRecoBW[nbins_cent+1], *hRecoBinByBinBW[nbins_cent+1], *hMeasBW[nbins_cent+1], *hGenBW[nbins_cent+1];
 TH1D *hReproduced[nbins_cent+1]; 

//added feb19
  TRandom3 *random = new TRandom3(0);
  TH1D *preUnfoldingSpec0 = new TH1D("histSpectraPre","", nbins_truth, boundaries_truth);
  TH1D *postUnfoldingSpec0 = new TH1D("histSpectraPost","", nbins_truth, boundaries_truth);
  //TH1D *unfoldedSpecOnePt = new TH1D("histOnePtBin","",nUnfoldingTrials, 0, nUnfoldingTrials);
  TH1D *postUnfoldingSpecOnePt[nbins_truth];
  double meanCS[nbins_truth] = {0.};
  double sigmaCS[nbins_truth] = {0.};	
  TH1D *preSpec0 = new TH1D("histPre","", nbins_truth, boundaries_truth);
 

  pbpb_Unfo->cd();

  for (int i=0;i<nbins_cent;i++) {
    // Do Bin-by-bin
 //   TH1D *hBinByBinCorRaw = (TH1D*)uhist[i]->hResponse->ProjectionY(); 
    TH1D *hBinByBinCorRaw = (TH1D*)uhist[i]->hMatrix->ProjectionY(); 
     hBinByBinCorRaw->Sumw2();
  //    if(isMC) uhist[i]->hMeas = (TH1D*) hBinByBinCorRaw->Clone(Form("hMeas_cent%d", i)); 
  //  TH1D *hMCGen           = (TH1D*)uhist[i]->hResponse->ProjectionX(); // gen
    TH1D *hMCGen           = (TH1D*)uhist[i]->hMatrix->ProjectionX(); // gen
    hMCGen->Sumw2(); 
 //   uhist[i]->hGen = (TH1D*) hMCGen->Clone(Form("hGen_cent%d", i)); 

    hBinByBinCorRaw->Divide(hMCGen);
    TF1 *f = new TF1("f","[0]+[1]*x");
 //   hBinByBinCorRaw->Fit("f","LL ","",recoJetPtCut,850);
    TH1D* hBinByBinCor = (TH1D*)hBinByBinCorRaw->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    uhist[i]->hRecoBinByBin = (TH1D*) uhist[i]->hMeas->Clone(Form("hRecoBinByBin_cent%d",i));
    uhist[i]->hRecoBinByBin->Divide(hBinByBinCor);
    TH1D *hPrior;//=(TH1D*) functionHist(fPow,uhist[i]->hMeas,Form("hPrior_cent%d",i));
    if(isMC) 
   //    hPrior=(TH1D*)hMCGen->Clone(Form("hPrior_cent%d", i));
       hPrior=(TH1D*)uhist[i]->hGen->Clone(Form("hPrior_cent%d", i));
    //   hPrior=(TH1D*)uhist[i]->hMeas->Clone(Form("hPrior_cent%d", i));
    else {
     // if(i==nbins_cent) 
     //    hPrior=(TH1D*)hMCGen->Clone(Form("hPrior_cent%d", i));
     //    hPrior=(TH1D*)uhist[i]->hGen->Clone(Form("hPrior_cent%d", i));
    // else {
         hPrior=(TH1D*)uhist[i]->hGen->Clone(Form("hPrior_cent%d", i));
    //     hPrior=(TH1D*)ftrue->GetHistogram();
      for(int ibin = 1 ; ibin <= uhist[i]->hGen->GetNbinsX(); ibin++){
           double binCenter = uhist[i]->hGen->GetBinCenter(ibin);
          hPrior->SetBinContent(ibin, ftrue->Eval(binCenter));
        }
  //    hPrior=(TH1D*)uhist[i]->hRecoBinByBin->Clone(Form("hPrior_cent%d", i));
 //        hPrior=(TH1D*)uhist[i]->hGen->Clone(Form("hPrior_cent%d", i));
  //     hPrior=(TH1D*)uhist[i]->hMeas->Clone(Form("hPrior_cent%d", i));
   //   }
     }
     removeZero(hPrior);
    for(int ibin=1 ; ibin < hPrior->GetNbinsX(); ibin++){
      cout <<" prior bin center =" <<hPrior->GetBinCenter(ibin) << " content =" << hPrior->GetBinContent(ibin) <<endl ;
  } 		
    // Do unfolding
    //if (isMC) uhist[i]->hMeas = (TH1D*)uhist[i]->hMatrix->ProjectionY()->Clone(Form("hMeas_cent%d",i));
    prior myPrior(uhist[i]->hMatrixFit,hPrior, 0);
  //  prior myPrior(uhist[i]->hMatrixFit,uhist[i]->hMeas,0);
  //  myPrior.unfold(uhist[i]->hMeas,1);
    myPrior.unfold(uhist[i]->hMeas,nBayesianIter);
//    hPrior = (TH1D*)uhist[i]->hGen->Clone("hPrior");
//    hPrior = (TH1D*)uhist[i]->hMeas->Clone(Form("hPrior_cent%d",i));
  //  TH1D *hReweighted = (TH1D*)(TH1D*)uhist[i]->hResponse->ProjectionY(Form("hReweighted_cent%d",i));
		
    bayesianUnfold myUnfoldingJECSys(uhist[i]->hMatrixFit,hPrior,0);
    myUnfoldingJECSys.unfold(uhist[i]->hMeasJECSys,nBayesianIter);
    bayesianUnfold myUnfoldingSmearSys(uhist[i]->hMatrixFit,hPrior,0);
    myUnfoldingSmearSys.unfold(uhist[i]->hMeasSmearSys,nBayesianIter);
    bayesianUnfold myUnfolding(uhist[i]->hMatrixFit,hPrior,0);
    myUnfolding.unfold(uhist[i]->hMeas,nBayesianIter);
    cout <<"Unfolding bin "<<i<<endl;
    
    delete hBinByBinCorRaw;
    delete hMCGen;


//added feb19
  preSpec0 = (TH1D*) uhist[i]->hMeas->Clone(Form("hMeasToBeUnfolded"));
  divideBinWidth(preSpec0);
  for(int j = 1; j <= uhist[i]->hMeas->GetNbinsX(); j++){
    meanCS[j-1] = uhist[i]->hMeas->GetBinContent(j);
    //sigmaCS[j-1] = uhist[i]->hMeas->GetBinError(j);
    sigmaCS[j-1] = uhist[i]->hGen->GetBinError(j);
  }
  for(int k = 1; k <= nUnfoldingTrials; k++){
    for(int j = 1; j <= nbins_truth; j++){
      preUnfoldingSpec0->SetBinContent(j, random->Gaus(meanCS[j-1], sigmaCS[j-1]));
      preUnfoldingSpec0->SetBinError(j, sigmaCS[j-1]/sqrt(nUnfoldingTrials));
      uhist[i]->hMeasPreUnfoldingSpec2D->SetBinContent(k, j, preUnfoldingSpec0->GetBinContent(j));
      uhist[i]->hMeasPreUnfoldingSpec2D->SetBinError(k, j, preUnfoldingSpec0->GetBinError(j));
    }
    bayesianUnfold myUnfoldingMulti(uhist[i]->hMatrixFit,hPrior,0);
    myUnfoldingMulti.unfold(preUnfoldingSpec0,nBayesianIter);    
    postUnfoldingSpec0 = (TH1D*) myUnfoldingMulti.hPrior->Clone(Form("Unfolded_Trial%d", k));
    divideBinWidth(postUnfoldingSpec0);
    for(int j = 1; j <= nbins_truth; j++){
      uhist[i]->hRecoUnfoldedSpec2D->SetBinContent(k, j, postUnfoldingSpec0->GetBinContent(j));
      uhist[i]->hRecoUnfoldedSpec2D->SetBinError(k, j, postUnfoldingSpec0->GetBinError(j));
    }
    delete postUnfoldingSpec0;
  }

  for(int j = 1; j <= nbins_truth; j++){
    //postUnfoldingSpecOnePt[j] = new TH1D(Form("postHistSpectra%d",j),"", 100, meanCS[j]-10.*sigmaCS[j], meanCS[j]+10.*sigmaCS[j]);
    postUnfoldingSpecOnePt[j] = new TH1D(Form("postHistSpectra%d",j),"", 100, 0, 1);
    for(int k = 1; k <= nUnfoldingTrials; k++){
      //unfoldedSpecOnePt->SetBinContent(i, uhist[i]->hRecoUnfoldedSpec2D->GetBinContent(i, j)); 
      postUnfoldingSpecOnePt[j]->Fill(uhist[i]->hRecoUnfoldedSpec2D->GetBinContent(k, j));
    }
    uhist[i]->hRecoUnfoldedSpec->SetBinContent(j, postUnfoldingSpecOnePt[j]->GetMean());
    uhist[i]->hRecoUnfoldedSpec->SetBinError(j, postUnfoldingSpecOnePt[j]->GetRMS());
    delete postUnfoldingSpecOnePt[j];
  }
 //end feb19
  

  // Iteration Systematics
    for (int j=2;j<=10;j++)
      {

	bayesianUnfold myUnfoldingSys(uhist[i]->hMatrixFit,myPrior.hPrior,0);
	myUnfoldingSys.unfold(uhist[i]->hMeas,j);
	uhist[i]->hRecoIterSys[j]  = (TH1D*) myUnfoldingSys.hPrior->Clone(Form("hRecoRAA_IterSys%d_cent%d",j,i));
      }

    uhist[i]->hReco         = (TH1D*) uhist[i]->hRecoIterSys[nBayesianIter]->Clone(Form("Unfolded_cent%i",i));
    uhist[i]->hRecoJECSys   = (TH1D*) myUnfoldingJECSys.hPrior->Clone(Form("UnfoldedJeCSys_cent%i",i));
    uhist[i]->hRecoSmearSys   = (TH1D*) myUnfoldingSmearSys.hPrior->Clone(Form("UnfoldedSmearSys_cent%i",i));
    uhist[i]->hRecoBinByBin->SetName(Form("UnfoldedBinByBin_cent%i",i));

//    cout <<"Unfolding bin 22222222"<<i<<endl;
    
    if (doToy) {
      TCanvas *cToy = new TCanvas("cToy","toy",600,600);
      cToy->cd();
      int nExp=1000;
      TH1D *hTmp[nbins_truth+1];
      TH1D *hTmp2[nbins_truth+1];
      for (int j=1;j<=nbins_truth;j++) {
	hTmp[j] = new TH1D(Form("hTmp%d",j),"",200,0,10.+uhist[i]->hReco->GetBinContent(j)*2);
	hTmp2[j] = new TH1D(Form("hTmp2%d",j),"",200,0,10.+uhist[i]->hRecoBinByBin->GetBinContent(j)*2);
      }
      for (int exp =0; exp<nExp; exp++) {
	TH1D *hToy = (TH1D*)uhist[i]->hMeas->Clone();   
	TH2D *hMatrixToy = (TH2D*)uhist[i]->hMatrixFit->Clone();
	hToy->SetName("hToy");
	if (exp%100==0) cout <<"Pseudo-experiment "<<exp<<endl;
	for (int j=1;j<=hToy->GetNbinsX();j++) {
	  double value = gRandom->Poisson(uhist[i]->hMeas->GetBinContent(j));
	  hToy->SetBinContent(j,value);
	}
				
	for (int j=1;j<=hMatrixToy->GetNbinsX();j++) {
	  for (int k=1;k<=hMatrixToy->GetNbinsY();k++) {
	    double value = gRandom->Gaus(uhist[i]->hMatrixFit->GetBinContent(j,k),uhist[i]->hMatrixFit->GetBinError(j,k));
	    hMatrixToy->SetBinContent(j,k,value);
	  }
	}

	prior myPriorToy(hMatrixToy,hToy,0.0);
	myPriorToy.unfold(hToy,1);
	bayesianUnfold myUnfoldingToy(hMatrixToy,myPriorToy.hPrior,0.0);
	myUnfoldingToy.unfold(hToy,nBayesianIter);
	TH1D *hRecoTmp = (TH1D*) myUnfoldingToy.hPrior->Clone();
				
	for (int j=1;j<=hRecoTmp->GetNbinsX();j++) {
	  hTmp[j]->Fill(hRecoTmp->GetBinContent(j));
	}
	delete hToy;
	delete hRecoTmp;
	delete hMatrixToy;
      }
      TF1 *fGaus = new TF1("fGaus","[0]*TMath::Gaus(x,[1],[2])");
      for (int j=1;j<=nbins_truth;j++)
	{

	  f->SetParameters(hTmp[j]->GetMaximum(),hTmp[j]->GetMean(),hTmp[j]->GetRMS());
				
	  if (hTmp[j]->GetMean()>0) {
	    hTmp[j]->Fit(fGaus,"LL Q ");
	    hTmp[j]->Fit(fGaus,"LL Q ");
	    uhist[i]->hReco->SetBinError(j,f->GetParameter(2));
	  }	       
	  f->SetParameters(hTmp2[j]->GetMaximum(),hTmp2[j]->GetMean(),hTmp2[j]->GetRMS());
	  if (hTmp2[j]->GetMean()>0) {
	    hTmp2[j]->Fit(fGaus,"LL Q ");
	    hTmp2[j]->Fit(fGaus,"LL Q ");
	    uhist[i]->hRecoBinByBin->SetBinError(j,f->GetParameter(2));
	  }	       
	  delete hTmp[j];
	  delete hTmp2[j];
	}
    }  //finish the doToy part

   //do binwidth normalization
/*   divideBinWidth(uhist[i]->hMeas);
   divideBinWidth(uhist[i]->hReco);
   divideBinWidth(uhist[i]->hRecoBinByBin);
   divideBinWidth(uhist[i]->hGen);
  */
   //   divideBinWidth(uhist[i]->hGen); 
    uhist[i]->hMeas->SetMarkerStyle(20);
    uhist[i]->hMeas->SetMarkerColor(2);
    uhist[i]->hReco->SetMarkerStyle(25);
    uhist[i]->hReco->SetName(Form("hReco_cent%d",i));
    
    uhist[i]->hReco->SetTitle("Baysian Unfolded");
    uhist[i]->hRecoBinByBin->SetTitle("Bin-by-bin Unfolded");
    uhist[i]->hReco->SetXTitle("p_{T} (GeV/c)");    
    uhist[i]->hReco->SetYTitle("Counts");    
    uhist[i]->hReco->GetXaxis()->SetNdivisions(505);
    uhist[i]->hReco ->GetYaxis()->SetTitleOffset(1.3);
    uhist[i]->hReco->GetXaxis()->SetTitleOffset(1.2);
    //uhist[i]->hReco->Draw("");    
    uhist[i]->hReco->SetAxisRange(recoJetPtCut,600,"X");
	    uhist[i]->hReco->SetAxisRange(recoJetPtCut,600, "X");
	    uhist[i]->hReco->GetXaxis()->SetRangeUser(recoJetPtCut,600);
	    uhist[i]->hReco->GetXaxis()->SetNdivisions(505);
	   // uhist[i]->hReco->DrawCopy("");   
	     
	    uhist[i]->hGen->SetLineWidth(2);
	    uhist[i]->hGen->SetLineColor(2);
	   // if(isMC)uhist[i]->hGen->Draw("hist same");
	   // uhist[i]->hReco->Draw("same");    
	    uhist[i]->hRecoBinByBin->SetMarkerStyle(28);
	    uhist[i]->hRecoBinByBin->SetMarkerColor(4);
	    uhist[i]->hRecoBinByBin->SetLineColor(4);
	  //  uhist[i]->hRecoBinByBin->DrawCopy("same");    

	    hReproduced[i] = (TH1D*)myUnfolding.hReproduced->Clone(Form("hReproduced_cent%d",i));
	    hReproduced[i]->SetMarkerColor(1);
	    hReproduced[i]->SetMarkerStyle(24);
	    //uhist[i]->hMeas->Draw("same");    

	    hRecoBW[i] = (TH1D*)uhist[i]->hReco->Clone(Form("hReco%d",i));
	    hRecoBinByBinBW[i] = (TH1D*)uhist[i]->hRecoBinByBin->Clone(Form("hRecoBinByBin%d",i));
	    hMeasBW[i] = (TH1D*)uhist[i]->hMeas->Clone(Form("hMeas%d",i));
	  //  if(isMC)
               hGenBW[i] = (TH1D*)uhist[i]->hGen->Clone(Form("hGen%d",i));

	    divideBinWidth(hRecoBW[i]);    
	  //  if(isMC)
                 divideBinWidth(hGenBW[i]); 
           //   if(i==0) 
           //      hGenBW[i]->Scale(1./(etamax-etamin));  
	    divideBinWidth(hRecoBinByBinBW[i]);    
	    divideBinWidth(hMeasBW[i]);    
	    divideBinWidth(hReproduced[i]);
	    hRecoBW[i]->SetTitle("Baysian Unfolded");
	    hRecoBinByBinBW[i]->SetTitle("Bin-by-bin Unfolded");
   }

  TCanvas * cPbPb = new TCanvas("cPbPb","Comparison",600,600);
//  cPbPb->Divide(nbins_cent,1);
  cPbPb->cd(1);
    for (int i=0;i<nbins_cent;i++) {
           cPbPb->cd(i+1);
           cPbPb->cd(i+1)->SetLogy();   
	    hMeasBW[i]->SetAxisRange(recoJetPtCut,600, "X");
	    hMeasBW[i]->GetXaxis()->SetRangeUser(recoJetPtCut,600);
	  //  hMeasBW[i]->GetYaxis()->SetRangeUser(1.e-10, 1.e-3);
            hMeasBW[i]->GetXaxis()->SetNdivisions(505);
	    hMeasBW[i]->DrawCopy("");
	  //  hRecoBW[i]->DrawCopy();
	  //  if(isMC)
               hGenBW[i]->DrawCopy("hist,same");
	    hRecoBinByBinBW[i]->DrawCopy("same");
	    hRecoBW[i]->DrawCopy("same");
	    


	    TLegend *leg = new TLegend(0.5,0.5,0.9,0.9);
	    leg->SetBorderSize(0);
	    leg->SetFillStyle(0);
	/*    if(i==0)leg->AddEntry(uhist[i]->hMeas,"PPb","");
	    if(i==1)leg->AddEntry(uhist[i]->hMeas,"PP",""); 
	   leg->AddEntry(uhist[i]->hMeas,"Measured","pl");
	    leg->AddEntry(uhist[i]->hReco,"Bayesian unfolded","pl");
	    leg->AddEntry(uhist[i]->hRecoBinByBin,"Bin-by-bin unfolded","pl");
	    if(isMC)leg->AddEntry(uhist[i]->hGen,"Generator level truth","l");
	 */   if(i==0)leg->AddEntry(hRecoBW[i],"PPb","");
	    if(i==1)leg->AddEntry(hRecoBW[i],"PP","");
	   leg->AddEntry(hMeasBW[i],"Measured","pl");
	    leg->AddEntry(hRecoBW[i],"Bayesian unfolded","pl");
	    leg->AddEntry(hRecoBinByBinBW[i],"Bin-by-bin unfolded","pl");
	   // if(isMC)
               leg->AddEntry(hGenBW[i],"Generator level truth","l");
	    leg->Draw("same");
       }
     	    cPbPb->Update();

    TCanvas * cPbPbMeas = new TCanvas("cPbPbMeas","Measurement",600,600);
  //  cPbPbMeas->Divide((nbins_cent),1);
    cPbPbMeas->cd(1);
         for (int i=0;i<nbins_cent;i++) {
            cPbPbMeas->cd(i+1);
	    cPbPbMeas->cd(i+1)->SetLogy();   
	  //  uhist[i]->hMeas->SetAxisRange(30,300,"X");
	    hMeasBW[i]->SetAxisRange(recoJetPtCut,600,"X");
	    hMeasBW[i]->GetXaxis()->SetRangeUser(recoJetPtCut,600);
	    hMeasBW[i]->GetXaxis()->SetNdivisions(505);
            // divideBinWidth(hMeasBW[i]);
	  //   divideBinWidth(hReproduced);
	  //  uhist[i]->hMeas->Draw();
	    hMeasBW[i]->DrawCopy();
	    hReproduced[i]->DrawCopy("same");

	    TLegend *leg2 = new TLegend(0.5,0.5,0.85,0.9);
	    leg2->SetBorderSize(0);
	    leg2->SetFillStyle(0);
	  // if(isMC) leg2->AddEntry(hReproduced,"PYTHIA+HIJING","");
	  //  else leg2->AddEntry(hReproduced,"pPb","");
	    if(i==0)leg2->AddEntry(hMeasBW[i],"PPb","");
	    if(i==1)leg2->AddEntry(hMeasBW[i],"PP","");
	    leg2->AddEntry(hMeasBW[i],"Measured","pl");
	    leg2->AddEntry(hReproduced[i],"Reproduced","pl");

	    leg2->Draw("same");
	  }	     
	    cPbPbMeas->Update();

	 // ======================= Unfolding closure in MC =========================================================
	 TCanvas * cRatio ;
	    TH1D * hReco[nbins_cent+1];
	    TH1D * hRecoBinByBin[nbins_cent+1];
	    TH1D * hMeas[nbins_cent+1];
	    TH1D * hGenScale[nbins_cent+1];
	    TLegend *leg[nbins_cent+1];
	  TLine *line = new TLine(recoJetPtCut,1,600,1);
		line->SetLineStyle(2);
		line->SetLineWidth(2);

	    for (int i=0;i<nbins_cent;i++) {
		hReco[i]          = (TH1D*)uhist[i]->hReco->Clone(Form("hRecoScaled_Cent%d", i));
		hRecoBinByBin[i]          = (TH1D*)uhist[i]->hRecoBinByBin->Clone(Form("hRecoBinByBinScaled_Cent%d", i));
		hMeas[i]          = (TH1D*)uhist[i]->hMeas->Clone(Form("hMeasScaled_Cent%d", i));
	//	if(isMC) 
                    hGenScale[i]          = (TH1D*)uhist[i]->hGen->Clone(Form("hGenScaled_Cent%d", i));
                divideBinWidth(hReco[i]);
               // if(isMC)
                    divideBinWidth(hGenScale[i]);
                divideBinWidth(hRecoBinByBin[i]);
                divideBinWidth(hMeas[i]);
	    }
	  if(isMC){
	//  cRatio = new TCanvas("cRatio","Ratio",1200,600);
	  cRatio = new TCanvas("cRatio","Ratio",600,600);
	 //   cRatio->Divide(nbins_cent,1);
	      for (int i=0;i<nbins_cent;i++) {
            /*    divideBinWidth(hReco[i]);
                divideBinWidth(hGen[i]);
                divideBinWidth(hRecoBinByBin[i]);
                divideBinWidth(hMeas[i]);
*/
		  hMeas[i]->Divide(hGenScale[i]);
		  hRecoBinByBin[i]->Divide(hGenScale[i]);
		  hReco[i]->Divide(hGenScale[i]);
		  cRatio->cd(i+1);

			//hRecoPP->SetAxisRange(90,300,"X");
			hReco[i]->SetAxisRange(0.5,1.5,"Y");
		//	hReco[i]->SetAxisRange(20.,600.,"X");
			hReco[i]->GetXaxis()->SetRangeUser(recoJetPtCut,600.);
			hReco[i]->GetXaxis()->SetNdivisions(505);
			hReco[i]->SetMarkerStyle(25);
			hReco[i] ->SetLineColor(1);
			hReco[i] ->SetMarkerColor(1);
			hMeas[i]->SetMarkerStyle(20);
			hMeas[i]->SetLineColor(2);
			hMeas[i]->SetMarkerColor(2);
			hRecoBinByBin[i]->SetMarkerStyle(28);
			hRecoBinByBin[i]->SetLineColor(4);
			hRecoBinByBin[i]->SetMarkerColor(4);

			makeHistTitle(hReco[i],"",Form("Jet p_{T} (GeV/c)"),Form("Reco / Truth"));
			hReco[i]->GetYaxis()->SetTitleOffset(1.3);
	                hReco[i]->GetXaxis()->SetNdivisions(505);
         		hReco[i]->GetXaxis()->SetTitleOffset(1.2);
			hReco[i]->DrawCopy("");
			hRecoBinByBin[i]->Draw("same");
			hMeas[i]->DrawCopy("same");
			line->Draw("same");
			leg[i] = myLegend(0.52,0.65,0.85,0.9);
			 if(i==0)leg[i]->AddEntry(hReco[i],"PYTHIA+HIJING","");
		//	 if(i==1)leg[i]->AddEntry(hReco[i],"PYTHIA","");
			leg[i]->AddEntry(hReco[i],"Bayesian","pl");
			leg[i]->AddEntry(hRecoBinByBin[i],"Bin-by-bin","pl");
			leg[i]->AddEntry(hMeas[i],"no unfolding","pl");
			leg[i]->Draw("same");
			putCMSPrel(0.2,0.83,0.04);
			drawText("Anti-k_{T} Particle Flow Jets   R = 0.3",0.2,0.23,20);
			drawText("CMS Simulation",0.5,0.4,20);
			drawText(Form("%.1f< #eta_{CM} <%.1f ",etamin, etamax),0.5,0.31,20);
	      }
	  }
	  else {
	     cRatio = new TCanvas("cRatio","Ratio",600,600);
	      cRatio->cd(1);

	    
	      for (int i=0;i<nbins_cent;i++) {
/*            hMeas[i]            ->Scale(1./208);
                  hRecoBinByBin[i]            ->Scale(1./208);
                  hReco[i]            ->Scale(1./208);
 
                  hMeas[i]->Divide(hGenScale[nbins_cent]);
                  hRecoBinByBin[i]->Divide(hGenScale[nbins_cent]);
                  hReco[i]->Divide(hGenScale[nbins_cent]);
*/
                  hRecoBinByBin[i]->Divide(hMeas[i]);
                  hReco[i]->Divide(hMeas[i]);

		  hReco[i]->SetAxisRange(0,3.,"Y");
		  hReco[i]->SetAxisRange(recoJetPtCut,600,"X");
		  hReco[i]->GetXaxis()->SetRangeUser(recoJetPtCut,600);
		  hReco[i]->GetXaxis()->SetNdivisions(505);
		  hReco[i]->SetMarkerStyle(25);
		  hReco[i] ->SetLineColor(1);
		  hReco[i] ->SetMarkerColor(1);
		  hMeas[i]->SetMarkerStyle(20);
		  hMeas[i]->SetLineColor(2);
		  hMeas[i]->SetMarkerColor(2);
		  hRecoBinByBin[i]->SetMarkerStyle(28);
		  hRecoBinByBin[i]->SetLineColor(4);
		  hRecoBinByBin[i]->SetMarkerColor(4);
		 // if(i==0){
		//  makeHistTitle(hMeas[i],"",Form("Jet p_{T} (GeV/c)"),Form("Jet R_{pA}"));
		  makeHistTitle(hReco[i],"",Form("Jet p_{T} (GeV/c)"),Form("Unfolded/Measured"));
		  hReco[i]->GetYaxis()->SetTitleOffset(1.2);
		  hReco[i]->GetXaxis()->SetTitleOffset(1.2);
		  hReco[i]->Draw("");
		      leg[i] = myLegend(0.62,0.75,0.85,0.9);
		      if(i==0)leg[i]->AddEntry(hReco[i],"PPb","");
		      if(i==1)leg[i]->AddEntry(hReco[i],"PP",""); 
		      leg[i]->AddEntry(hReco[i],"Bayesian","pl");
		      leg[i]->AddEntry(hRecoBinByBin[i],"Bin-by-bin","pl");
		   //   leg[i]->AddEntry(hMeas[i],"no unfolding","pl");
		      leg[i]->Draw("same");
		//  hReco[i]->Draw("same");
		  hRecoBinByBin[i]->Draw("same");
		//  hReco[i]->Draw("same");
	      //    }
	       }
		 line->Draw("same");

		  putCMSPrel(0.2,0.83,0.04);
		 // drawText(Form("#intL dt = %.f #mub^{-1}",150.),0.2,0.78,20);
		  drawText("Anti-k_{T} PF Jets R = 0.3",0.2,0.33,20);
		  drawText(Form("%.1f< #eta_{CM} <%.1f ",etamin, etamax),0.2,0.27,20);

	  }
	    cRatio->Update();
	 
	  pbpb_Unfo->Write();

	  SysData systematics;
	  //TLine *line = new TLine(60,1,250,1);

	  // Iteration systematics
	//  TCanvas *cIterSys = new TCanvas("cIterSys","cIterSys",1200,600);
	  TCanvas *cIterSys = new TCanvas("cIterSys","cIterSys",600,600);
	//  cIterSys->Divide(2,1);
	/*  cIterSys->cd(2);
	  TH1D *hRecoIterSysPP[100];
	  TH1D *hRebinPP_tmp         = rebin(uhist[nbins_cent]->hReco, (char*)"hRebinPP_tmp");
	  TLegend *legBayesianIterPP = myLegend(0.4,0.7,0.9,0.9);
	  legBayesianIterPP->AddEntry("","PP","");
		 
	  for (int j=2;j<7;j++) {
	    hRecoIterSysPP[j] = rebin(uhist[nbins_cent]->hRecoIterSys[j],Form("hRecoIterSysPP_IterSys%d",j));
	    hRecoIterSysPP[j]->SetLineColor(colorCode[j-2]);
	    hRecoIterSysPP[j]->SetMarkerColor(colorCode[j-2]);
	    hRecoIterSysPP[j]->SetMarkerStyle(20);
	    hRecoIterSysPP[j]->Divide(hRebinPP_tmp);
	    if (j==2){
	      makeHistTitle(hRecoIterSysPP[j],(char*)"",(char*)"Jet p_{T} (GeV/c)",(char*)"Ratio (Unfolded / Nominal)");
	      hRecoIterSysPP[j]->SetTitleOffset(1.3,"Y");
	      hRecoIterSysPP[j]->SetTitleOffset(1.2,"X");
	     // if(isMC)
               hRecoIterSysPP[j]->SetAxisRange(0.5,1.5,"Y");
	    //  else hRecoIterSysPP[j]->SetAxisRange(0.,2.,"Y");
	     // hRecoIterSysPP[j]->SetAxisRange(20,500,"X");
	      hRecoIterSysPP[j]->GetXaxis()->SetRangeUser(recoJetPtCut,600);
	      hRecoIterSysPP[j]->GetXaxis()->SetNdivisions(505);
	      hRecoIterSysPP[j]->DrawCopy(); 
	    } else {
	      hRecoIterSysPP[j]->DrawCopy("same");
	    }
		 
	    checkMaximumSys(systematics.hSysIter[nbins_cent],hRecoIterSysPP[j],0,1.);
	    legBayesianIterPP->AddEntry(hRecoIterSysPP[j],Form("Iteration %d",j),"pl");     
	  }

	  legBayesianIterPP->Draw("same");
	  line->Draw("same");
	  drawEnvelope(systematics.hSysIter[nbins_cent],(char*)"hist same");

*/
	  cIterSys->cd(1);
	  TH1D *hRecoIterSysPbPb[100];
	  TH1D *hRebinPbPb_tmp         = rebin(uhist[0]->hReco, (char*)"hRebinPbPb_tmp");
	  TLegend *legBayesianIterPbPb = myLegend(0.4,0.7,0.9,0.9);
	  legBayesianIterPbPb->AddEntry("","PPb","");
	  for (int j=2;j<7;j++) {
	    hRecoIterSysPbPb[j] = rebin(uhist[0]->hRecoIterSys[j],Form("hRecoIterSysPbPb_IterSys%d",j));
	    hRecoIterSysPbPb[j]->SetLineColor(colorCode[j-2]);
	    hRecoIterSysPbPb[j]->SetMarkerColor(colorCode[j-2]);
	    hRecoIterSysPbPb[j]->SetMarkerStyle(20);
	    hRecoIterSysPbPb[j]->Divide(hRebinPbPb_tmp);
	    if (j==2){
	      makeHistTitle(hRecoIterSysPbPb[j],(char*)"",(char*)"Jet p_{T} (GeV/c)",(char*)"Ratio (Unfolded / Nominal)");
	      hRecoIterSysPbPb[j]->SetTitleOffset(1.3,"Y");
	      hRecoIterSysPbPb[j]->SetTitleOffset(1.2,"X");
	    //  if(isMC)
                 hRecoIterSysPbPb[j]->SetAxisRange(0.5,1.5,"Y");
	    //  else hRecoIterSysPbPb[j]->SetAxisRange(0.,2.,"Y");
	 //     hRecoIterSysPbPb[j]->SetAxisRange(20,500,"X");
	      hRecoIterSysPbPb[j]->GetXaxis()->SetRangeUser(recoJetPtCut,600);
	      hRecoIterSysPbPb[j]->GetXaxis()->SetNdivisions(505);
	      hRecoIterSysPbPb[j]->DrawCopy(); 
	    } else {
	      hRecoIterSysPbPb[j]->DrawCopy("same");
	    }
		 
	    checkMaximumSys(systematics.hSysIter[0],hRecoIterSysPbPb[j],0,1.);
	    legBayesianIterPbPb->AddEntry(hRecoIterSysPbPb[j],Form("Iteration %d",j),"pl");     
	  }
	  legBayesianIterPbPb->Draw("same");
	  line->Draw("same");
	  drawEnvelope(systematics.hSysIter[0],(char*)"hist same");

	 cIterSys->Update();

	  TString data ;
	  if(isMC) data="MC";
	  else data="Data";
	 TString anaType ;
	  if(doBjets) anaType="Bjet";
	  else anaType="Inclusive";

	  if(SavePlot){
	    cMatrix->SaveAs(Form("plots/%s%s%sResponseMatrixEtaBin%.f_%.f.gif", data.Data(), anaType.Data(), algoName[algo], etamin*10,etamax*10));
	    cPbPb->SaveAs(Form("plots/%s%s%sRawPtMatrixJetSpectraEtaBin%.f_%.f.gif",  data.Data(), anaType.Data(), algoName[algo],etamin*10,etamax*10));
	    cRatio->SaveAs(Form("plots/%s%s%sRawPtMatrixJetRatioEtaBin%.f_%.f.gif",  data.Data(), anaType.Data(), algoName[algo],etamin*10,etamax*10));
	    cIterSys->SaveAs(Form("plots/%s%s%sRawPtMatrixIterationSysEtaBin%.f_%.f.gif",  data.Data(), anaType.Data(), algoName[algo],etamin*10,etamax*10));
	   }

	}





