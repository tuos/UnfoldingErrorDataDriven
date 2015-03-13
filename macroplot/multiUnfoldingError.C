#include "TH1D.h"
#include "TH2D.h"
#include <TF1.h>
#include <TMath.h>
#include "TRandom3.h"
#include "TFile.h"
using namespace std;
void multiUnfoldingError(){

  static const double boundaries_truth[] = {3, 4, 5, 7, 9, 12, 15, 18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,1000};
  static const int nbins_truth = sizeof(boundaries_truth)/sizeof(Double_t)-1;
  const int nUnfoldingTrials = 1000;

  TFile *infile = new TFile("/home/tuos/workingpA/v1/Unfold/outputs/PPb_UnfoPriorGen_akPu3PFOfficialMCNoIDCut_MCJECv8_jtpt20_RecoPtScaleUpEtaBin-10_10_Inc_MultiUnfolding_v6.root","read");

TH1D* histSpectraUnfolded;
histSpectraUnfolded = (TH1D*)infile->Get(Form("histSpectraUnfolded0"));
//histSpectraUnfolded->Draw();

  TCanvas *c1 = new TCanvas("c1","c1",1,1,550,460);
  c1->SetFillColor(10);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderSize(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.17);
  c1->SetBottomMargin(0.15);
  c1->SetTopMargin(0.03);
  c1->SetRightMargin(0.03);
  //c1->Divide(2,1,0,0);
  gStyle->SetOptStat(0);
  c1->SetTicks(-1);
  //c1->SetLogx();

 TH1D* hist = new TH1D("hist","",200,49.,500);
 hist->SetXTitle("p_{T} [GeV/c]");
 hist->SetYTitle("Percent Error (Stat.) of Jet Cross Section");
 hist->SetMinimum(0.00001);
 hist->SetMaximum(0.03);
 hist->GetXaxis()->CenterTitle(1);
 hist->GetYaxis()->CenterTitle(1);
 hist->GetYaxis()->SetTitleOffset(1.6);
 hist->GetXaxis()->SetTitleOffset(1.2);
 hist->GetXaxis()->SetTitleSize(0.056);
 hist->GetYaxis()->SetTitleSize(0.05);
 hist->GetXaxis()->SetLabelSize(0.05);
 hist->GetYaxis()->SetLabelSize(0.05);
 hist->Draw();

  TH1D *ratio1 = new TH1D("hist1","", nbins_truth, boundaries_truth);
  TH1D *ratio2 = new TH1D("hist2","", nbins_truth, boundaries_truth);
  TH1D *ratio3 = new TH1D("hist3","", nbins_truth, boundaries_truth);
  
  TH1D *hMeas = (TH1D*)infile->Get("hMeas0");
  TH1D *hReco = (TH1D*)infile->Get("hReco0");
  for(int j = 1; j <= nbins_truth; j++){
   if(hMeas->GetBinContent(j)>0){
    ratio1->SetBinContent(j, hMeas->GetBinError(j)/hMeas->GetBinContent(j));
    ratio2->SetBinContent(j, hReco->GetBinError(j)/hReco->GetBinContent(j));
    ratio3->SetBinContent(j, histSpectraUnfolded->GetBinError(j)/histSpectraUnfolded->GetBinContent(j));
    //cout<<"j="<<j<<" pt="<<hMeas->GetBinCenter(j)<<", err="<<hMeas->GetBinError(j)/hMeas->GetBinContent(j)<<endl;
   }
  }

ratio1->SetMarkerStyle(24);
ratio2->SetMarkerStyle(25);
ratio3->SetMarkerStyle(28);
ratio1->SetMarkerColor(1);
ratio2->SetMarkerColor(4);
ratio3->SetMarkerColor(2);

  ratio1->SetLineColor(1);
  ratio2->SetLineColor(4);
  ratio3->SetLineColor(2);

ratio1->Draw("psame");
ratio2->Draw("psame");
ratio3->Draw("psame");

    TLegend *leg0 = new TLegend(0.25,0.665,0.44,0.895);
    leg0->SetFillColor(10);
    leg0->SetBorderSize(0.035);
    leg0->SetTextFont(42);
    leg0->SetTextSize(0.048);
    leg0->AddEntry(ratio1,"Error before Unfolding","p");
    leg0->AddEntry(ratio2,"Default Error after Unfolding","p");
    leg0->AddEntry(ratio3,"Data Driven Error after Unfolding","p");
    leg0->Draw();

    TLatex *tex1= new TLatex(92.01,0.0165,"Anti-k_{T} Particle Flow Jets: R=0.3, |#eta_{CM}|<1.0");
    tex1->SetTextColor(1);
    tex1->SetTextSize(0.045);
    tex1->SetTextFont(42);
    tex1->Draw();
    TLatex *tex1= new TLatex(252.01,0.0135,"Number of trials = 1000");
    tex1->SetTextColor(1);
    tex1->SetTextSize(0.045);
    tex1->SetTextFont(42);
    tex1->Draw();

c1->SaveAs("percentErrorDataDrivenN1000.pdf");

}


