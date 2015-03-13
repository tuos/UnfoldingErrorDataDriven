#ifndef __generate_smooth_matrix__H
#define __generate_smooth_matrix__H

#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TLine.h"

double _epsilon = 1e-6; // speed up calculations with acceptable loss of precision
const int nk = 4; // number of kernel parameters (excluding pt, eta)
const int Nevents = 1.e5 ;

double resolution[100];

TF1 *_kernel = 0;
int maxEntry=5;
int iFit=0;
const double pi=acos(-1.);
const double pi2=2*pi -1;
const char *fopt="R+";
const int knpx=5000;
float fitmin=0.00;
float fitmax=2.00;
const Double_t jetpt[]={32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,1000};
const int nJetPtBin = sizeof(jetpt)/sizeof(Double_t)-1 ;

TString coll;

TH1 * normalizeByBinWidth(TH1 *histo) {
    for(int i = 1; i <= histo->GetNbinsX(); i++) {
        float content = histo->GetBinContent(i);
        float error = histo->GetBinError(i);
        histo->SetBinContent(i,content/histo->GetBinWidth(i));
        histo->SetBinError(i,error/histo->GetBinWidth(i));
    }
    return histo ;
}

TH2D *generateSmoothMatrix(int iter, TH2D *hmatrix, double genPtmin=20., double recoPtmin=20., double etalo=-3., double etah=3., TString collision="PPb", bool binBbin_or_bayes=1, int _debug=1){

  TH1D * hgen = (TH1D*)hmatrix->ProjectionX(Form("hMatrixProjectionGen"));
//    hgen->Scale(1./(etah-etalo));
//    normalizeByBinWidth(hgen);
  TH1D * hpt = (TH1D*)hmatrix->ProjectionY(Form("hMatrixProjectionReco"));
//    hpt->Scale(1./(etah-etalo));
    normalizeByBinWidth(hpt);
    // fit on the generator level ture distributions
    TF1 * ftrue = new TF1("ftrue", "[0]*exp([1]/x)*pow(x,[2])*pow(1-x*cosh([4])/4000.,[3])", 10, 1000) ;
    ftrue->ReleaseParameter(0);
    ftrue->ReleaseParameter(1);
    ftrue->ReleaseParameter(2);
    ftrue->ReleaseParameter(3);
//    ftrue->SetParLimits(0,1.e5,2.e11);     //! A
//    ftrue->SetParLimits(1,-200,0.);  //! exp
//    ftrue->SetParLimits(2,-10, 0); //! power
//    ftrue->SetParLimits(3,0, 20); //! eta constrain power
//    ftrue->SetParLimits(4,-3., 3.); //! eta
    ftrue->SetNpx(8000);

    // initial fit of the NLO curve to a histogram
    if(collision=="PPb"){
        if(etalo==-1.5 && etah ==-1.0) ftrue->SetParameters(1.30469e+11,  -51.459,   -5.77699,  7.11111,   -1.25);
        if(etalo==-1.0 && etah ==-0.5) ftrue->SetParameters(1.16979e+12 , -81.477 ,  -6.18545 , 4.05446,   -0.75);
        if(etalo==-0.5  && etah ==0.5) ftrue->SetParameters(2.71508e+05 , -3.38264e+01 ,  -5.34095e+00 , 1.24268e+01,   -5.00000e-01);
        if(etalo==0.5  && etah ==1.0) ftrue->SetParameters(7.70751e+11 , -75.3034 ,  -6.10602 , 4.62163,   0.75);
        if(etalo==1.0  && etah ==1.5) ftrue->SetParameters(1.14359e+11 , -52.6724 ,  -5.73266 , 7.42525,   1.25);
        if(etalo==1.5  && etah ==2.0) ftrue->SetParameters(1.03152e+03 , 7.69961e+00 ,  -4.32889e+00 , 1.83222e+01,   1.50000e+00);
        if(etalo==2.0  && etah ==2.5) ftrue->SetParameters(6.67874e-03 , 1.03588e+02 ,  -1.63869e+00 , 2.93514e+01,   2.00000e+00);
        if(etalo==-1.0  && etah ==1.0) ftrue->SetParameters(5.87637e+07 , -9.27923e+01 ,  -6.28103e+00 , 2.58283e+00,   0.0);
//        if(etalo==-1.5 && etah ==-1.0) ftrue->SetParameters(9096.18 , 8.01576 ,  -4.6605 , 13.1916,   1.25);
//        if(etalo==-1.0 && etah ==-0.5) ftrue->SetParameters(97959 , -20.3457 ,  -5.12859 , 12.563,   0.75);
//        if(etalo==-0.5  && etah ==0.5) ftrue->SetParameters(252631 , -32.82 ,  -5.30642 , 13.3043,   0);
//        if(etalo==0.5  && etah ==1.0) ftrue->SetParameters(97420.4 , -20.258 ,  -5.12755 , 12.5692,   -0.75);
//        if(etalo==1.0  && etah ==1.5) ftrue->SetParameters(9073.59 , 8.04636 ,  -4.66001 , 13.194,   -1.25);
//        if(etalo==1.5  && etah ==2.0) ftrue->SetParameters(402.079 , 39.2633 ,  -4.01103 , 13.8887,   -1.75);
//        if(etalo==2.0  && etah ==2.5) ftrue->SetParameters(125.387 , 35.0107 ,  -3.68867 , 13.4849,   -2.25);
//        if(etalo==-1.0  && etah ==1.0) ftrue->SetParameters(252631 , -32.82 ,  -5.30642 , 13.3043,   0);
    }
    if(collision=="PbP"){
        if(etalo==-2.5 && etah ==-2.0) ftrue->SetParameters(1.23531e+03 , 2.19167e+00 ,  -4.39614e+00 , 8.76154e+00,   -2.50000e+00);
        if(etalo==-2.0 && etah ==-1.5) ftrue->SetParameters(1.70343e+04 , -2.03677e+01 ,  -4.95367e+00 , 8.30510e+00,   -2.00000e+00);
        if(etalo==-1.5  && etah ==-1.0) ftrue->SetParameters(1.16317e+04 , -9.18673e+00 ,  -4.88392e+00 , 9.80253e+00,   -1.50000e+00);
        if(etalo==-1.0  && etah ==-0.5) ftrue->SetParameters(6.94513e+04 , -2.71453e+01 ,  -5.23892e+00, 1.01414e+01,   -1.00000e+00);
        if(etalo==-0.5  && etah ==0.5) ftrue->SetParameters(2.62864e+05 , -3.46184e+01 ,  -5.34838e+00 , 1.23353e+01,   -5.00000e-01);
        if(etalo==0.5  && etah ==1.0) ftrue->SetParameters(8.31729e+04 , -3.29895e+01 ,  -5.25295e+00 , 1.50652e+01,   5.00000e-01);
        if(etalo==1.0  && etah ==1.5) ftrue->SetParameters(1.95006e+03 , 9.93962e+00 ,  -4.51668e+00 , 1.90073e+01,   1.00000e+00);
        if(etalo==-1.0  && etah ==1.0) ftrue->SetParameters(5.79426e+05 , -3.60963e+01 ,  -5.37458e+00 , 8.94894e+00,   -1.00000e+00);
    }
//    if(collision=="PbP"){
//        if(etalo==-2.5 && etah ==-2.0) ftrue->SetParameters(125.387 , 35.0107 ,  -3.68867 , 13.4849,   -2.25);
//        if(etalo==-2.0 && etah ==-1.5) ftrue->SetParameters(402.079 , 39.2633 ,  -4.01103 , 13.8887,   -1.75);
//        if(etalo==-1.5  && etah ==-1.0) ftrue->SetParameters(9073.59 , 8.04636 ,  -4.66001 , 13.194,   -1.25);
//        if(etalo==-1.0  && etah ==-0.5) ftrue->SetParameters(97420.4 , -20.258 ,  -5.12755 , 12.5692,   -0.75);
//        if(etalo==-0.5  && etah ==0.5) ftrue->SetParameters(252631 , -32.82 ,  -5.30642 , 13.3043,   0);
//        if(etalo==0.5  && etah ==1.0) ftrue->SetParameters(97959 , -20.3457 ,  -5.12859 , 12.563,   0.75);
//        if(etalo==1.0  && etah ==1.5) ftrue->SetParameters(9096.18 , 8.01576 ,  -4.6605 , 13.1916,   1.25);
//        if(etalo==-1.0  && etah ==1.0) ftrue->SetParameters(252631 , -32.82 ,  -5.30642 , 13.3043,   0);
//    }

    if(hgen->Integral()>0){
//        ftrue->ReleaseParameter(0);
//        ftrue->ReleaseParameter(1);
//        ftrue->ReleaseParameter(2);
//        ftrue->ReleaseParameter(3);
//        ftrue->SetParLimits(0,5.e6,2.e12);     //! A
//        ftrue->SetParLimits(1,-200,0.);  //! exp
//        ftrue->SetParLimits(2,-10, 0); //! power
//        ftrue->SetParLimits(3,0, 20); //! eta constrain power
//        ftrue->SetParLimits(4,-3., 3.); //! eta
//        ftrue->SetNpx(2000);
//        ftrue->SetParameter(0,hgen->GetEntries());
//        ftrue->SetParameter(1,-18);
//        ftrue->SetParameter(2,-5.2);
//        ftrue->SetParameter(3,8.9);
////        ftrue->FixParameter(4,(etalo+etah)/2.);
//        ftrue->FixParameter(4,etalo);
//        hgen->Fit("ftrue",fopt,"",30,800);

        if(etalo==2.0)
            ftrue->SetParameter(0,hgen->Integral()*2);
            else
                ftrue->SetParameter(0,hgen->GetEntries()*2);
 //
  //      ftrue->FixParameter(4,(etalo+etah)/2.);
        ftrue->FixParameter(4,etalo);
        if(TMath::Abs(etalo)>1.5)
        hgen->Fit("ftrue","R","NO",30,500);
        else
            hgen->Fit("ftrue","R","NO",30,800);
        double FitSys[5];
        for(int i = 0 ; i < 5 ; i++ ){
            FitSys[i] = 0;
        //    FitSys[i] = hgen->GetFunction("ftrue")->GetParError(i);
            ftrue->SetParameter(i,hgen->GetFunction("ftrue")->GetParameter(i)+FitSys[i]);

        }
//        ftrue->SetParameter(0,hgen->GetFunction("ftrue")->GetParameter(0));
//        ftrue->SetParameter(1,hgen->GetFunction("ftrue")->GetParameter(1));
//        ftrue->SetParameter(2,hgen->GetFunction("ftrue")->GetParameter(2));
//        ftrue->SetParameter(3,hgen->GetFunction("ftrue")->GetParameter(3));
//        ftrue->SetParameter(4,hgen->GetFunction("ftrue")->GetParameter(4));
//        
        ftrue->SetParError(0,hgen->GetFunction("ftrue")->GetParError(0));
        ftrue->SetParError(1,hgen->GetFunction("ftrue")->GetParError(1));
        ftrue->SetParError(2,hgen->GetFunction("ftrue")->GetParError(2));
        ftrue->SetParError(3,hgen->GetFunction("ftrue")->GetParError(3));
        ftrue->SetParError(4,hgen->GetFunction("ftrue")->GetParError(4));
        
        ftrue->SetChisquare(hgen->GetFunction("ftrue")->GetChisquare());
        ftrue->SetNDF(hgen->GetFunction("ftrue")->GetNDF());
        cout << "Chi2 =" << hgen->GetFunction("ftrue")->GetChisquare() << " /NDF =" << hgen->GetFunction("ftrue")->GetNDF() <<endl ;
    }

    TF1 *fres = new TF1("fres", "sqrt([0]*[0]+pow(([1]/sqrt(x)),2)+pow([2]/x, 2))", 10, 1000);
    fres->SetParNames("C", "S", "N");
    fres->ReleaseParameter(0);
    fres->ReleaseParameter(1);
    fres->ReleaseParameter(2);
//    fres->SetParLimits(0,0.,0.09);     //! s0
//    fres->SetParLimits(1,-1.5,1.5);  //! S1
//    fres->SetParLimits(2,-1.,5.); //! S2
    fres->SetNpx(8000);

    if(coll=="PPb"){
        if(etalo==-1.5 && etah ==-1.0) fres->SetParameters(0.0609197,  0.906949 ,  2.48503);
        if(etalo==-1.0 && etah ==-0.5) fres->SetParameters(0.0714619,  0.985635 ,  1.74154e-05);
        if(etalo==-0.5 && etah ==0.5) fres->SetParameters(3.33388e-02,  1.03408e+00 ,  -1.23498e-04);
        if(etalo==0.5 && etah ==1.0) fres->SetParameters(0.061536,  0.861556 ,  -1.85362e-06);
        if(etalo==1.0 && etah ==1.5) fres->SetParameters(0.0639121,  0.877887 ,  0.000157601);
        if(etalo==1.5 && etah ==2.0) fres->SetParameters(7.46450e-02,  8.55510e-01 ,  -2.53529e+00);
        if(etalo==2.0 && etah ==2.5) fres->SetParameters(8.62815e-02,  4.45425e-01 ,  4.07510e+00);
        if(etalo==-1.0 && etah ==1.0) fres->SetParameters(0.0631308,  1.01648 ,  -2.21389e-05);
    }
    if(coll=="PbP"){
        if(etalo==-2.5 && etah ==-2.0)  fres->SetParameters(7.13885e-02,  5.11680e-01 ,  4.27961e+00);
        if(etalo==-2.0 && etah ==-1.5)  fres->SetParameters(6.87628e-02,  8.92007e-01 ,  2.40102e+00);
        if(etalo==-1.5  && etah ==-1.0) fres->SetParameters(5.64397e-02,  8.80112e-01 ,  7.96925e-06);
        if(etalo==-1.0  && etah ==-0.5)  fres->SetParameters(6.78418e-02,  7.94595e-01 ,  3.53689e-04);
        if(etalo==-0.5  && etah ==0.5)  fres->SetParameters(5.39973e-02,  8.92564e-01 ,  6.89275e-05);
        if(etalo==0.5  && etah ==1.0) fres->SetParameters(7.81749e-02,  7.87037e-01 ,  3.54545e+00);
        if(etalo==1.0  && etah ==1.5) fres->SetParameters( 8.84454e-02,  3.83879e-01 ,  5.19901e+00);
        if(etalo==-1.0  && etah ==1.0)  fres->SetParameters(6.17235e-02,  8.87176e-01 ,  -1.48551e-06);
    }

    TF1 * fjes = new TF1("fjes", "[0]+[1]/x+[2]/pow(x,2)", 10, 1000) ;
    fjes->SetParNames("Mean", "p1", "p2");
    fjes->SetNpx(2000);

    if(coll=="PPb"){
        if(etalo==-1.5 && etah ==-1.0) fjes->SetParameters(0.997719,  -1.04062,   58.3865);
        if(etalo==-1.0 && etah ==-0.5) fjes->SetParameters( 0.993326,  -0.302366,   43.0027);
        if(etalo==-0.5 && etah ==0.5) fjes->SetParameters(1.00851e+00,  -3.04442e+00,   1.24786e+02);
        if(etalo==0.5 && etah ==1.0) fjes->SetParameters(0.998289,  -0.547141,   29.7031);
        if(etalo==1.0 && etah ==1.5) fjes->SetParameters(0.999787,  -0.504515,   33.0348);
        if(etalo==1.5 && etah ==2.0) fjes->SetParameters(1.00024e+00,  -1.19350e+00,   7.03205e+01);
        if(etalo==2.0 && etah ==2.5) fjes->SetParameters(9.92744e-01,  -5.09833e-01,   4.50533e+01);
        if(etalo==-1.0 && etah ==1.0) fjes->SetParameters(0.999928,  -1.62684,   60.6889);
    }
    if(coll=="PbP"){
        if(etalo==-2.5 && etah ==-2.0)  fjes->SetParameters(1.05954e+00,  -9.73692e+00,   3.19695e+02);
        if(etalo==-2.0 && etah ==-1.5)  fjes->SetParameters(1.00127e+00,  -9.62021e-01,   4.77740e+01);
        if(etalo==-1.5  && etah ==-1.0) fjes->SetParameters(1.00251e+00,  -4.18485e-01,   2.92556e+01);
        if(etalo==-1.0  && etah ==-0.5)  fjes->SetParameters(1.00742e+00,  -1.59507e+00,   6.30174e+01);
        if(etalo==-0.5  && etah ==0.5)  fjes->SetParameters(1.00958e+00,  -2.47912e+00,   9.51166e+01);
        if(etalo==0.5  && etah ==1.0) fjes->SetParameters(1.00682e+00,  -2.19965e+00,   1.04539e+02);
        if(etalo==1.0  && etah ==1.5) fjes->SetParameters(1.01220e+00,  -4.16544e+00,   1.70781e+02);
        if(etalo==-1.0  && etah ==1.0)  fjes->SetParameters(1.00918e+00,  -2.42773e+00,   9.93835e+01);
    }
    
//    TF1 *f1 = new TF1("f1","(([0]/(2*pi*[1]*[1]))*exp(-1.*((x-[2])*(x-[2])/(2*[1]*[1])))*exp([3]/x)*pow(x,[4])*pow(1-x*cosh([6])/4000.,[5]))",30,800);
//    TF1 *fgauss = new TF1("fgauss","([0]/(2*pi*[1]*[1]))*exp(-1.*((x-[2])*(x-[2])/(2*[1]*[1])))",10,1000);
    TF1 *f1 = new TF1("f1","[0]*exp([1]/x)*pow(x,[2])*pow(1-x*cosh([4])/3500.,[3])",10,1000);
    const int nbinspt = hmatrix->GetNbinsX();
    TH1D * recoptgenbin[nbinspt];
    double ptresolution[nbinspt];
    double ptscale[nbinspt];
    TH1D *hMean, *hSigma;
    hMean = new TH1D("hMean", "hMean", nJetPtBin, jetpt);
    hMean->Sumw2();
    hSigma = new TH1D("hSigma", "hSigma", nJetPtBin, jetpt);
    hSigma->Sumw2();
    
    TF1  * fitReco[nbinspt];
//    f1->ReleaseParameter(0);
//    f1->ReleaseParameter(1);
//    f1->ReleaseParameter(2);
//    f1->ReleaseParameter(3);
//    f1->ReleaseParameter(4);
//    if(hpt->Integral()>0){
//        double ptreco = hpt->GetMean();
//        double ptreco1 = hpt->GetBinCenter(hpt->GetMinimumBin());
//        double ptreco2 = hpt->GetBinCenter(hpt->GetMaximumBin());
//        double pt ;
//        if(ptreco1<ptreco2) pt = fjes->GetX(ptreco,ptreco1,ptreco2);
//        else
//            pt = fjes->GetX(ptreco,ptreco2,ptreco1);
//        double res = fres->Eval(pt)*pt;
//        
//        f1->SetParameter(0,hpt->GetEntries()*2);
//        f1->SetParameter(1,ftrue->GetParameter(1));
//        f1->SetParameter(2,ftrue->GetParameter(2));
//        f1->SetParameter(3,ftrue->GetParameter(3));
//        f1->FixParameter(4,(etalo+etah)/2.);
//    //    f1->FixParameter(4,etalo);
//        hpt->Fit("f1","R","NO",30,800);
//        f1->SetParameter(0,hpt->GetFunction("f1")->GetParameter(0));
//        f1->SetParameter(1,hpt->GetFunction("f1")->GetParameter(1));
//        f1->SetParameter(2,hpt->GetFunction("f1")->GetParameter(2));
//        f1->SetParameter(3,hpt->GetFunction("f1")->GetParameter(3));
//        f1->SetParameter(4,hpt->GetFunction("f1")->GetParameter(4));
//        f1->SetParError(0,hpt->GetFunction("f1")->GetParError(0));
//        f1->SetParError(1,hpt->GetFunction("f1")->GetParError(1));
//        f1->SetParError(2,hpt->GetFunction("f1")->GetParError(2));
//        f1->SetParError(3,hpt->GetFunction("f1")->GetParError(3));
//        f1->SetParError(4,hpt->GetFunction("f1")->GetParError(4));
//        
//        f1->SetChisquare(hpt->GetFunction("f1")->GetChisquare());
//        f1->SetNDF(hpt->GetFunction("f1")->GetNDF());
//        
//        cout << "Chi2 =" << hpt->GetFunction("f1")->GetChisquare() << " /NDF =" << hpt->GetFunction("f1")->GetNDF() <<endl ;
//    }

    for(int ipt = 0 ; ipt < nJetPtBin ; ipt++){
        recoptgenbin[ipt]  = new TH1D(Form("RecoPtGenPtRatioBin%.f_%.f", jetpt[ipt], jetpt[ipt+1]),Form("RecoPtBin%.f_%.f", jetpt[ipt], jetpt[ipt+1]), 500, 0., 5.);
        recoptgenbin[ipt]->Sumw2();
        double lbin = (double)hMean->GetXaxis()->GetBinLowEdge(ipt);
        double hbin = (double)hMean->GetXaxis()->GetBinLowEdge(ipt+1);
        
        for(int ip=hmatrix->GetXaxis()->FindBin(jetpt[ipt]);ip<hmatrix->GetXaxis()->FindBin(jetpt[ipt+1]);++ip){
            double binCenter = hmatrix->GetXaxis()->GetBinCenter(ip);
            for(int ibin =1 ; ibin < hmatrix->GetNbinsY()+1;++ibin){
                double binCenterY = hmatrix->GetYaxis()->GetBinCenter(ibin);
                double binContent = hmatrix->GetBinContent(ip, ibin);
                recoptgenbin[ipt]->Fill(binCenterY/binCenter, binContent);
            }
        }
    }
    for(int ipt = 0 ; ipt < nJetPtBin ; ipt++){
        if(recoptgenbin[ipt]->Integral()>0){
            hMean ->SetBinContent(ipt,recoptgenbin[ipt]->GetMean());
            hMean ->SetBinError(ipt,recoptgenbin[ipt]->GetMeanError());
            hSigma->SetBinContent(ipt,recoptgenbin[ipt]->GetRMS());
            hSigma ->SetBinError(ipt,recoptgenbin[ipt]->GetRMSError());
        }
    }
    if(hMean->Integral()>0){
        if(TMath::Abs(etalo)>1.5)
            hMean->Fit("fjes","", "No",30,500);
        else
            hMean->Fit("fjes","", "No",30,800);

    //    hMean->Fit("fjes","", "No", 30., 800.);
        fjes->SetParameter(0,hMean->GetFunction("fjes")->GetParameter(0));
        fjes->SetParameter(1,hMean->GetFunction("fjes")->GetParameter(1));
        fjes->SetParameter(2,hMean->GetFunction("fjes")->GetParameter(2));
        fjes->SetParError(0,hMean->GetFunction("fjes")->GetParError(0));
        fjes->SetParError(1,hMean->GetFunction("fjes")->GetParError(1));
        fjes->SetParError(2,hMean->GetFunction("fjes")->GetParError(2));
        cout << " JES Chi2 =" << hMean->GetFunction("fjes")->GetChisquare() << " /NDF =" << hMean->GetFunction("fjes")->GetNDF() <<endl ;

    }
    if(hSigma->Integral()>0){
        if(TMath::Abs(etalo)>1.5)
            hSigma->Fit("fres","", "No",30,500);
        else
            hSigma->Fit("fres","", "No",30,800);

        fres->SetParameter(0,hSigma->GetFunction("fres")->GetParameter(0));
        fres->SetParameter(1,hSigma->GetFunction("fres")->GetParameter(1));
        fres->SetParameter(2,hSigma->GetFunction("fres")->GetParameter(2));
        fres->SetParError(0,hSigma->GetFunction("fres")->GetParError(0));
        fres->SetParError(1,hSigma->GetFunction("fres")->GetParError(1));
        fres->SetParError(2,hSigma->GetFunction("fres")->GetParError(2));
        cout << " JER Chi2 =" << hSigma->GetFunction("fres")->GetChisquare() << " /NDF =" << hSigma->GetFunction("fres")->GetNDF() <<endl ;
        
    }

  // Create smeared theory curve
  double maxpt = 3950./cosh((etalo+etah)/2.);

 if (_debug)
    cout << "Calculate forward smearing and unfold hpt" << endl << flush;

    // Calculate smearing matrix
  if (_debug)
    cout << "Generating smearing matrix T..." << endl << flush;
  double tmp_eps = _epsilon;

  // NB: GetArray only works if custom x binning
  //outdir->cd();

  // Deduce range and binning for true and measured spectra
  vector<double> vx; // true
  vector<double> vy; // measured
  for (int i = 1; i != hgen->GetNbinsX()+1; ++i) {

    double x = hgen->GetBinCenter(i);
    double x1 = hgen->GetBinLowEdge(i);
    double x2 = hgen->GetBinLowEdge(i+1);
    double y = hgen->GetBinContent(i);

     if (x>=genPtmin && y>0) {
      if (vx.size()==0)vx.push_back(x1);
       vx.push_back(x2);
    }
  } //for true pt binning
  for (int i = 1; i != hpt->GetNbinsX()+1; ++i) {
    double x = hpt->GetBinCenter(i);
    double x1 = hpt->GetBinLowEdge(i);
    double x2 = hpt->GetBinLowEdge(i+1);
    double y = hpt->GetBinContent(i);

    if (x1>=recoPtmin && y>0) {
      if (vy.size()==0) vy.push_back(x1);
       vy.push_back(x2);
    }
  } // for measured pt binning
  
  // copy over relevant part of hpt
  TH1D *hreco = new TH1D(Form("hreco_jet_%d",iter),";p_{T,reco} (GeV)",
vy.size()-1,&vy[0]);
  for (int i = 1; i != hreco->GetNbinsX()+1; ++i) {
    int j = hpt->FindBin(hreco->GetBinCenter(i));
    double dpt = hpt->GetBinWidth(j);
    hreco->SetBinContent(i, hpt->GetBinContent(j)/dpt);
    hreco->SetBinError(i, hpt->GetBinError(j)/dpt);
  }
//  for (int i = 1; i != hgen->GetNbinsX()+1; ++i) {
//    int j = hpt->FindBin(hgen->GetBinCenter(i));
//    double dpt = hpt->GetBinWidth(j);
//    hgen->SetBinContent(i, hgen->GetBinContent(j)/dpt);
//    hgen->SetBinError(i, hgen->GetBinError(j)/dpt);
//  }

//    TH2D *mt = new TH2D(Form("mt_jet_%d",iter),"mt;p_{T,gen};p_{T,reco}", 200, 0., 800., 200, 0., 800.);
    TH2D *mt = new TH2D(Form("mt_jet_%d",iter),"mt;p_{T,gen};p_{T,reco}", nbins_truth, boundaries_truth,nbins_rec, boundaries_rec);
  TH1D *mx = new TH1D(Form("mx_jet_%d",iter),"mx;p_{T,gen};#sigma/dp_{T}",
nbins_truth, boundaries_truth);
  TH1D *my = new TH1D(Form("my_jet_%d",iter),"my;p_{T,reco};#sigma/dp_{T}",
                      nbins_rec, boundaries_rec); //vy.size()-1, &vy[0]);

  // From http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html
  // For 1-dimensional true and measured distribution bins Tj and Mi,
  // the response matrix element Rij gives the fraction of events
  // from bin Tj that end up measured in bin Mi.
    cout << "x bins = " << mt->GetNbinsX()<<endl ;
  for (int i = mt->GetXaxis()->FindBin(genPtmin); i < mt->GetNbinsX()+1; i++) {

    double binlo =  mt->GetXaxis()->GetBinLowEdge(i);
    double binhi = mt->GetXaxis()->GetBinWidth(i)+binlo;
//      double xsec = ftrue->Integral(binlo, binhi)/(binhi-binlo);
//     double xsec = ftrue->Integral(binlo, binhi);
      for(int ievt = 0 ; ievt < Nevents ; ievt++){
          double genpt = gRandom->Uniform(binlo,binhi);
 //         double genpt = ftrue->GetRandom(binlo,binhi);
           double cs = ftrue->Eval(genpt);
          int ibin = mt->GetXaxis()->FindBin(genpt);
//          double jes = fjes->Eval(genpt);
          double jes ;
//          if(genpt < 60.) jes = hMean->GetBinContent(ibin);
//          else
              jes = fjes->Eval(genpt);
          
          double resol = gRandom->Gaus(jes,fres->Eval(genpt*jes));
//          double resol = gRandom->Gaus(1.0,fres->Eval(genpt));
          double recpt = genpt*resol;
     //     int jbin = mt->GetYaxis()->FindBin(recpt);
          double xsec = ftrue->Integral(mt->GetXaxis()->GetBinLowEdge(ibin), mt->GetXaxis()->GetBinLowEdge(ibin+1))/(mt->GetXaxis()->GetBinLowEdge(ibin+1)-mt->GetXaxis()->GetBinLowEdge(ibin));
          
          // cout<<genpt<<" "<<recpt<<" "<<xsec<<endl;
//          mx->Fill(genpt,xsec/Nevents);
//          my->Fill(recpt,xsec/Nevents);
//          mt->Fill(genpt,recpt,xsec/Nevents);
          mx->Fill(genpt,cs/Nevents);
          my->Fill(recpt,cs/Nevents);
          mt->Fill(genpt,recpt,cs/Nevents);
      }
    // double res = fres->Integral(ptgen1, ptgen2)/(ptgen2-ptgen1)*ptgen;

//      fgauss->SetParameter(0, hpt->Integral());
//      fgauss->SetParameter(1, res);
//      fgauss->SetParameter(2, ptgen);
     // cout << fgauss->Integral()<<endl;
//    if (_debug )
//        cout << "ptgen1: "<< ptgen1 << " ptgen2: "<< ptgen2 << " ygen: "<< ygen << " ptgen: "<< ptgen << endl;
//     for (int j = 1; j != mt->GetNbinsY()+1; ++j) {
//         double ptreco1 = mt->GetYaxis()->GetBinLowEdge(j);
//         double ptreco2 = mt->GetYaxis()->GetBinLowEdge(j+1);
//         double ptreco = mt->GetYaxis()->GetBinCenter(j);
//         double smearf = TMath::Gaus(ptreco, ptgen, res, kTRUE); ///TMath::Gaus(ptreco, ptgen, res, kTRUE)->Integral();
////         double smearf ;
////         if(ptreco < ptgen) smearf = fgauss->Integral(0, ptreco)/fgauss->Integral(0,-1);
////         else smearf = fgauss->Integral(ptreco, -1)/fgauss->Integral(0,-1);
//  //       cout << "smear factors =" << smearf <<endl ;
// //           mt->SetBinContent(i, j, (ftrue->Eval(ptgen)*(ptreco2-ptreco1)*(ptgen2-ptgen1)), smearf);
//         mt->SetBinContent(i, j, (ftrue->Eval(ptgen)*smearf*(ptreco2-ptreco1)*(ptgen2-ptgen1)));
// //        mt->SetBinContent(i, j, (f1->Eval(ptreco)*smearf*(ptreco2-ptreco1)));
// //        mt->SetBinContent(i, j, (hgen->GetBinContent(i)*smearf*(ptreco2-ptreco1)*(ptgen2-ptgen1)));
//         if (_debug){
//             cout << "ptgen: "<< ptgen << " ptreco: "<< ptreco << " ygen: "<< ygen << " yreco: "<< f1->Integral(ptreco1, ptreco2)/(ptreco2-ptreco1) << " smeared: "<< ftrue->Eval(ptgen)*smearf*(ptreco2-ptreco1) << endl;
//             cout << "print original now!!!" <<endl ;
//             cout << "ptgen: "<< hmatrix->GetXaxis()->GetBinCenter(i) << " ptreco: "<< hmatrix->GetYaxis()->GetBinCenter(j) << " matrix z = : "<< hmatrix->GetBinContent(i,j)  << endl;
//         }
//
//          //  mt->SetBinContent(i, j, f1->Eval(ptreco)*(ptreco2-ptreco1));
//       //  if(fitReco[i]->Eval(ptreco)/ftrue->Eval(ptgen))
//    //     cout << "ptgen =" << ptgen << "ptreco =" << ptreco << " y =" << f1->Eval(ptreco)*(ptreco2-ptreco1) <<endl;
//    //     }
//      } // for j
  } // for i

//    TF1 *fgauss = new TF1("fgauss","([0]/(2*pi*[1]*[1]))*exp(-1.*((x-[2])*(x-[2])/(2*[1]*[1])))",10,1000);
//
//    for(int i = 0 ; i < Nevents ; i++){
//        double ptgen = ftrue->GetRandom(genPtmin, 1000);
////        double ptgen = ftrue->GetRandom();
//        double cs = ftrue->Eval(ptgen);
//        double res = fres->Eval(ptgen)*ptgen;
//        double jes = fjes->Eval(ptgen)*ptgen;
////        fgauss->SetParameter(0, 1.);
////        fgauss->SetParameter(1, res);
////        fgauss->SetParameter(2, jes);
////        double ptreco = fgauss->GetRandom(10., 900);
//
// //       double ptreco = gRandom->Gaus(ptgen, res);
//        double ptreco = gRandom->Gaus(jes, res);
////        double prab = TMath::Gaus(ptreco,ptgen, res, kTRUE);
////        if(ptreco<recoPtmin) continue ;
// //       mt->Fill(ptgen, ptreco, cs*prab);
//        mt->Fill(ptgen, ptreco);
//    }
    cout << "filling matrix using functions done!!" <<endl ;

//    for(int i = 0 ; i < 10 ; i++){
//        cout << "i =" << i <<endl ;
//        double ptgen = ftrue->GetRandom(genPtmin, 600);
//        cout << "ptgen =" << ptgen <<endl ;
//        mx->Fill(ptgen);
//        double res = fres->Eval(ptgen)*ptgen;
//        cout << "res =" << res <<endl ;
//        double ptreco = gRandom->Gaus(ptgen, res);
//        cout << "ptreco =" << ptreco <<endl ;
//        my->Fill(ptreco);
//    }

//  for (int j = 1; j != mt->GetNbinsX()+1; ++j) {
//
//    double ptgen1 = min(4000./cosh((etalo+etah)/2.), mt->GetXaxis()->GetBinLowEdge(j));
//    double ptgen2 = min(4000./cosh((etalo+etah)/2.), mt->GetXaxis()->GetBinLowEdge(j+1));
//    double ygen = ftrue->Integral(ptgen1, ptgen2);
//    mx->SetBinContent(j, ygen);
//  }
//  
//  for (int i = 1; i != mt->GetNbinsY()+1; ++i) {
//      double ptreco1 = mt->GetYaxis()->GetBinLowEdge(i);
//      double ptreco2 = mt->GetYaxis()->GetBinLowEdge(i+1);
//      double dpt = mt->GetYaxis()->GetBinWidth(i);
////      double yreco = fnlos->Integral(ptreco1, ptreco2);
//      
//    double yreco(0);
//    for (int j = 1; j != mt->GetNbinsX()+1; ++j) {
//      yreco += mt->GetBinContent(j, i);
//    }
//
//    my->SetBinContent(i, yreco);
//  } // for i
  
  TH2D *mtu = (TH2D*)mt->Clone(Form("mtu_jet_%d",iter));
  for (int i = 1; i != mt->GetNbinsX()+1; ++i) {
    for (int j = 1; j != mt->GetNbinsY()+1; ++j) {
      if (mx->GetBinContent(i)!=0) {
//          mtu->SetBinContent(i, j, mt->GetBinContent(i,j) / mx->GetBinContent(j));
        mtu->SetBinContent(i, j, mt->GetBinContent(i,j));
      }
    } // for j
  } // for i


  // For BinByBin and SVD, need square matrix
  TH2D *mts(0);
  TH1D *mxs(0);


  mts = new TH2D(Form("mts_jet_%d",iter),"mts;p_{T,reco};p_{T,gen}",
vy.size()-1, &vy[0], vy.size()-1, &vy[0]);
  mxs = new TH1D(Form("mxs_jet_%d",iter),"mxs;p_{T,gen};#sigma/dp_{T}",
vy.size()-1, &vy[0]);
  for (int i = 1; i != mts->GetNbinsX()+1; ++i) {
    for (int j = 1; j != mts->GetNbinsY()+1; ++j) {
      double x = mts->GetBinCenter(i);
double y = mts->GetBinCenter(j);
int i2 = mt->GetXaxis()->FindBin(x);
int j2 = mt->GetYaxis()->FindBin(y);
mts->SetBinContent(i, j, mt->GetBinContent(i2, j2));
mts->SetBinError(i, j, mt->GetBinError(i2, j2));
      }
    }
    for (int i = 1; i != mxs->GetNbinsX()+1; ++i) {
      double x = mxs->GetBinCenter(i);
      int i2 = mx->FindBin(x);
      mxs->SetBinContent(i, mx->GetBinContent(i2));
      mxs->SetBinError(i, mx->GetBinError(i2));
    }


//RooUnfoldResponse *uResp;
//if(binBbin_or_bayes) uResp = new RooUnfoldResponse(my, mx, mt);
// else uResp = uResp = new RooUnfoldResponse(my, mxs, mts);

_epsilon = tmp_eps;

if(binBbin_or_bayes) return mtu;
 else return mts;

    if (_debug)
        cout << "done." << endl << flush;
    
}

#endif
