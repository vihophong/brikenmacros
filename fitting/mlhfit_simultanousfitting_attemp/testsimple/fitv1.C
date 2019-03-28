/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
///  'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #312
///
///  Performing fits in multiple (disjoint) ranges in one or more dimensions
///
/// \macro_output
/// \macro_code
/// \author 07/2008 - Wouter Verkerke


#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooCategory.h"
#include "TRandom3.h"
#include "TFile.h"
#include "fitF.cxx"
#include "fitFp1n.cxx"
#include "fitFp2n.cxx"
#include "RooGaussian.h"
#include "RooConstVar.h"


#include "RooRandom.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooAddition.h"
#include "RooDataHist.h"
#include "RooPoisson.h"
#include "RooPlot.h"

#include "RooNLLVar.h"

#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "RooStats/ProofConfig.h"
#include "RooStats/ToyMCSampler.h"
#include "RooProfileLL.h"


using namespace RooFit ;
using namespace RooStats;

Double_t p_deadtime=0.;
Double_t p_timerange=10.0;
Double_t nsigma=5.;


void getbackground(Double_t* par, char* infile)
{
    TFile *f=TFile::Open(infile);
    TTree* tree;
    f->GetObject("treeb",tree);
    Double_t nall=tree->Draw("x",Form("x<%f&&x>%f",-p_deadtime,-p_timerange),"goff");
    par[0]=nall;
    par[1]=nall/nsigma;
    par[2]=nall+nall*nsigma;


    TTree* treep1n;
    f->GetObject("treep1n",treep1n);
    nall=treep1n->Draw("x",Form("x<%f&&x>%f",-p_deadtime,-p_timerange),"goff");
    par[3]=nall;
    par[4]=nall/nsigma;
    par[5]=nall+nall*nsigma;

    TTree* treep2n;
    f->GetObject("treep2n",treep2n);
    nall=treep2n->Draw("x",Form("x<%f&&x>%f",-p_deadtime,-p_timerange),"goff");
    par[6]=nall;
    par[7]=nall/nsigma;
    par[8]=nall+nall*nsigma;

    f->Close();
}

void simfit(char* infile)
{
    //! define variables
    RooRealVar x("x","x",p_deadtime,p_timerange) ;
    RooRealVar* p[3];//hl.p1n,p2n

    Double_t hl=log(2)/0.12;
    Double_t hlplus=log(2)/0.01;
    Double_t hlminus=log(2)/1.;

    p[0]=new RooRealVar(Form("p%d",0),Form("p%d",0),hl,hlminus,hlplus);
    p[1]=new RooRealVar(Form("p%d",1),Form("p%d",1),0.2,0.,1.);
    p[2]=new RooRealVar(Form("p%d",2),Form("p%d",2),0.25,0.,1.);
    // sig pdf
    fitF totdecaymodel("totdecaymodel","totdecaymodel",x,*p[0]);
    fitFp1n decay1nmodel("decay1nmodel","decay1nmodel",x,*p[0],*p[1]);
    fitFp2n decay2nmodel("decay2nmodel","decay2nmodel",x,*p[0],*p[1],*p[2]);
    // bkg pdf
    RooPolynomial bkgtotdecaymodel("bkgtotdecaymodel","bkgtotdecaymodel",x);
    RooPolynomial bkgdecay1nmodel("bkgdecay1nmodel","bkgdecay1nmodel",x);
    RooPolynomial bkgdecay2nmodel("bkgdecay2nmodel","bkgdecay2nmodel",x);

    // number of singal for extened pdf

    Double_t bkgpar[9];
    getbackground(bkgpar,infile);
    RooRealVar nsigtotdecay("nsigtotdecay","number of nsignal events",15000,1000,150000);
    RooRealVar nbkgtotdecay("nbkgtotdecay","number of bkgs events",bkgpar[0],bkgpar[1],bkgpar[2]);

    RooRealVar nbkgdecay1n("nbkgdecay1n","number of bkgs events",bkgpar[3],bkgpar[4],bkgpar[5]);
    RooRealVar nbkgdecay2n("nbkgdecay2n","number of bkgs events",bkgpar[6],bkgpar[7],bkgpar[8]);



    // set constant background from backword timing
    //nbkgtotdecay.setConstant(kTRUE);
    //nbkgdecay1n.setConstant(kTRUE);
    //nbkgdecay2n.setConstant(kTRUE);


    // add pdf
    RooAddPdf final_totdecaymodel("final_totdecaymodel","final_totdecaymodel",RooArgList(totdecaymodel,bkgtotdecaymodel),RooArgList(nsigtotdecay,nbkgtotdecay));

    RooAddPdf final_decay1nmodel("final_decay1nmodel","final_decay1nmodel",RooArgList(decay1nmodel,bkgdecay1nmodel),RooArgList(nsigtotdecay,nbkgdecay1n));

    RooAddPdf final_decay2nmodel("final_decay2nmodel","final_decay2nmodel",RooArgList(decay2nmodel,bkgdecay2nmodel),RooArgList(nsigtotdecay,nbkgdecay2n));

    // Import tree to total decay data
    TFile *f=TFile::Open(infile);
    TTree* tree;
    f->GetObject("treeb",tree);
    RooDataSet* data=new RooDataSet("data","data",x,Import(*tree)) ;
    data->Print() ;


    // Import tree to p1n data
    TTree* treep1n;
    f->GetObject("treep1n",treep1n);
    RooDataSet* datap1n=new RooDataSet("data","data",x,Import(*treep1n)) ;
    datap1n->Print() ;

    // Import tree to p2n data
    TTree* treep2n;
    f->GetObject("treep2n",treep2n);
    RooDataSet* datap2n=new RooDataSet("data","data",x,Import(*treep2n)) ;
    datap2n->Print() ;



    nsigtotdecay.setVal(7.55537e+03);
    nbkgdecay1n.setVal(8.55097e+04);
    nsigtotdecay.setConstant();
    nbkgdecay1n.setConstant();


    // perform fit individually
    RooFitResult* fitres;
    //fitres=final_totdecaymodel.fitTo(*data,Extended(),NumCPU(24),Save()) ;

    //nsigtotdecay.setConstant(kTRUE);

    RooFitResult* fitresp1n;
    fitresp1n=final_decay1nmodel.fitTo(*datap1n,Extended(),NumCPU(24),Save()) ;

    RooFitResult* fitresp2n;
    //fitresp2n=final_decay2nmodel.fitTo(*datap2n,Extended(),NumCPU(24),Save()) ;


    TCanvas *c1 = new TCanvas("c1","c1",1200, 800);
    c1->Divide(3,1);

    RooPlot* xframe0 = x.frame(Title("fit")) ;
    data->plotOn(xframe0,Binning(200)) ;
    final_totdecaymodel.plotOn(xframe0) ;

    RooPlot* xframe1 = x.frame(Title("fit")) ;
    datap1n->plotOn(xframe1,Binning(200)) ;
    final_decay1nmodel.plotOn(xframe1) ;

    RooPlot* xframe2 = x.frame(Title("fit")) ;
    datap2n->plotOn(xframe2,Binning(200)) ;
    final_decay2nmodel.plotOn(xframe2) ;


    c1->cd(1);
    gPad->SetLeftMargin(0.15) ; xframe0->GetYaxis()->SetTitleOffset(1.4) ; xframe0->Draw() ;
    c1->cd(2);
    gPad->SetLeftMargin(0.15) ; xframe1->GetYaxis()->SetTitleOffset(1.4) ; xframe1->Draw() ;
    c1->cd(3);
    gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.4) ; xframe2->Draw() ;

}
