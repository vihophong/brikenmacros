#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooCategory.h"
#include "TRandom3.h"
#include "TFile.h"
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


#include "fitF.cxx"
#include "fitFskewed.cxx"


#include <fstream>
using namespace RooFit ;
using namespace RooStats;


void fitsplitnormal()
{
    RooRealVar x("x","x",0,1) ;
    RooRealVar mean("mean","mean",0.5,0.,1.) ;
    RooRealVar sig1("sig1","sig1",0.05,0.,1.) ;
    RooRealVar sig2("sig2","sig2",0.2,0.,1.) ;

    fitF model("totdecaymodel","totdecaymodel",x,mean,sig1,sig2);

    TTree* tree;
    TFile *f=TFile::Open("test.root");
    f->GetObject("tree",tree);

    //RooDataSet* data = model.generate(x,100000) ;
    RooDataSet* data=new RooDataSet("data","data",x,Import(*tree));
    data->Print() ;

    RooFitResult* fitres;
    fitres=model.fitTo(*data,NumCPU(2),Save()) ;
    TCanvas *c1 = new TCanvas("c1","c1",1200, 800);
    RooPlot* xframe0 = x.frame(Title("fit")) ;
    data->plotOn(xframe0,Binning(500)) ;
    model.plotOn(xframe0) ;
    xframe0->Draw() ;
}

void fitskewed()
{
    RooRealVar x("x","x",0,3) ;
    RooRealVar mean("mean","mean",0.5,0.,1.) ;
    RooRealVar sig1("sig1","sig1",0.05,0.,1.) ;
    RooRealVar sig2("sig2","sig2",0.2,0.,1.) ;

    fitFskewed model("totdecaymodel","totdecaymodel",x,mean,sig1,sig2);

    TTree* tree;
    TFile *f=TFile::Open("test.root");
    f->GetObject("tree",tree);

    //RooDataSet* data = model.generate(x,100000) ;
    RooDataSet* data=new RooDataSet("data","data",x,Import(*tree));
    data->Print() ;

    RooFitResult* fitres;
    fitres=model.fitTo(*data,NumCPU(2),Save()) ;
    TCanvas *c1 = new TCanvas("c1","c1",1200, 800);
    RooPlot* xframe0 = x.frame(Title("fit")) ;
    data->plotOn(xframe0,Binning(500)) ;
    model.plotOn(xframe0) ;
    xframe0->Draw() ;
}
void fitasymdis()
{
  fitsplitnormal();
  //fitskewed();
}
