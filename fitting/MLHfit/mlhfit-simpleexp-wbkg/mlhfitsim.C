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
#include "fitFbkg.cxx"
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

#include <fstream>

using namespace RooFit ;
using namespace RooStats;

Double_t p_deadtime=0;
Double_t p_timerange=10;

using namespace std;

void mlhfitsim()
{
    //! define variables
    RooRealVar x("x","x",p_deadtime,p_timerange) ;
    RooRealVar xbkg("x","x",-p_timerange,-p_deadtime) ;


    RooRealVar* p[7];//hl.p1n,p2n

    Double_t hl=TMath::Log(2)/0.128;
    Double_t hlplus=TMath::Log(2)/0.01;
    Double_t hlminus=TMath::Log(2)/1.;

    p[0]=new RooRealVar(Form("p%d",0),Form("p%d",0),hl,hlminus,hlplus);
    p[1]=new RooRealVar(Form("p%d",1),Form("p%d",1),0.5,0.,1.);
    p[2]=new RooRealVar(Form("p%d",2),Form("p%d",2),0.25,0.,1.);

    p[3]=new RooRealVar(Form("p%d",3),Form("p%d",3),0.5,0.,1.);//p1nbkg ratio
    p[4]=new RooRealVar(Form("p%d",4),Form("p%d",4),0.5,0.,1.);//p2nbkg ratio



    //p[1]->setConstant(kTRUE);
    //p[2]->setConstant(kTRUE);


    //! define discrete variable y
    RooCategory y("y","y");
    y.defineType("0neu",0);
    y.defineType("1neu",1);
    y.defineType("2neu",2);


    // bkg pdf
    fitFbkg bkgmodel("bkgmodel","bkgmodel",xbkg,y,*p[3],*p[4]);


    // sig pdf
    fitF totdecaymodel("totdecaymodel","totdecaymodel",x,y,*p[0],*p[1],*p[2]);

    // Import tree to total decay data
    TFile *f=TFile::Open("outhist.root");
    TTree* tree;
    f->GetObject("tree",tree);

    RooDataSet* data=new RooDataSet("data","data",RooArgSet(x,y),Import(*tree)) ;
    data->Print() ;

    RooDataSet* data2=new RooDataSet("data","data",RooArgSet(xbkg,y),Import(*tree)) ;
    data2->Print() ;

    // Fit
    RooFitResult* fitres;
    fitres=bkgmodel.fitTo(*data2,NumCPU(2),Save()) ;

    p[3]->setConstant();
    p[4]->setConstant();

    Double_t ntotal=tree->Draw("","x>0","goff");
    Double_t nbkg=tree->Draw("","x<0","goff");
    ntotal=ntotal-nbkg;

    p[5]=new RooRealVar(Form("p%d",5),Form("p%d",5),nbkg,0.,nbkg*3);//nbkg
    p[6]=new RooRealVar(Form("p%d",6),Form("p%d",6),ntotal,0.,ntotal*3);//nsig

    p[5]->setConstant();

    RooAddPdf final_pdf("final_pdf","final pdf",RooArgList(totdecaymodel,bkgmodel),RooArgList(*p[6],*p[5]));

    fitres=final_pdf.fitTo(*data,NumCPU(2),Save()) ;



    TCanvas *c1 = new TCanvas("c1","c1",1200, 800);
    c1->Divide(2,3);
    c1->cd(1);
    RooPlot* xframe0 = x.frame(Title("0 neutron fit")) ;
    data->plotOn(xframe0,Cut("y==y::0neu"),Binning(500)) ;
    final_pdf.plotOn(xframe0,Slice(y,"0neu")) ;
    xframe0->Draw() ;
    c1->cd(3);
    RooPlot* xframe1 = x.frame(Title("1 neutron fit")) ;
    data->plotOn(xframe1,Cut("y==y::1neu"),Binning(500)) ;
    final_pdf.plotOn(xframe1,Slice(y,"1neu")) ;
    xframe1->Draw() ;
    c1->cd(5);
    RooPlot* xframe2 = x.frame(Title("2 neutron fit")) ;
    data->plotOn(xframe2,Cut("y==y::2neu"),Binning(500)) ;
    final_pdf.plotOn(xframe2,Slice(y,"2neu")) ;
    xframe2->Draw() ;



    c1->cd(2);
    RooPlot* xframe3 = xbkg.frame(Title("0 neutron fit")) ;
    data2->plotOn(xframe3,Cut("y==y::0neu"),Binning(500)) ;
    bkgmodel.plotOn(xframe3,Slice(y,"0neu")) ;
    xframe3->Draw() ;
    c1->cd(4);
    RooPlot* xframe4 = xbkg.frame(Title("1 neutron fit")) ;
    data2->plotOn(xframe4,Cut("y==y::1neu"),Binning(500)) ;
    bkgmodel.plotOn(xframe4,Slice(y,"1neu")) ;
    xframe4->Draw() ;
    c1->cd(6);
    RooPlot* xframe5 = xbkg.frame(Title("2 neutron fit")) ;
    data2->plotOn(xframe5,Cut("y==y::2neu"),Binning(500)) ;
    bkgmodel.plotOn(xframe5,Slice(y,"2neu")) ;
    xframe5->Draw() ;




    ofstream str("out.txt",ios::app);
    str<<p[0]->getVal()<<"\t"<<p[1]->getVal()<<"\t"<<p[2]->getVal()<<endl;

}
