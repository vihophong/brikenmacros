
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
#include "fitFbkg.cxx"


#include <fstream>

using namespace RooFit ;
using namespace RooStats;

const Int_t ncpu = 16;
Long64_t nentrieslimit = -5000;//set negative value for fitting to all entries

Double_t p_deadtime=0.08;
Double_t p_timerange=10;

using namespace std;
void mlhfitsim()
{
    // Import tree to total decay data
    TFile *f=TFile::Open("outhist.root");
    TTree* tree;
    f->GetObject("tree",tree);
    if (nentrieslimit>0) tree->SetEntries(nentrieslimit);

    RooRealVar x("x","x",p_deadtime,p_timerange) ;
    RooRealVar xbkg("x","x",-p_timerange,-p_deadtime) ;

    //! define discrete variable y
    RooCategory y("y","y");
    y.defineType("0neu",0);
    y.defineType("1neu",1);
    y.defineType("2neu",2);


    RooRealVar bkg1nratio("bkg1nratio","bkg1nratio",0.5,0,1) ;
    RooRealVar bkg2nratio("bkg2nratio","bkg2nratio",0.5,0,1) ;

    RooRealVar slope1("slope1","slope1",0.,-0.1,0.1) ;
    RooRealVar slope2("slope2","slope2",0.,-0.1,0.1) ;
    RooRealVar slope3("slope3","slope3",0.,-0.1,0.1) ;


    //! bkg pdf
    fitFbkg bkgmodel("bkgmodel","bkgmodel",xbkg,y,bkg1nratio,bkg2nratio,slope1,slope2,slope3);

    //! prepare data set for bkg fit
    RooDataSet* data2=new RooDataSet("data","data",RooArgSet(xbkg,y),Import(*tree)) ;
    data2->Print() ;

    //! fit background
    bkgmodel.fitTo(*data2,NumCPU(ncpu),Save()) ;


    RooRealVar slope1pos("slope1pos","slope1pos",-slope1.getVal(),-0.1,0.1) ;
    RooRealVar slope2pos("slope2pos","slope2pos",-slope2.getVal(),-0.1,0.1) ;
    RooRealVar slope3pos("slope3pos","slope3pos",-slope3.getVal(),-0.1,0.1) ;

    //! set background ratio and slope as constant
    bkg1nratio.setConstant();
    bkg2nratio.setConstant();

    slope1pos.setConstant();
    slope2pos.setConstant();
    slope3pos.setConstant();

    //! bkg pdf positive
    fitFbkg bkgmodelpos("bkgmodelpos","bkgmodelpos",x,y,bkg1nratio,bkg2nratio,slope1pos,slope2pos,slope3pos);

    RooRealVar ineueff("ineueff","ineueff",0.62,0.,1.) ;
    ineueff.setConstant(kTRUE);


    //! read input
    makepath("In134short.txt");
    //makepath("In134full.txt");
    //makepath("Cd130full_isomer.txt");

    for (int i=0;i<knri*4;i++){
        //if (isparmsfix[i]!=2){
            cout<<"Parms for unbinned fit, index1 = "<<i<<" : ";
            cout<<riname[i]<<"\t"<<parms[i]<<"\t"<<parmserr[i]<<"\t"<<parmsmin[i]<<"\t"<<parmsmax[i]<<"\t"<<isparmsfix[i]<<endl;
        //}
    }

    RooRealVar* p[151];
    // decay parameters
    for (int i=0;i<knri*4;i++){
        p[i]=new RooRealVar(Form("p%d",i),Form("p%d",i),parms[i],parmsmin[i],parmsmax[i]);
        if (isparmsfix[i]!=0){
            cout<<"fixed value p"<<i<<"\tval="<<parms[i]<<endl;
            p[i]->setConstant(kTRUE);
        }else{
            cout<<"variable value p"<<i<<"\tval "<<parms[i]<<"\tmin="<<parmsmin[i]<<"\tmax="<<parmsmax[i]<<endl;
        }
    }

    //! initial activity set to 1
    p[knri*4]=new RooRealVar(Form("p%d",knri*4),Form("p%d",knri*4),1,0,2);
    p[knri*4]->setConstant();

    char tempchar1[1000];
    sprintf(tempchar1,"hdecay");
    TH1F* hdecay=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecay1nbwd");
    TH1F* hdecay1nbwd=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecaygt0nbwd");
    TH1F* hdecaygt0nbwd=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecay2nbwd");
    TH1F* hdecay2nbwd=(TH1F*) gDirectory->Get(tempchar1);

    //! Calculate random coincidence paramters
    Double_t n1nbwd=(Double_t) hdecay1nbwd->GetEntries();
    Double_t gt0nbwd=(Double_t) hdecaygt0nbwd->GetEntries();
    Double_t n2nbwd=(Double_t) hdecay2nbwd->GetEntries();
    Double_t nball=(Double_t) hdecay->GetEntries();

    Double_t randcoinf1n=n1nbwd/nball;
    Double_t randcoinfgt0n=gt0nbwd/nball;
    Double_t randcoinf2n=n2nbwd/nball;

    // random coincicence parameters
    p[knri*4+1]=new RooRealVar(Form("p%d",knri*4+1),Form("p%d",knri*4+1),randcoinf1n,randcoinf1n/3,randcoinf1n*3);
    p[knri*4+2]=new RooRealVar(Form("p%d",knri*4+2),Form("p%d",knri*4+2),randcoinfgt0n,randcoinfgt0n/3,randcoinfgt0n*3);
    p[knri*4+3]=new RooRealVar(Form("p%d",knri*4+3),Form("p%d",knri*4+3),randcoinf2n,randcoinf2n/3,randcoinf2n*3);
    p[knri*4+1]->setConstant(kTRUE);
    p[knri*4+2]->setConstant(kTRUE);
    p[knri*4+3]->setConstant(kTRUE);

    cout<<"randcoinf1n="<<randcoinf1n<<endl;
    cout<<"randcoinfgt0n="<<randcoinfgt0n<<endl;
    cout<<"randcoinf2n="<<randcoinf2n<<endl;


    for (int i=knri*4+4;i<151;i++){
        p[i]=new RooRealVar(Form("p%d",i),Form("p%d",i),1,0,2);
        p[i]->setConstant();
    }


    //! sig pdf
    fitF totdecaymodel("totdecaymodel","totdecaymodel",x,y,ineueff,p);

    //! set background/signal counts
    Double_t nnsig=tree->Draw("",Form("x>%f&&x<%f",p_deadtime,p_timerange),"goff");
    Double_t nnbkg=tree->Draw("",Form("x<%f&&x>%f",-p_deadtime,-p_timerange),"goff");
    nnsig=nnsig-nnbkg;

    RooRealVar nbkg("nbkg","nbkg",nnbkg,0,nnbkg*10) ;
    RooRealVar nsig("nsig","nsig",nnsig,0,nnsig*10) ;
    nbkg.setConstant();

    RooAddPdf final_pdf("final_pdf","final pdf",RooArgList(totdecaymodel,bkgmodelpos),RooArgList(nsig,nbkg));


    cout<<"\n\n\n********************************\nAttaching data and performing the fit\n********************************\n\n\n"<<endl;

    //! Toy MC dataset (not working!, just for testing)
    //RooDataSet* data = final_pdf.generate(RooArgSet(x,y),1) ;


    //! Prepare data set for fitting forward correlated data
    RooDataSet* data=new RooDataSet("data","data",RooArgSet(x,y),Import(*tree)) ;
    data->Print() ;

    cout<<"\n********************************\n"<<endl;
    // Test Fit
    RooFitResult* fitres;
    fitres=final_pdf.fitTo(*data,NumCPU(ncpu),Save()) ;

    cout<<"\n\n\n********************************\nFitting DONE!, now plotting\n********************************\n\n\n"<<endl;


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

}
