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

const Int_t ncpu = 8;
Long64_t nentrieslimit = 10000;//set negative value for fitting to all entries

Double_t p_deadtime=0.;
Double_t p_timerange=10;

using namespace std;
void mlhfitsim()
{
    Int_t fitopt=0;
    // Import tree to total decay data
    TFile *f=TFile::Open("fittrees/Cd130.root");
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

    slope1pos.setError(slope1.getError());
    slope2pos.setError(slope2.getError());
    slope3pos.setError(slope3.getError());


    RooRealVar ineueff("ineueff","ineueff",0.62,0.,1.) ;

    //! read input
    makepath("Cd130full_isomer.txt");
    //makepath("In134short.txt");

    for (int i=0;i<knri*4;i++){
            cout<<"Parms for unbinned fit, index1 = "<<i<<" : ";
            cout<<riname[i]<<"\t"<<parms[i]<<"\t"<<parmserr[i]<<"\t"<<parmsmin[i]<<"\t"<<parmsmax[i]<<"\t"<<isparmsfix[i]<<endl;
    }

    RooRealVar* p[151];
    //! decay parameters
    for (int i=0;i<knri*4;i++){
        p[i]=new RooRealVar(Form("p%d",i),Form("p%d",i),parms[i],parmsmin[i],parmsmax[i]);
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
    p[knri*4+1]->setError(randcoinf1n*sqrt(1/n1nbwd+1/nball));
    p[knri*4+2]->setError(randcoinfgt0n*sqrt(1/gt0nbwd+1/nball));
    p[knri*4+3]->setError(randcoinf2n*sqrt(1/n2nbwd+1/nball));


    //! set the rest of unused paramters
    for (int i=knri*4+4;i<151;i++){
        p[i]=new RooRealVar(Form("p%d",i),Form("p%d",i),1,0,2);
        p[i]->setConstant();
    }

    //! set background/signal counts for final pdf
    Double_t nnsig=tree->Draw("",Form("x>%f",p_deadtime),"goff");
    Double_t nnbkg=tree->Draw("",Form("x<%f",-p_deadtime),"goff");
    nnsig=nnsig-nnbkg;

    RooRealVar nbkg("nbkg","nbkg",nnbkg,0,nnbkg*10) ;
    RooRealVar nsig("nsig","nsig",nnsig,0,nnsig*10) ;

    nbkg.setError(sqrt(nnbkg));



    //! set constant parameters and constrains

    //! define external contrains
    RooArgSet externalconstrains;
    RooGaussian* pconstr[knri*4+4];
    RooGaussian* neueffconstr=new RooGaussian("neueffconstr","neueffconstr",ineueff,RooConst(ineueff.getVal()),RooConst(0.01));

    RooGaussian* bkg1nratioconstr=new RooGaussian("bkg1nratioconstr","bkg1nratioconstr",bkg1nratio,RooConst(bkg1nratio.getVal()),RooConst(bkg1nratio.getError()));
    RooGaussian* bkg2nratioconstr=new RooGaussian("bkg2nratioconstr","bkg2nratioconstr",bkg2nratio,RooConst(bkg2nratio.getVal()),RooConst(bkg2nratio.getError()));

    RooGaussian* slope1posconstr=new RooGaussian("slope1posconstr","slope1posconstr",slope1pos,RooConst(slope1pos.getVal()),RooConst(slope1pos.getError()));
    RooGaussian* slope2posconstr=new RooGaussian("slope2posconstr","slope2posconstr",slope2pos,RooConst(slope2pos.getVal()),RooConst(slope2pos.getError()));
    RooGaussian* slope3posconstr=new RooGaussian("slope3posconstr","slope3posconstr",slope3pos,RooConst(slope3pos.getVal()),RooConst(slope3pos.getError()));


    RooGaussian* nbkgconstr=new RooGaussian("nbkgconstr","nbkgconstr",nbkg,RooConst(nbkg.getVal()),RooConst(nbkg.getError()));



     for (int i=0;i<knri*4;i++)
         pconstr[i]=new RooGaussian(Form("p%dconstr",i),Form("p%dconstr",i),*p[i],RooConst(parms[i]),RooConst(parmserr[i]));

     pconstr[knri*4+1]=new RooGaussian(Form("p%dconstr",knri*4+1),Form("p%dconstr",knri*4+1),*p[knri*4+1],RooConst(p[knri*4+1]->getVal()),RooConst(p[knri*4+1]->getError()));
     pconstr[knri*4+2]=new RooGaussian(Form("p%dconstr",knri*4+2),Form("p%dconstr",knri*4+2),*p[knri*4+2],RooConst(p[knri*4+2]->getVal()),RooConst(p[knri*4+2]->getError()));
     pconstr[knri*4+3]=new RooGaussian(Form("p%dconstr",knri*4+3),Form("p%dconstr",knri*4+3),*p[knri*4+3],RooConst(p[knri*4+3]->getVal()),RooConst(p[knri*4+3]->getError()));

    // decay parameters
    for (int i=0;i<knri*4;i++){
        // stuff for isomers
        Bool_t flag_continue=false;
        if (flag_sum_isomer_ratio){
            if (i>knri*3-1){
                for (Int_t j=0;j<nisomers;j++){
                    if ((i-knri*3) == groundstate[j]){
                        p[i]->setConstant();
                        flag_continue=true;
                    }
                }
            }
        }
        if (flag_continue) continue;

        //set constant or constrains

        if (isparmsfix[i]!=0){
            if (fitopt==0) {
                //if (i<20)
                p[i]->setConstant();
                cout<<"set Constant for parameter p"<<i<<"\tmean="<<parms[i]<<"\tstd="<<parmserr[i]<<endl;
            }else{
                if (parmserr[i]==0){
                    cout<<"set Constant for parameter p"<<i<<"\tmean="<<parms[i]<<"\tstd="<<parmserr[i]<<endl;
                    p[i]->setConstant();
                }else{
                    if (isparmsfix[i]==1) {
                        //p[i]->setConstant();
                        externalconstrains.add(*pconstr[i]);
                        cout<<"set constrain model for p"<<i<<"\tmean="<<parms[i]<<"\tstd="<<parmserr[i]<<endl;
                    }else{
                        p[i]->setConstant();
                        cout<<"set Constant for parameter p"<<i<<"\tmean="<<parms[i]<<"\tstd="<<parmserr[i]<<endl;
                    }
                }
            }
        }else{
            cout<<"variable value p"<<i<<"\tval "<<parms[i]<<"\tmin="<<parmsmin[i]<<"\tmax="<<parmsmax[i]<<endl;
        }
    }

    //! Error propagation using external contrains for other parameters


    if (fitopt==1){
         externalconstrains.add(*nbkgconstr);
         externalconstrains.add(*bkg1nratioconstr);
         externalconstrains.add(*bkg2nratioconstr);


         externalconstrains.add(*slope1posconstr);
         externalconstrains.add(*slope2posconstr);
         externalconstrains.add(*slope3posconstr);

         externalconstrains.add(*neueffconstr);

         externalconstrains.add(*pconstr[knri*4+1]);
         externalconstrains.add(*pconstr[knri*4+2]);
         externalconstrains.add(*pconstr[knri*4+3]);
     }else{
         // constant
         nbkg.setConstant();
         bkg1nratio.setConstant();
         bkg2nratio.setConstant();

         slope1pos.setConstant();
         slope2pos.setConstant();
         slope3pos.setConstant();

         ineueff.setConstant(kTRUE);

         p[knri*4+1]->setConstant(kTRUE);
         p[knri*4+2]->setConstant(kTRUE);
         p[knri*4+3]->setConstant(kTRUE);
     }





    //! pdf models
    //! bkg pdf positive
    fitFbkg bkgmodelpos("bkgmodelpos","bkgmodelpos",x,y,bkg1nratio,bkg2nratio,slope1pos,slope2pos,slope3pos);
    //! sig pdf
    fitF totdecaymodel("totdecaymodel","totdecaymodel",x,y,ineueff,p);
    RooAddPdf final_pdf("final_pdf","final pdf",RooArgList(totdecaymodel,bkgmodelpos),RooArgList(nsig,nbkg));


    cout<<"\n\n\n********************************\nAttaching data and performing the fit\n********************************\n\n\n"<<endl;

    //! Toy MC dataset (not working!, just for testing)
    //RooDataSet* data = final_pdf.generate(RooArgSet(x,y),100) ;

    //! Prepare data set for fitting forward correlated data
    RooDataSet* data=new RooDataSet("data","data",RooArgSet(x,y),Import(*tree)) ;
    data->Print() ;

    cout<<"\n********************************\n"<<endl;
    // Test Fit
    RooFitResult* fitres;

    if (fitopt==0)
        fitres=final_pdf.fitTo(*data,NumCPU(ncpu),Save()) ;
    else
        fitres=final_pdf.fitTo(*data,ExternalConstraints(externalconstrains),NumCPU(ncpu),Save()) ;

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

    TFile* outfileroot=new TFile("test.root","recreate");
    c1->Write();
    outfileroot->Close();
}
