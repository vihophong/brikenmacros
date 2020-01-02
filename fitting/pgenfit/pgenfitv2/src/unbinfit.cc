//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * Copyright@2019 Vi Ho Phong, email: phong@ribf.riken.jp           *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications.                    *
// ********************************************************************
//
/// \file unbinfit.cc
/// \brief Implementation of the unbinfit class

#include "unbinfit.hh"
#include <iostream>

#include "RooDataSet.h"


#include "fitF.hh"

#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooPlot.h"

#include "RooAddPdf.h"
#include "TString.h"


using namespace RooFit;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unbinfit::unbinfit()
{
    fnentrieslimit=5000;

    p_deadtime=0.05;
    p_timerange=10;

    ncpu=2;

    finputData=new char[1000];
    finputParms=new char[1000];

    path=new decaypath;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unbinfit::~unbinfit()
{
    delete finputData;
    delete finputParms;
    delete path;
    delete tree;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::Init(char* inputParms, char* inputData)
{
    sprintf(finputParms,"%s",inputParms);
    sprintf(finputData,"%s",inputData);
    std::clog<< __PRETTY_FUNCTION__ <<"read input files:"<<
               std::endl<<finputParms<<
               std::endl<<finputData<<std::endl;

    path->Init(finputParms);
    path->makePath();
    path->printMember();
    path->printPath();
    path->writePath();
    path->drawPath((char*)TString("outdecayroutes.root").Data());

    // Import tree to total decay data
    TFile *f=TFile::Open(finputData);
    f->GetObject("tree",tree);
    if (fnentrieslimit>0) tree->SetEntries(fnentrieslimit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::Run(char* outputFileName)
{
    //! sig pdf
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
    RooDataSet* data2=new RooDataSet("data2","data2",RooArgSet(xbkg,y),Import(*tree)) ;
    data2->Print() ;

    //! fit background
    bkgmodel.fitTo(*data2,NumCPU(ncpu),Save()) ;

    //! set positive background parms
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

    RooRealVar ineueff("ineueff","ineueff",0.6441523116,0.,1.) ;
    ineueff.setConstant(kTRUE);


    //! set all decay parameters
    RooAbsReal* p[path->getNMember()*4+4];
    RooRealVar* pvar[path->getNMember()*4+4];

    for (Int_t i=0;i<path->getNMember()*4+4;i++){
        p[i]=new RooRealVar(Form("p%i",i),Form("p%i",i),0.5,0,1);
        pvar[i]=(RooRealVar*) p[i];
        pvar[i]->setConstant();
    }

    // decay parameters
    for (int i=0;i<path->getNMember();i++){
        p[i]=new RooRealVar(Form("p%d",i),Form("p%d",i),path->getMember(i)->decay_lamda,path->getMember(i)->decay_lamdalow,path->getMember(i)->decay_lamdaup);
        pvar[i]=(RooRealVar*) p[i];
        if (path->getMember(i)->is_decay_hl_fix!=0){
            cout<<"fixed decay rate p"<<i<<"\tval="<<path->getMember(i)->decay_lamda<<endl;
            pvar[i]->setConstant(kTRUE);
        }else{
            cout<<"variable decay rate p"<<i<<"\tval "<<path->getMember(i)->decay_lamda<<"\tmin="<<path->getMember(i)->decay_lamdalow<<"\tmax="<<path->getMember(i)->decay_lamdaup<<endl;
        }

        p[path->getNMember()+i]=new RooRealVar(Form("p%d",path->getNMember()+i),Form("p%d",path->getNMember()+i),path->getMember(i)->decay_p1n,path->getMember(i)->decay_p1nlow,path->getMember(i)->decay_p1nup);
        pvar[path->getNMember()+i]=(RooRealVar*) p[path->getNMember()+i];
        if (path->getMember(i)->is_decay_p1n_fix!=0){
            cout<<"fixed P1n p"<<path->getNMember()+i<<"\tval="<<path->getMember(i)->decay_p1n<<endl;
            pvar[path->getNMember()+i]->setConstant(kTRUE);
        }else{
            cout<<"variable P1n p"<<path->getNMember()+i<<"\tval "<<path->getMember(i)->decay_p1n<<"\tmin="<<path->getMember(i)->decay_p1nlow<<"\tmax="<<path->getMember(i)->decay_p1nup<<endl;
        }

        p[path->getNMember()*2+i]=new RooRealVar(Form("p%d",path->getNMember()*2+i),Form("p%d",path->getNMember()*2+i),path->getMember(i)->decay_p2n,path->getMember(i)->decay_p2nlow,path->getMember(i)->decay_p2nup);
        pvar[path->getNMember()*2+i]=(RooRealVar*) p[path->getNMember()*2+i];
        if (path->getMember(i)->is_decay_p2n_fix!=0){
            cout<<"fixed P2n p"<<path->getNMember()*2+i<<"\tval="<<path->getMember(i)->decay_p2n<<endl;
            pvar[path->getNMember()*2+i]->setConstant(kTRUE);
        }else{
            cout<<"variable P2n p"<<path->getNMember()*2+i<<"\tval "<<path->getMember(i)->decay_p2n<<"\tmin="<<path->getMember(i)->decay_p2nlow<<"\tmax="<<path->getMember(i)->decay_p2nup<<endl;
        }

        p[path->getNMember()*3+i]=new RooRealVar(Form("p%d",path->getNMember()*3+i),Form("p%d",path->getNMember()*3+i),path->getMember(i)->population_ratio,path->getMember(i)->population_ratiolow,path->getMember(i)->population_ratioup);
        pvar[path->getNMember()*3+i]=(RooRealVar*) p[path->getNMember()*3+i];
        if (path->getMember(i)->is_population_ratio_fix!=0){
            cout<<"fixed production ratio p"<<path->getNMember()*3+i<<"\tval="<<path->getMember(i)->population_ratio<<endl;
            pvar[path->getNMember()*3+i]->setConstant(kTRUE);
        }else{
            cout<<"variable production ratio p"<<path->getNMember()*3+i<<"\tval "<<path->getMember(i)->population_ratio<<"\tmin="<<path->getMember(i)->population_ratiolow<<"\tmax="<<path->getMember(i)->population_ratioup<<endl;
        }
    }

    //! initial activity set to 1
    p[path->getDecayPath()->nri*4]=new RooRealVar(Form("p%d",path->getDecayPath()->nri*4),Form("p%d",path->getDecayPath()->nri*4),1,0,2);
    pvar[path->getDecayPath()->nri*4]=(RooRealVar*) p[path->getDecayPath()->nri*4];
    pvar[path->getDecayPath()->nri*4]->setConstant();

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
    p[path->getDecayPath()->nri*4+1]=new RooRealVar(Form("p%d",path->getDecayPath()->nri*4+1),Form("p%d",path->getDecayPath()->nri*4+1),randcoinf1n,randcoinf1n/3,randcoinf1n*3);
    p[path->getDecayPath()->nri*4+2]=new RooRealVar(Form("p%d",path->getDecayPath()->nri*4+2),Form("p%d",path->getDecayPath()->nri*4+2),randcoinfgt0n,randcoinfgt0n/3,randcoinfgt0n*3);
    p[path->getDecayPath()->nri*4+3]=new RooRealVar(Form("p%d",path->getDecayPath()->nri*4+3),Form("p%d",path->getDecayPath()->nri*4+3),randcoinf2n,randcoinf2n/3,randcoinf2n*3);
    pvar[path->getDecayPath()->nri*4+1]=(RooRealVar*) p[path->getDecayPath()->nri*4+1];
    pvar[path->getDecayPath()->nri*4+2]=(RooRealVar*) p[path->getDecayPath()->nri*4+2];
    pvar[path->getDecayPath()->nri*4+3]=(RooRealVar*) p[path->getDecayPath()->nri*4+3];
    pvar[path->getDecayPath()->nri*4+1]->setConstant(kTRUE);
    pvar[path->getDecayPath()->nri*4+2]->setConstant(kTRUE);
    pvar[path->getDecayPath()->nri*4+3]->setConstant(kTRUE);


    fitF totdecaymodel("totdecaymodel","totdecaymodel",x,y,ineueff,p);

    //! set background/signal counts
    Double_t nnsig=tree->Draw("",Form("x>%f&&x<%f",p_deadtime,p_timerange),"goff");
    Double_t nnbkg=tree->Draw("",Form("x<%f&&x>%f",-p_deadtime,-p_timerange),"goff");
    nnsig=nnsig-nnbkg;

    RooRealVar nbkg("nbkg","nbkg",nnbkg,nnbkg/3,nnbkg*3) ;
    RooRealVar nsig("nsig","nsig",nnsig,nnsig/3,nnsig*3) ;
    nbkg.setConstant();

    RooAddPdf final_pdf("final_pdf","final pdf",RooArgList(totdecaymodel,bkgmodelpos),RooArgList(nsig,nbkg));

////! Toy MC dataset (just for testing)
//    RooDataSet* datamc = final_pdf.generate(RooArgSet(x,y),10000) ;
//    datamc->Print() ;

    //! Prepare data set for fitting forward correlated data
    RooDataSet* data=new RooDataSet("data","data",RooArgSet(x,y),Import(*tree)) ;
    data->Print() ;

    RooFitResult* fitres;
    fitres=final_pdf.fitTo(*data,NumCPU(ncpu),Save());

    TCanvas *c1 = new TCanvas("c1","c1",1200, 800);
    c1->Divide(3,1);
    c1->cd(1);
    RooPlot* xframe0 = x.frame(Title("0 neutron fit")) ;
    data->plotOn(xframe0,Cut("y==y::0neu"),Binning(500)) ;
    final_pdf.plotOn(xframe0,Slice(y,"0neu")) ;
    xframe0->Draw() ;
    c1->cd(2);
    RooPlot* xframe1 = x.frame(Title("1 neutron fit")) ;
    data->plotOn(xframe1,Cut("y==y::1neu"),Binning(500)) ;
    final_pdf.plotOn(xframe1,Slice(y,"1neu")) ;
    xframe1->Draw() ;
    c1->cd(3);
    RooPlot* xframe2 = x.frame(Title("2 neutron fit")) ;
    data->plotOn(xframe2,Cut("y==y::2neu"),Binning(500)) ;
    final_pdf.plotOn(xframe2,Slice(y,"2neu")) ;
    xframe2->Draw() ;
    c1->SaveAs(outputFileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
