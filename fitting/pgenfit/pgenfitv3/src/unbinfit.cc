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
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "TMath.h"


using namespace RooFit;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unbinfit::unbinfit()
{
    ffitopt=2;// only 0, 1 and 2

    fnentrieslimit=300000;

    p_deadtime=0.08;
    p_timerange=10;

    ncpu=8;

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
    //!*****************************************
    //! Prepare X,y
    //! *****************************************

    // sig pdf
    RooRealVar x("x","x",p_deadtime,p_timerange) ;
    RooRealVar xbkg("x","x",-p_timerange,-p_deadtime) ;

    // define discrete variable y
    RooCategory y("y","y");
    y.defineType("0neu",0);
    y.defineType("1neu",1);
    y.defineType("2neu",2);


    //!*****************************************
    //! Prepare and fit negative background function
    //! *****************************************
    // bkg parameters
    RooRealVar bkg1nratio("bkg1nratio","bkg1nratio",0.5,0,1) ;
    RooRealVar bkg2nratio("bkg2nratio","bkg2nratio",0.5,0,1) ;

    RooRealVar slope1("slope1","slope1",0.,-0.1,0.1) ;
    RooRealVar slope2("slope2","slope2",0.,-0.1,0.1) ;
    RooRealVar slope3("slope3","slope3",0.,-0.1,0.1) ;

    // bkg pdf
    fitFbkg bkgmodel("bkgmodel","bkgmodel",xbkg,y,bkg1nratio,bkg2nratio,slope1,slope2,slope3);
    // prepare data set for bkg fit
    RooDataSet* data2=new RooDataSet("data2","data2",RooArgSet(xbkg,y),Import(*tree)) ;
    data2->Print() ;
    // fit background
    bkgmodel.fitTo(*data2,NumCPU(ncpu),Save()) ;

    //!*****************************************
    //! Prepare positive background function
    //! *****************************************

    // set positive background parms
    RooRealVar slope1pos("slope1pos","slope1pos",-slope1.getVal(),-0.1,0.1) ;
    RooRealVar slope2pos("slope2pos","slope2pos",-slope2.getVal(),-0.1,0.1) ;
    RooRealVar slope3pos("slope3pos","slope3pos",-slope3.getVal(),-0.1,0.1) ;

    // bkg pdf positive
    fitFbkg bkgmodelpos("bkgmodelpos","bkgmodelpos",x,y,bkg1nratio,bkg2nratio,slope1pos,slope2pos,slope3pos);

    //!*****************************************
    //! Initialize all parameters
    //! *****************************************
    // Declare parent neutron detection efficiency
    RooRealVar ineueff("ineueff","ineueff",0.6441523116,0.63,0.65) ;

    // Declare all parameters
    RooAbsReal* p[path->getNMember()*5+4];
    RooRealVar* pvar[path->getNMember()*5+4];

    // Initialize decay parameters
    for (int i=0;i<path->getNMember();i++){
        p[i]=new RooRealVar(Form("p%d",i),Form("p%d",i),path->getMember(i)->decay_lamda,path->getMember(i)->decay_lamdalow,path->getMember(i)->decay_lamdaup);
        pvar[i]=(RooRealVar*) p[i];

        p[path->getNMember()+i]=new RooRealVar(Form("p%d",path->getNMember()+i),Form("p%d",path->getNMember()+i),path->getMember(i)->decay_p1n,path->getMember(i)->decay_p1nlow,path->getMember(i)->decay_p1nup);
        pvar[path->getNMember()+i]=(RooRealVar*) p[path->getNMember()+i];

        p[path->getNMember()*2+i]=new RooRealVar(Form("p%d",path->getNMember()*2+i),Form("p%d",path->getNMember()*2+i),path->getMember(i)->decay_p2n,path->getMember(i)->decay_p2nlow,path->getMember(i)->decay_p2nup);
        pvar[path->getNMember()*2+i]=(RooRealVar*) p[path->getNMember()*2+i];

        p[path->getNMember()*3+i]=new RooRealVar(Form("p%d",path->getNMember()*3+i),Form("p%d",path->getNMember()*3+i),path->getMember(i)->population_ratio,path->getMember(i)->population_ratiolow,path->getMember(i)->population_ratioup);
        pvar[path->getNMember()*3+i]=(RooRealVar*) p[path->getNMember()*3+i];

        p[path->getNMember()*4+i]=new RooRealVar(Form("p%d",path->getNMember()*4+i),Form("p%d",path->getNMember()*4+i),0.9999999,0,1);
        pvar[path->getNMember()*4+i]=(RooRealVar*) p[path->getNMember()*4+i];
    }

    // initialize initial activity and set fixed to 1
    p[path->getNMember()*5]=new RooRealVar(Form("p%d",path->getNMember()*5),Form("p%d",path->getNMember()*5),1,0,2);
    pvar[path->getNMember()*5]=(RooRealVar*) p[path->getNMember()*5];
    pvar[path->getNMember()*5]->setConstant();


    // Get histograms to calculate random coincidence parameters
    char tempchar1[1000];
    sprintf(tempchar1,"hdecay");
    TH1F* hdecay=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecay1nbwd");
    TH1F* hdecay1nbwd=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecaygt0nbwd");
    TH1F* hdecaygt0nbwd=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecay2nbwd");
    TH1F* hdecay2nbwd=(TH1F*) gDirectory->Get(tempchar1);

    // Calculate random coincidence paramters
    Double_t n1nbwd=(Double_t) hdecay1nbwd->GetEntries();
    Double_t gt0nbwd=(Double_t) hdecaygt0nbwd->GetEntries();
    Double_t n2nbwd=(Double_t) hdecay2nbwd->GetEntries();
    Double_t nball=(Double_t) hdecay->GetEntries();

    Double_t randcoinf1n=n1nbwd/nball;
    Double_t randcoinfgt0n=gt0nbwd/nball;
    Double_t randcoinf2n=n2nbwd/nball;

    Double_t randcoinf1nerr=n1nbwd/nball*TMath::Sqrt(1/n1nbwd+1/nball);
    Double_t randcoinfgt0nerr=gt0nbwd/nball*TMath::Sqrt(1/gt0nbwd+1/nball);
    Double_t randcoinf2nerr=n2nbwd/nball*TMath::Sqrt(1/n2nbwd+1/nball);


    // Initialize random coincicence parameters
    p[path->getNMember()*5+1]=new RooRealVar(Form("p%d",path->getNMember()*5+1),Form("p%d",path->getNMember()*5+1),randcoinf1n,0,1);
    p[path->getNMember()*5+2]=new RooRealVar(Form("p%d",path->getNMember()*5+2),Form("p%d",path->getNMember()*5+2),randcoinfgt0n,0,1);
    p[path->getNMember()*5+3]=new RooRealVar(Form("p%d",path->getNMember()*5+3),Form("p%d",path->getNMember()*5+3),randcoinf2n,0,1);
    pvar[path->getNMember()*5+1]=(RooRealVar*) p[path->getNMember()*5+1];
    pvar[path->getNMember()*5+2]=(RooRealVar*) p[path->getNMember()*5+2];
    pvar[path->getNMember()*5+3]=(RooRealVar*) p[path->getNMember()*5+3];

    // set background/signal counts
    Double_t nnsig=tree->Draw("",Form("x>%f&&x<%f",p_deadtime,p_timerange),"goff");
    Double_t nnbkg=tree->Draw("",Form("x<%f&&x>%f",-p_deadtime,-p_timerange),"goff");
    nnsig=nnsig-nnbkg;

    RooRealVar nbkg("nbkg","nbkg",nnbkg,nnbkg/3,nnbkg*3) ;
    RooRealVar nsig("nsig","nsig",nnsig,nnsig/3,nnsig*3) ;


    //!******************************************
    //! Decide whether paramters are fixed
    //! *****************************************
    // decay parameters
    for (int i=0;i<path->getNMember();i++){
        if (ffitopt==2) {
            if (path->getMember(i)->is_decay_p1n_fix==2)
                pvar[path->getNMember()+i]->setConstant(kTRUE);
            if (path->getMember(i)->is_decay_p2n_fix==2)
                pvar[path->getNMember()*2+i]->setConstant(kTRUE);
            if (path->getMember(i)->is_population_ratio_fix==2)
                pvar[path->getNMember()*3+i]->setConstant(kTRUE);
            //if (path->getMember(i)->is_neueff_fcator_fix!=0)
                pvar[path->getNMember()*4+i]->setConstant(kTRUE);
        }else{
            if (path->getMember(i)->is_decay_hl_fix!=0)
                pvar[i]->setConstant(kTRUE);
            if (path->getMember(i)->is_decay_p1n_fix!=0)
                pvar[path->getNMember()+i]->setConstant(kTRUE);
            if (path->getMember(i)->is_decay_p2n_fix!=0)
                pvar[path->getNMember()*2+i]->setConstant(kTRUE);
            if (path->getMember(i)->is_population_ratio_fix!=0)
                pvar[path->getNMember()*3+i]->setConstant(kTRUE);
            //if (path->getMember(i)->is_neueff_fator_fix!=0)
                pvar[path->getNMember()*4+i]->setConstant(kTRUE);
        }
    }

    // others parameters
    slope1pos.setConstant();
    slope2pos.setConstant();
    slope3pos.setConstant();
    if (ffitopt==0){
        // neutron detection efficiency
        ineueff.setConstant(kTRUE);
        // backgrounds
        nbkg.setConstant(kTRUE);
        // background ratio and slope
        bkg1nratio.setConstant();
        bkg2nratio.setConstant();

        // random coincicence parameters
        pvar[path->getNMember()*5+1]->setConstant(kTRUE);//randcoinf1n
        pvar[path->getNMember()*5+2]->setConstant(kTRUE);//randcoinfgt0n
        pvar[path->getNMember()*5+3]->setConstant(kTRUE);//randcoinf2n
    }

    //!******************************************
    //! Define roogaussian for error propagation
    //! *****************************************
    RooArgSet externalconstrains;

    // initilze all roogausian
    RooGaussian* neueffconstr=new RooGaussian("neueffconstr","neueffconstr",ineueff,RooConst(ineueff.getVal()),RooConst(0.1));
    RooGaussian* nbkgconstr=new RooGaussian("nbkgconstr","nbkgconstr",nbkg,RooConst(nbkg.getVal()),RooConst(TMath::Sqrt(nbkg.getVal())));
    RooGaussian* bkg1nratiocnstr=new RooGaussian("bkg1nratiocnstr","bkg1nratiocnstr",bkg1nratio,RooConst(bkg1nratio.getVal()),RooConst(bkg1nratio.getError()));
    RooGaussian* bkg2nratiocnstr=new RooGaussian("bkg2nratiocnstr","bkg2nratiocnstr",bkg2nratio,RooConst(bkg2nratio.getVal()),RooConst(bkg2nratio.getError()));

    // decay paramters
    RooGaussian* pconstr[path->getNMember()*5+4];
    for (int i=0;i<path->getNMember();i++){
        pconstr[i]=new RooGaussian(Form("p%dconstr",i),Form("p%dconstr",i),*p[i],RooConst(pvar[i]->getVal()),RooConst(path->getMember(i)->decay_lamdaerr));
        pconstr[path->getNMember()+i]=new RooGaussian(Form("p%dconstr",path->getNMember()+i),Form("p%dconstr",path->getNMember()+i),*p[path->getNMember()+i],RooConst(pvar[path->getNMember()+i]->getVal()),RooConst(path->getMember(i)->decay_p1nerr));
        pconstr[path->getNMember()*2+i]=new RooGaussian(Form("p%dconstr",path->getNMember()*2+i),Form("p%dconstr",path->getNMember()*2+i),*p[path->getNMember()*2+i],RooConst(pvar[path->getNMember()*2+i]->getVal()),RooConst(path->getMember(i)->decay_p2nerr));
        pconstr[path->getNMember()*3+i]=new RooGaussian(Form("p%dconstr",path->getNMember()*3+i),Form("p%dconstr",path->getNMember()*3+i),*p[path->getNMember()*3+i],RooConst(pvar[path->getNMember()*3+i]->getVal()),RooConst(path->getMember(i)->population_ratioerr));
        pconstr[path->getNMember()*4+i]=new RooGaussian(Form("p%dconstr",path->getNMember()*4+i),Form("p%dconstr",path->getNMember()*4+i),*p[path->getNMember()*4+i],RooConst(pvar[path->getNMember()*4+i]->getVal()),RooConst(0.1));
    }
    // random coincidence paramters
    pconstr[path->getNMember()*5+1]=new RooGaussian(Form("p%dconstr",path->getNMember()*5+1),Form("p%dconstr",path->getNMember()*5+1),*p[path->getNMember()*5+1],RooConst(p[path->getNMember()*5+1]->getVal()),RooConst(randcoinf1nerr));
    pconstr[path->getNMember()*5+2]=new RooGaussian(Form("p%dconstr",path->getNMember()*5+2),Form("p%dconstr",path->getNMember()*5+2),*p[path->getNMember()*5+2],RooConst(p[path->getNMember()*5+2]->getVal()),RooConst(randcoinfgt0nerr));
    pconstr[path->getNMember()*5+3]=new RooGaussian(Form("p%dconstr",path->getNMember()*5+3),Form("p%dconstr",path->getNMember()*5+3),*p[path->getNMember()*5+3],RooConst(p[path->getNMember()*5+3]->getVal()),RooConst(randcoinf2nerr));

    // add to contrains set
    if (ffitopt==2) {
        for (int i=0;i<path->getNMember()*5;i++) {
            if (!pvar[i]->isConstant()&&i!=0&&i!=path->getNMember()&&i!=path->getNMember()*2) // let decay parameters of parent nuclei free
                externalconstrains.add(*pconstr[i]);
        }
    }


    // others parameters
    if (ffitopt!=0){
        // neutron detection efficiency and bkg parms
        externalconstrains.add(*neueffconstr);
        externalconstrains.add(*nbkgconstr);
        externalconstrains.add(*bkg1nratiocnstr);
        externalconstrains.add(*bkg2nratiocnstr);
        //random coincicence parameters
        externalconstrains.add(*pconstr[path->getNMember()*5+1]);
        externalconstrains.add(*pconstr[path->getNMember()*5+2]);
        externalconstrains.add(*pconstr[path->getNMember()*5+3]);
    }

    //!******************************************
    //! Construct final fit model
    //! *****************************************
    fitF totdecaymodel("totdecaymodel","totdecaymodel",x,y,ineueff,p);
    RooAddPdf final_pdf("final_pdf","final pdf",RooArgList(totdecaymodel,bkgmodelpos),RooArgList(nsig,nbkg));


    //!******************************************
    //! Toy MC dataset (just for testing)
    //! *****************************************
//    RooDataSet* datamc = final_pdf.generate(RooArgSet(x,y),10000) ;
//    datamc->Print() ;

    //!******************************************
    //! Prepare data set for fitting forward correlated data
    //! *****************************************
    RooDataSet* data=new RooDataSet("data","data",RooArgSet(x,y),Import(*tree)) ;
    data->Print() ;

    //!******************************************
    //! Perform the fit
    //! *****************************************
    RooFitResult* fitres;

    if (ffitopt==0)
        fitres=final_pdf.fitTo(*data,NumCPU(ncpu),Save(kTRUE));
    else
        fitres=final_pdf.fitTo(*data,ExternalConstraints(externalconstrains),NumCPU(ncpu),Save(kTRUE));


    //!******************************************
    //! Outputs
    //! *****************************************
    TFile* fout=new TFile(outputFileName,"recreate");
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
    c1->Write();
    fout->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
