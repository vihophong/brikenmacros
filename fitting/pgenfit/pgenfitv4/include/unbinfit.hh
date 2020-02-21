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

//! \file unbinfit.hh
//! \brief Definition of the unbinfit class

#ifndef unbinfit_h
#define unbinfit_h 1

//! Unbinned fitting class.

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include "RooDataSet.h"

#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooPlot.h"

#include "RooAddPdf.h"
#include "TString.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFitResult.h"
#include "TMath.h"

#include "TStopwatch.h"
#include "TRandom3.h"

#include "common.hh"
#include "fitF.hh"

#include "decaypath.hh"

class unbinfit
{
  public:
    unbinfit();
    virtual ~unbinfit();
    void Init(char* inputParms, char* inputData);
    void SetParameters();

    void setOutputFile(char* outputData){sprintf(foutputData,"%s",outputData);}
    void setStartTime(double deadtime){p_deadtime=deadtime;}
    void Run();
    void generateRoofitEvaluate();

    void fitBackground();

    void initFitParameters();
    void setNormalFit();//decide parameter is fix or not
    void setExernalContrainFit();

    void prepareData();
    void prepareMonteCarloData(int nevents);
    void doFit();

    void plotResults();
    void writeResults();
    void writeOutputTree(){foutputtree->Write();}
    void closeOutputFile(){fout->Close();}

    //! stuff for MC set
    void setNumberOfMC(Int_t ntimes){fnMC=ntimes;}

    void setCentralParameters();
    void bookOutputTree();
    void getParameters();
    void printCurrentParameters();
    void setValParameters();
    void generateMC();
    void writeResultsMC();


 private:
    void setModel();
    char* finputParms;
    char* finputData;
    char foutputData[500];
    TFile* fout;
    Long64_t fnentrieslimit;
    Double_t p_deadtime;
    Double_t p_timerange;

    Int_t ncpu;

    Int_t ffitopt;

    decaypath* fdecaypath;

    //! 2 dimensions parameters
    RooRealVar *x;
    RooCategory *y;

    //! background parameters
    RooRealVar* xbkg;
    RooRealVar* bkg1nratio;
    RooRealVar* bkg2nratio;

    RooRealVar* slope1pos;
    RooRealVar* slope2pos;
    RooRealVar* slope3pos;

    fitFbkg* bkgmodelpos;

    //! fit paramters
    // Declare all parameters
    RooAbsReal* p[kmaxparms];//fdecaypath->getNMember()*5+4];
    RooRealVar* pvar[kmaxparms];//[fdecaypath->getNMember()*5+4];
    RooRealVar* nbkg;
    RooRealVar* nsig;

    //! constrains set
    RooArgSet *externalconstrains;


    //! constrains set
    fitF* totdecaymodel;
    RooAddPdf* final_pdf;


    //! data sets
    TTree* tree;
    RooDataSet* databkg; //background data
    RooDataSet* data;

    //! fit results
    RooFitResult* fitres;

    //! stuffs for MC generation
    TRandom3* rseed;
    Int_t fnMC;
    TTree* foutputtree;

    Double_t pVal[kmaxparms];
    Double_t pCentralVal[kmaxparms];
    Double_t pValError[kmaxparms];
    Int_t ispVary[kmaxparms];
    Int_t ipVal[kmaxparms];

    Double_t nbkgVal;
    Double_t nbkgCentralVal;
    Double_t nbkgError;
    Double_t bkg1nratioVal;
    Double_t bkg1nratioCentralVal;
    Double_t bkg1nratioError;
    Double_t bkg2nratioVal;
    Double_t bkg2nratioCentralVal;
    Double_t bkg2nratioError;
    Double_t slope1posVal;
    Double_t slope2posVal;
    Double_t slope3posVal;
    Double_t slope1posCentralVal;
    Double_t slope2posCentralVal;
    Double_t slope3posCentralVal;
    Double_t slope1posError;
    Double_t slope2posError;
    Double_t slope3posError;

    Int_t fitStatus;
    Int_t fitCovQual;
    Int_t fitNumInvalidNLL;
    Double_t fitEdm;
    Double_t fitMinNll;

    Double_t nsigVal;
    Double_t nsigCentralVal;
    Double_t nsigError;

    TStopwatch *fStopWatch;//for debuging
    Double_t fMCGenTime;
    Double_t fFitTime;
};

#endif

