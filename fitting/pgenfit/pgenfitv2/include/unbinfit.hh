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

#include "decaypath.hh"

class unbinfit
{
  public:
    unbinfit();
    virtual ~unbinfit();
    void Init(char* inputParms, char* inputData);
    void Run(char *outputFileName);

 private:
    char* finputParms;
    char* finputData;
    Long64_t fnentrieslimit;
    Double_t p_deadtime;
    Double_t p_timerange;

    Int_t ncpu;

    decaypath* path;

    TTree* tree;
};

#endif

