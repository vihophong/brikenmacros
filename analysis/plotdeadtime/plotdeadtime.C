#include <iostream>
#include <iomanip>
#include <string>
#include <signal.h>
#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TCutG.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TStopwatch.h"
#include "TClonesArray.h"
#include "TVectorD.h"
#include "TSystem.h"
#include <fstream>

//#include "DataStructNew.h"
//#include "DataStructNewLinkDef.h"
#define MaxNRI 1000

void getdeadtime(char* listfile,Double_t* deadtimein,TH2F* h2){
    for (int i=0;i<6;i++) deadtimein[i]=0.;
    std::ifstream ifs(listfile);
    string filelist[1000];
    Int_t nfiles=0;
    while (!ifs.eof()){
        ifs>>filelist[nfiles];
        //cout<<filelist[nfiles]<<endl;
        nfiles++;
    }
    nfiles=nfiles-1;
    //cout<<"There are "<<nfiles<<" files in total!"<<endl;
    for (Int_t i=0;i<nfiles;i++){
        TFile *_file0 = TFile::Open(filelist[i].c_str());
        TVectorD* deadtimecontainer;
        deadtimecontainer=(TVectorD*) gDirectory->Get(Form("deadtime"));
        Double_t* deadtimearray=deadtimecontainer->GetMatrixArray();
        for (Int_t i=0;i<6;i++)
            deadtimein[i]+=deadtimearray[i];
        h2->Fill(i,(Double_t)deadtimearray[0]/(Double_t)deadtimearray[1]*100);
        _file0->Close();
    }
}


void chaineach(TChain* ch, char* listfile,char* pid){
  std::ifstream ifs(listfile);
  string filelist[1000];
  Int_t nfiles=0;
  while (!ifs.eof()){
      ifs>>filelist[nfiles];
      //cout<<filelist[nfiles]<<endl;
      nfiles++;
  }
  nfiles=nfiles-1;
  //cout<<"There are "<<nfiles<<" files in total!"<<endl;
  for (Int_t i=0;i<nfiles;i++){
      char tempchar2[1000];
      sprintf(tempchar2,"%s/tree%s",filelist[i].c_str(),pid);
      ch->Add(tempchar2);
  }
}

