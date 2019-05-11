#include "TChain.h"
#include "TLatex.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TString.h"

void chainpid(char* listfile){
  char pid[500];
  sprintf(pid,"");
  char tempchar1[1000];
  sprintf(tempchar1,"tree%s",pid);
  TChain* ch = new TChain(tempchar1);
  std::ifstream ifs(listfile);
  string filelist[1000];

  Int_t nfiles=0;
  while (!ifs.eof()){
      ifs>>filelist[nfiles];
      cout<<filelist[nfiles]<<endl;
      nfiles++;
  }
  nfiles=nfiles-1;
  cout<<"There are "<<nfiles<<" files in total!"<<endl;

  for (Int_t i=0;i<nfiles;i++){
      char tempchar2[1000];
      sprintf(tempchar2,"%s/tree%s",filelist[i].c_str(),pid);
      ch->Add(tempchar2);
  }
}
