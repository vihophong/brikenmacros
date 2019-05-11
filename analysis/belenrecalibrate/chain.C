#include "TChain.h"
#include "TLatex.h"
void chaingamma(char* listfile){
  char pid[500];
  sprintf(pid,"gamma");
  char tempchar1[1000];
  sprintf(tempchar1,"%s",pid);
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
      sprintf(tempchar2,"%s/%s",filelist[i].c_str(),pid);
      ch->Add(tempchar2);
  }
  ch->SetName("gmm");
}


void chainanc(char* listfile){
  char pid[500];
  sprintf(pid,"anc");
  char tempchar1[1000];
  sprintf(tempchar1,"%s",pid);
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
      sprintf(tempchar2,"%s/%s",filelist[i].c_str(),pid);
      ch->Add(tempchar2);
  }
}

void chainneu(char* listfile){
  char pid[500];
  sprintf(pid,"neutron");
  char tempchar1[1000];
  sprintf(tempchar1,"%s",pid);
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
      sprintf(tempchar2,"%s/%s",filelist[i].c_str(),pid);
      ch->Add(tempchar2);
  }
}
