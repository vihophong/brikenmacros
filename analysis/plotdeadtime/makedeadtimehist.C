#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TVectorD.h"
#include <iostream>
using namespace std;

void makedeadtimehist(char* listfile){
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

  TCanvas* c1=new TCanvas("c1","c1",700,1400);

  TFile *fileso[nfiles];
  TH1F *histsnpulser[nfiles];
  //TH2F* hist2dpulser=new TH2F("h2dtpulser","h2dtpulser",nfiles,0,nfiles,2000,0,10);

  TH1F *hists[nfiles];
  TH2F* hist2d=new TH2F("h2","h2",nfiles,0,nfiles,2000,0,20);

  TH1F *hists2[nfiles];
  TH2F* hist2d2=new TH2F("h22","h22",nfiles,0,nfiles,2000,0,20);

  //TVectorD* deadtimecontainer[nfiles];

  for (Int_t i=0;i<nfiles;i++){
      fileso[i]=TFile::Open((char*)filelist[i].c_str());
      hists[i]=(TH1F*)gDirectory->Get("h1deadtime");
      hists2[i]=(TH1F*)gDirectory->Get("h1deadtime3");
      histsnpulser[i]=(TH1F*)gDirectory->Get("h1dtpulser");
      //deadtimecontainer[i]=(TVectorD*)gDirectory->Get("deadtime");
      //Double_t * arr=deadtimecontainer[i]->GetMatrixArray();
      //Double_t daqdeadtime=100*(1-(Double_t)histsnpulser[i]->GetEntries()/(arr[6]*10));
      //cout<<i<<"\t"<<filelist[i]<<daqdeadtime<<endl;
      //hist2dpulser->Fill(i,daqdeadtime);
      hists[i]->SetName(Form("hdeadtime%d",i));
      hists2[i]->SetName(Form("hdeadtime2%d",i));

      for (Int_t j=0;j<2000;j++){
          hist2d->SetBinContent(i+1,j,hists[i]->GetBinContent(j+1));
          hist2d2->SetBinContent(i+1,j,hists2[i]->GetBinContent(j+1));
      }
  }
  c1->Divide(1,2);
  //c1->cd(1);
  //hist2dpulser->Draw("colz");
  c1->cd(1);
  hist2d->SetTitle("Deadtime using pulser and start-stop time");
  hist2d->Draw("colz");
  c1->cd(2);
  hist2d2->SetTitle("Deadtime using pulser and deadtimepulser");
  hist2d2->Draw("colz");
}


void makedeadtimehistcorr(char* listfile){
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

  TCanvas* c1=new TCanvas("c1","c1",700,1400);

  TFile *fileso[nfiles];
  //TH2F* hist2dpulser=new TH2F("h2dtpulser","h2dtpulser",nfiles,0,nfiles,2000,0,10);

  TH1F *hists[nfiles];
  TH2F *histsD2[nfiles];
  TH2F* hist2d=new TH2F("h2","h2",nfiles,0,nfiles,2000,0,20);

  TH1F *hists2[nfiles];
  TH1F *hists3[nfiles];

  //TVectorD* deadtimecontainer[nfiles];

  for (Int_t i=0;i<nfiles;i++){
      fileso[i]=TFile::Open((char*)filelist[i].c_str());
      histsD2[i]=(TH2F*)gDirectory->Get("h2deadtime");

      hists[i]=(TH1F*) histsD2[i]->ProjectionY(Form("hdeadtime%d",i),2,24);
      hists2[i]=(TH1F*) histsD2[i]->ProjectionY(Form("prj2_%d",i),65,131);
      hists3[i]=(TH1F*) histsD2[i]->ProjectionY(Form("prj3_%d",i),133,141);

      hists[i]->Add(hists2[i]);
      hists[i]->Add(hists3[i]);

      for (Int_t j=0;j<2000;j++){
          hist2d->SetBinContent(i+1,j,hists[i]->GetBinContent(j+1));
      }
  }
  hist2d->SetTitle("Deadtime using pulser and start-stop time");
  hist2d->Draw("colz");
}
