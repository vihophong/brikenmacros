#include "TChain.h"
#include "TLatex.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TString.h"

void chainallimplant(char* listfile){
  char pid[500];
  sprintf(pid,"implant");
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
void chaineachimplant(char* listfile,char* pid){
  char tempchar1[1000];
  sprintf(tempchar1,"treeimp%s",pid);
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
      sprintf(tempchar2,"%s/treeimp%s",filelist[i].c_str(),pid);
      ch->Add(tempchar2);
  }
}

void chainallneutron(char* listfile){
  char pid[500];
  sprintf(pid,"neutron");
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

void chainimplant(char* listfile,char* pidfile){
    std::ifstream ifspid(pidfile);
    Int_t nri=0;
    Int_t nbinszet,nbinsaoq;
    Double_t zetrange[2];
    Double_t aoqrange[2];
    ifspid>>nbinsaoq>>aoqrange[0]>>aoqrange[1]>>nbinszet>>zetrange[0]>>zetrange[1];
    ifspid>>nri;

    Int_t enablepid[nri];
    Int_t enablepid2[nri];

    TString nameri[nri];
    TString latexnametri[nri];
    Double_t parmsri[nri][7];
    TString tempriname,tempria;
    TCutG* cutg[nri];
    TLatex* pidtag[nri];

    Double_t halflife[nri];

    for (Int_t i=0;i<nri;i++){
        ifspid>>enablepid[i]>>enablepid2[i]>>tempria>>tempriname>>halflife[i];
        for(Int_t j=0;j<7;j++) ifspid>>parmsri[i][j];
        nameri[i]=tempriname+tempria;
        latexnametri[i]=TString("^{")+tempria+TString("}"+tempriname);
        //cout<<nameri[i]<<"\t"<<latexnametri[i]<<"\t"<<halflife[i];
        for(Int_t j=0;j<7;j++) cout<<"\t"<<parmsri[i][j];
        //cout<<endl;
        pidtag[i]=new TLatex(parmsri[i][0],parmsri[i][1]+0.2,latexnametri[i]);
        pidtag[i]->SetTextSize(0.025);
        pidtag[i]->SetTextColor(2);
        chaineachimplant(listfile,(char*)nameri[i].Data());
    }
    chaineachimplant(listfile,(char*)TString("").Data());
}
void chaineach(char* listfile,char* pid){
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

void chain(char* listfile,char* pidfile){
    std::ifstream ifspid(pidfile);
    Int_t nri=0;
    Int_t nbinszet,nbinsaoq;
    Double_t zetrange[2];
    Double_t aoqrange[2];
    ifspid>>nbinsaoq>>aoqrange[0]>>aoqrange[1]>>nbinszet>>zetrange[0]>>zetrange[1];
    ifspid>>nri;

    Int_t enablepid[nri];
    Int_t enablepid2[nri];

    TString nameri[nri];
    TString latexnametri[nri];
    Double_t parmsri[nri][7];
    TString tempriname,tempria;
    TCutG* cutg[nri];
    TLatex* pidtag[nri];

    Double_t halflife[nri];

    for (Int_t i=0;i<nri;i++){
        ifspid>>enablepid[i]>>enablepid2[i]>>tempria>>tempriname>>halflife[i];
        for(Int_t j=0;j<7;j++) ifspid>>parmsri[i][j];
        nameri[i]=tempriname+tempria;
        latexnametri[i]=TString("^{")+tempria+TString("}"+tempriname);
        //cout<<nameri[i]<<"\t"<<latexnametri[i]<<"\t"<<halflife[i];
        for(Int_t j=0;j<7;j++) cout<<"\t"<<parmsri[i][j];
        //cout<<endl;
        pidtag[i]=new TLatex(parmsri[i][0],parmsri[i][1]+0.2,latexnametri[i]);
        pidtag[i]->SetTextSize(0.025);
        pidtag[i]->SetTextColor(2);
        chaineach(listfile,(char*)nameri[i].Data());
    }
    chaineach(listfile,(char*)TString("").Data());
}

void plotPID(char* listfile,char* pidfile){;
    std::ifstream ifspid(pidfile);
    Int_t nri=0;
    Int_t nbinszet,nbinsaoq;
    Double_t zetrange[2];
    Double_t aoqrange[2];
    ifspid>>nbinsaoq>>aoqrange[0]>>aoqrange[1]>>nbinszet>>zetrange[0]>>zetrange[1];
    ifspid>>nri;

    Int_t enablepid[nri];
    Int_t enablepid2[nri];

    TString nameri[nri];
    TString latexnametri[nri];
    Double_t parmsri[nri][7];
    TString tempriname,tempria;
    TCutG* cutg[nri];
    TLatex* pidtag[nri];

    Double_t halflife[nri];

    for (Int_t i=0;i<nri;i++){
        ifspid>>enablepid[i]>>enablepid2[i]>>tempria>>tempriname>>halflife[i];
        for(Int_t j=0;j<7;j++) ifspid>>parmsri[i][j];
        nameri[i]=tempriname+tempria;
        latexnametri[i]=TString("^{")+tempria+TString("}"+tempriname);
        cout<<nameri[i]<<"\t"<<latexnametri[i]<<"\t"<<halflife[i];
        for(Int_t j=0;j<7;j++) cout<<"\t"<<parmsri[i][j];
        cout<<endl;

        pidtag[i]=new TLatex(parmsri[i][0],parmsri[i][1]+0.2,latexnametri[i]);
        pidtag[i]->SetTextSize(0.025);
        pidtag[i]->SetTextColor(2);
    }
    std::ifstream ifs(listfile);
    string filelist[1000];
    Int_t nfiles=0;
    while (!ifs.eof()){
        ifs>>filelist[nfiles];
        cout<<filelist[nfiles]<<endl;
        nfiles++;
    }
/*
                                           Int_t ncutpts=20;// number of cut points
                                           for (Int_t i=0;i<nri;i++){
                                               cutg[i]=new TCutG(nameri[i],ncutpts);
                                               cutg[i]->SetTitle(nameri[i]);
                                               cutg[i]->SetVarX("decayaoq");
                                               cutg[i]->SetVarY("decay.zet");

                                               Double_t theta=0;
                                               Double_t x1=parmsri[i][0];
                                               Double_t y1=parmsri[i][1];
                                               Double_t r1=parmsri[i][2];
                                               Double_t r2=parmsri[i][3];

                                               Double_t phi1 = TMath::Min(0,360);
                                               Double_t phi2 = TMath::Max(0,360);
                                               Double_t kPI = 3.14159265358979323846;
                                               Double_t angle,dx,dy;

                                               Double_t x[ncutpts], y[ncutpts];
                                               Double_t dphi = (phi2-phi1)*kPI/(180*ncutpts);
                                               Double_t ct   = TMath::Cos(kPI*theta/180);
                                               Double_t st   = TMath::Sin(kPI*theta/180);

                                               for (Int_t j=0;j<ncutpts;j++) {
                                                  angle = phi1*kPI/180 + Double_t(j)*dphi;
                                                  dx    = r1*TMath::Cos(angle);
                                                  dy    = r2*TMath::Sin(angle);
                                                  x[j]  = x1 + dx*ct - dy*st;
                                                  y[j]  = y1 + dx*st + dy*ct;
                                                  //cout<<x[j]<<"<"<<y[j]<<endl;
                                                  cutg[i]->SetPoint(j,x[j],y[j]);
                                               }
                                               cutg[i]->SetPoint(ncutpts,x[0],y[0]);
                                           }
*/

    TFile* f0=new TFile(filelist[0].c_str());
    for (Int_t i=0;i<nri;i++){
         cutg[i]=(TCutG* )f0->Get(nameri[i]);
    }

    TChain* ch = new TChain("tree");
    nfiles=nfiles-1;
    cout<<"There are "<<nfiles<<" files in total!"<<endl;

    for (Int_t i=0;i<nfiles;i++){
        char tempchar2[1000];
        sprintf(tempchar2,"%s/tree",filelist[i].c_str());
        ch->Add(tempchar2);
    }
    TCanvas* c1=new TCanvas("pid","pid",900,700);
    ch->Draw(Form("zet:aoq>>h1(%d,%f,%f,%d,%f,%f)",nbinsaoq,aoqrange[0],aoqrange[1],nbinszet,zetrange[0],zetrange[1]),"zet>0&&aoq>0&&ion_z!=5&&ion_z>=0","colz");
    for (Int_t i=0;i<nri;i++){
         cutg[i]->Draw("same");
         pidtag[i]->Draw("same");
    }
}

void plotPIDimp(char* listfile,char* pidfile){;
    std::ifstream ifspid(pidfile);
    Int_t nri=0;
    Int_t nbinszet,nbinsaoq;
    Double_t zetrange[2];
    Double_t aoqrange[2];
    ifspid>>nbinsaoq>>aoqrange[0]>>aoqrange[1]>>nbinszet>>zetrange[0]>>zetrange[1];
    ifspid>>nri;

    Int_t enablepid[nri];
    Int_t enablepid2[nri];

    TString nameri[nri];
    TString latexnametri[nri];
    Double_t parmsri[nri][7];
    TString tempriname,tempria;
    TCutG* cutg[nri];
    TLatex* pidtag[nri];

    Double_t halflife[nri];

    for (Int_t i=0;i<nri;i++){
        ifspid>>enablepid[i]>>enablepid2[i]>>tempria>>tempriname>>halflife[i];
        for(Int_t j=0;j<7;j++) ifspid>>parmsri[i][j];
        nameri[i]=tempriname+tempria;
        latexnametri[i]=TString("^{")+tempria+TString("}"+tempriname);
        cout<<nameri[i]<<"\t"<<latexnametri[i]<<"\t"<<halflife[i];
        for(Int_t j=0;j<7;j++) cout<<"\t"<<parmsri[i][j];
        cout<<endl;

        pidtag[i]=new TLatex(parmsri[i][0],parmsri[i][1]+0.2,latexnametri[i]);
        pidtag[i]->SetTextSize(0.025);
        pidtag[i]->SetTextColor(2);
    }
    std::ifstream ifs(listfile);
    string filelist[1000];
    Int_t nfiles=0;
    while (!ifs.eof()){
        ifs>>filelist[nfiles];
        cout<<filelist[nfiles]<<endl;
        nfiles++;
    }
/*
                                           Int_t ncutpts=20;// number of cut points
                                           for (Int_t i=0;i<nri;i++){
                                               cutg[i]=new TCutG(nameri[i],ncutpts);
                                               cutg[i]->SetTitle(nameri[i]);
                                               cutg[i]->SetVarX("decayaoq");
                                               cutg[i]->SetVarY("decay.zet");

                                               Double_t theta=0;
                                               Double_t x1=parmsri[i][0];
                                               Double_t y1=parmsri[i][1];
                                               Double_t r1=parmsri[i][2];
                                               Double_t r2=parmsri[i][3];

                                               Double_t phi1 = TMath::Min(0,360);
                                               Double_t phi2 = TMath::Max(0,360);
                                               Double_t kPI = 3.14159265358979323846;
                                               Double_t angle,dx,dy;

                                               Double_t x[ncutpts], y[ncutpts];
                                               Double_t dphi = (phi2-phi1)*kPI/(180*ncutpts);
                                               Double_t ct   = TMath::Cos(kPI*theta/180);
                                               Double_t st   = TMath::Sin(kPI*theta/180);

                                               for (Int_t j=0;j<ncutpts;j++) {
                                                  angle = phi1*kPI/180 + Double_t(j)*dphi;
                                                  dx    = r1*TMath::Cos(angle);
                                                  dy    = r2*TMath::Sin(angle);
                                                  x[j]  = x1 + dx*ct - dy*st;
                                                  y[j]  = y1 + dx*st + dy*ct;
                                                  //cout<<x[j]<<"<"<<y[j]<<endl;
                                                  cutg[i]->SetPoint(j,x[j],y[j]);
                                               }
                                               cutg[i]->SetPoint(ncutpts,x[0],y[0]);
                                           }
*/

    TFile* f0=new TFile(filelist[0].c_str());
    for (Int_t i=0;i<nri;i++){
         cutg[i]=(TCutG* )f0->Get(nameri[i]);
    }

    TChain* ch = new TChain("treeimpall");
    nfiles=nfiles-1;
    cout<<"There are "<<nfiles<<" files in total!"<<endl;

    for (Int_t i=0;i<nfiles;i++){
        char tempchar2[1000];
        sprintf(tempchar2,"%s/treeimpall",filelist[i].c_str());
        ch->Add(tempchar2);
    }
    TCanvas* c1=new TCanvas("pid","pid",900,700);
    ch->Draw(Form("zet:aoq>>h1(%d,%f,%f,%d,%f,%f)",nbinsaoq,aoqrange[0],aoqrange[1],nbinszet,zetrange[0],zetrange[1]),"zet>0&&aoq>0&&fz!=5","colz");
    for (Int_t i=0;i<nri;i++){
         if (enablepid2[i]){
            pidtag[i]->SetTextSizePixels(25);
            cutg[i]->Draw("same");
            pidtag[i]->Draw("same");
         }
    }
}


void chaindeadtimetree(char* listfile){
  char tempchar1[1000];
  sprintf(tempchar1,"treedeadtime");
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
      sprintf(tempchar2,"%s/treedeadtime",filelist[i].c_str());
      ch->Add(tempchar2);
  }
}


void palette()
{
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.05);
  gStyle->SetPadRightMargin(0.15);
  TH2D *h = new TH2D("h", "h", 100, -5., 5., 100, -5., 5.);
  for(unsigned int i=0; i<10000000; i++){
    h->Fill(gRandom->Gaus(0., 1.), gRandom->Gaus(0., 1.));
  }
  TCanvas *c1 = new TCanvas("c1", "c", 1000,0,1200, 900+24);

  h->Draw("colz");
  /*
  Int_t col[12] = {1,51,52,53,56,81,93,99,100,103,106,110};
  TH2D *hh[12];

  for(int i=0;i<12;i++){
    c1->cd(i+1);
    gPad->SetLogz(1);
    if (i!=0) gStyle->SetPalette(col[i]);

    hh[i] = (TH2D*)h->Clone();
    if(i==0) hh[i]->SetTitle(Form("Default"));
    else     hh[i]->SetTitle(Form("SetPalette(%d)",col[i]));
    hh[i]->Draw("colz");
    gPad->Update();
  }
  c1->SaveAs("palette_log.png");
  */
}




