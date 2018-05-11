#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <fstream>
#include <iostream>
void checkstat(char* infile,Double_t low, Double_t high)
{
    std::ifstream ifs(infile);
    Double_t thr[2000];
    Int_t nlines=0;
    while (!ifs.eof()){
        ifs>>thr[nlines];
        nlines++;
    }
    TH1F* h1=new TH1F("h1","h1",50,low,high);
    for (int i=0;i<nlines-1;i++){
        h1->Fill(thr[i]);
    }
    h1->Draw();
}

void checkstathl(char* infile,Double_t low, Double_t high)
{
    std::ifstream ifs(infile);
    Double_t thr[2000];
    Int_t nlines=0;
    while (!ifs.eof()){
        ifs>>thr[nlines];
        nlines++;
    }
    TH1F* h1=new TH1F("h1","h1",50,log(2)/high,log(2)/low);
    for (int i=0;i<nlines-1;i++){
        h1->Fill(log(2)/thr[i]);
    }
    h1->Draw();
}
