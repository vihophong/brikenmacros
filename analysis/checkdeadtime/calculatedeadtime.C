#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TString.h"
#include "TMath.h"
#include <fstream>

void calculatedeadtime(char* infile="out_deadtime_160221.txt")
{

    TGraph* gr1=new TGraph(infile,"%lg %*lg %*lg %lg %*lg %*lg %*lg %*lg");
    TGraph* gr2=new TGraph(infile,"%lg %*lg %*lg %*lg %lg %*lg %*lg %*lg");
    TGraph* gr3=new TGraph(infile,"%lg %*lg %*lg %*lg %*lg %lg %*lg %*lg");
    TGraph* gr4=new TGraph(infile,"%lg %*lg %*lg %*lg %*lg %*lg %lg %*lg");
    TGraph* gr5=new TGraph(infile,"%lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
    Double_t* gr1_y=gr1->GetY();
    Double_t* gr2_y=gr2->GetY();
    Double_t* gr3_y=gr3->GetY();
    Double_t* gr4_y=gr4->GetY();
    Double_t* gr5_y=gr5->GetY();

    Long64_t t_deadtime=0;
    Long64_t t_total=0;
    for (Int_t i=0;i<gr1->GetN();i++){
        t_deadtime+=gr1_y[i];
        t_total+=gr2_y[i];
    }
    cout<<(Double_t)t_deadtime/(Double_t)t_total*100<<"\t"<<TMath::Mean(gr3->GetN(),gr3_y)<<"\t"
       <<TMath::Mean(gr4->GetN(),gr4_y)
       <<"\t"<<TMath::Mean(gr5->GetN(),gr5_y)<<endl;
}
