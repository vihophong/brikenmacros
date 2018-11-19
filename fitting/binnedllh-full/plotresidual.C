#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

void plotresidual(char* infilename,Double_t lower=0.05,Double_t upper=10)
{
    TFile* infile=TFile::Open(infilename);
    char tempchar1[1000];
    sprintf(tempchar1,"residual_decay");
    TH1F* hdecayres=(TH1F*) gDirectory->Get(tempchar1);

    sprintf(tempchar1,"residual_decay1neu");
    TH1F* hdecay1nres=(TH1F*) gDirectory->Get(tempchar1);

    sprintf(tempchar1,"residual_decay2neu");
    TH1F* hdecay2nres=(TH1F*) gDirectory->Get(tempchar1);


    TCanvas* c2=new TCanvas("c2","c2",900,1200);
    c2->Divide(1,3);
    c2->cd(1);
    hdecayres->GetXaxis()->SetRangeUser(lower,upper);
    hdecayres->GetXaxis()->SetTitle("Time (s)");
    hdecayres->GetYaxis()->SetTitle("Counts");
    hdecayres->GetXaxis()->SetTitleSize(0.05);
    hdecayres->GetYaxis()->SetTitleSize(0.05);
    hdecayres->SetMarkerStyle(20);
    hdecayres->SetMarkerSize(0.8);
    hdecayres->GetXaxis()->SetLabelSize(0.05);
    hdecayres->GetYaxis()->SetLabelSize(0.05);
    hdecayres->Draw("P0L");
    c2->cd(2);
    hdecay1nres->GetXaxis()->SetRangeUser(lower,upper);
    hdecay1nres->GetXaxis()->SetTitle("Time (s)");
    hdecay1nres->GetYaxis()->SetTitle("Counts");
    hdecay1nres->GetXaxis()->SetTitleSize(0.05);
    hdecay1nres->GetYaxis()->SetTitleSize(0.05);
    hdecay1nres->SetMarkerStyle(20);
    hdecay1nres->SetMarkerSize(0.8);
    hdecay1nres->GetXaxis()->SetLabelSize(0.05);
    hdecay1nres->GetYaxis()->SetLabelSize(0.05);
    hdecay1nres->Draw("P0L");
    c2->cd(3);
    hdecay2nres->GetXaxis()->SetRangeUser(lower,upper);
    hdecay2nres->GetXaxis()->SetTitle("Time (s)");
    hdecay2nres->GetYaxis()->SetTitle("Counts");
    hdecay2nres->GetXaxis()->SetTitleSize(0.05);
    hdecay2nres->GetYaxis()->SetTitleSize(0.05);
    hdecay2nres->SetMarkerStyle(20);
    hdecay2nres->SetMarkerSize(0.8);
    hdecay2nres->GetXaxis()->SetLabelSize(0.05);
    hdecay2nres->GetYaxis()->SetLabelSize(0.05);
    hdecay2nres->Draw("P0L");

    sprintf(tempchar1,"c1");
    TCanvas* cfit=(TCanvas*) gDirectory->Get(tempchar1);
    cfit->Draw();
}
