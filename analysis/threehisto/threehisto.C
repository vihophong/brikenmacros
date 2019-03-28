#include "TH1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TFile.h"
#include <fstream>

void threehisto(TH1F* h1,TH1F* h2,TH1F* h3)
{
    TCanvas * c3 = new TCanvas("Title","Title",
                               10,10,800,800);
    c3->Divide(1,3,0,0);

    c3->cd(1);
    gPad->SetBottomMargin(0.001);
    gPad->SetTopMargin(0.01);
    gPad->SetRightMargin(0.01);
    h1->Draw("P0");

    c3->cd(2);
    gPad->SetBottomMargin(0.001);
    gPad->SetTopMargin(0.001);
    gPad->SetRightMargin(0.01);
    h2->Draw("P0");

    c3->cd(3);
    gPad->SetBottomMargin(0.1);
    gPad->SetTopMargin(0.001);
    gPad->SetRightMargin(0.01);
    h3->Draw("P0");
}
