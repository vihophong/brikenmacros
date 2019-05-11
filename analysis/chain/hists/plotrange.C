#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TROOT.h"
void plotrange(Int_t beg,Int_t end)
{
    Int_t npad=end-beg+1;
    TCanvas* c1=new TCanvas("c1","c1",1000,700);
    c1->Divide(npad/2,2);

    TH1F* hist[1000];
    for (Int_t i=beg;i<end+1;i++){
        c1->cd(i-beg+1);
        hist[i]=(TH1F*) gDirectory->Get(Form("h%d",i));
        hist[i]->GetXaxis()->SetRangeUser(0,100);
        hist[i]->Draw();
        c1->cd(i-beg+1)->SetLogy();
    }
}
