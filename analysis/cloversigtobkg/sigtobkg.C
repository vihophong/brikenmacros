//** FIND signal to background ratio for single peak****//

#include "TTree.h"
#include "TFile.h"

#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"

#include "TCut.h"
#include "TCutG.h"
#include "TLatex.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TCanvas.h"

TH1F* background_compton(TH1F* h,Double_t min,Double_t max) {

    Int_t binmin=h->GetXaxis()->FindBin(min);
    Int_t binmax=h->GetXaxis()->FindBin(max);
     Double_t * source = new Double_t[binmax-binmin+1];
     TH1F *d1 = new TH1F("d1","",h->GetXaxis()->GetNbins(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());

    h->Draw("L");
   TSpectrum *s = new TSpectrum();
   for (Int_t i = 0; i < binmax-binmin+1; i++) source[i]=h->GetBinContent(binmin + i);
   s->Background(source,binmax-binmin+1,10,1,3,kTRUE,5,kTRUE);
   for (Int_t i = 0; i < h->GetXaxis()->GetNbins(); i++) {
       if (i+1>=binmin&&i+1<=binmax)
       d1->SetBinContent(i + 1,source[i+1-binmin]);
   }
      d1->SetLineColor(kRed);
      d1->Draw("SAME L");
      return d1;
}
Double_t IntegralMod(TH1F* hin,Double_t min, Double_t max){

    Int_t binmin=hin->GetXaxis()->FindBin(min);
    Int_t binmax=hin->GetXaxis()->FindBin(max);
    Double_t mkx[2];
    Double_t mky[2];
    mkx[0]=min;
    mkx[1]=max;
    mky[0]=hin->GetBinContent(binmin);
    mky[1]=hin->GetBinContent(binmax);


    TGraph* gr2=new TGraph(2,mkx,mky);
    gr2->SetMarkerStyle(22);
    gr2->SetMarkerColor(4);
    gr2->SetLineColor(4);
    gr2->SetMarkerSize(1.5);
    gr2->SetLineWidth(2);
    gr2->Draw("P SAME");
    cout<<hin->Integral(binmin,binmax)<<endl;
    return (Double_t) hin->Integral(binmin,binmax);
}

void sigtobkg(TH1F* hin,Double_t bgrlow,Double_t bgrhi, Double_t pkrlow,Double_t pkrhi, Double_t cen_in=0,Double_t sig_in=0)
{
    TH1F* hbg=background_compton(hin,bgrlow,bgrhi);
    TH1F* hpk=(TH1F*)hin->Clone();
    hpk->Add(hbg,-1);
    hpk->Fit("gaus","QR+","goff",pkrlow,pkrhi);
    hin->Draw();
    hbg->Draw("same");
    hpk->GetFunction("gaus")->SetLineColor(3);
    hpk->SetLineColor(1);
    hpk->Draw("same");
    hpk->GetFunction("gaus")->Draw("same");
    Double_t cen=hpk->GetFunction("gaus")->GetParameter(1);
    Double_t sig=hpk->GetFunction("gaus")->GetParameter(2);
    cout<<cen<<"/"<<sig<<endl;
    if (!(cen_in==0&&sig_in==0)){
        cen=cen_in;
        sig=sig_in;
    }
    Double_t ntotal=IntegralMod(hin,cen-sig*5,cen+sig*5);
    Double_t nbg=IntegralMod(hbg,cen-sig*5,cen+sig*5);

    cout<<"\n\n*******\nE\tsignal\tbackground\tsig/bg"<<endl;
    cout<<cen<<"\t"<<ntotal-nbg<<"\t"<<nbg<<"\t"<<(ntotal-nbg)/nbg<<"\t"<<endl;
}
