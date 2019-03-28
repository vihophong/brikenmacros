#include "TSpectrum.h"
#include "TH1.h"
#include "TFile.h"
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
#include <algorithm>

#include <TROOT.h>

TCanvas* c1;
void calibspec(Int_t ch=8) {
   gROOT->Reset();
   gROOT->Clear();


   Double_t max=2000;
   Double_t min=500;

   Double_t source[10000];
   TFile *f = new TFile("spec.root");
   c1=new TCanvas("c1","c1",900,700);
   c1->Divide(2,1);
   c1->cd(1);
   TH1F *h=(TH1F*) f->Get(Form("h%d",ch));
   TH1F* d1=(TH1F*)h->Clone();
   d1->Reset();
   h->Draw("L");
   TSpectrum *s = new TSpectrum();
   Int_t k=0;
   Int_t minbin,maxbin;



   for (Int_t i = 0; i < h->GetNbinsX(); i++) {
       if (h->GetXaxis()->GetBinCenter(i+1)>min&&h->GetXaxis()->GetBinCenter(i+1)<max){
           source[k]=h->GetBinContent(i + 1);
           if (k==0) minbin=i+1;
           maxbin=k+1;
           k++;
       }
   }

   s->Background(source,h->GetNbinsX(),15,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder8,kTRUE,
   TSpectrum::kBackSmoothing5,kTRUE);

   for (Int_t i = 0; i < h->GetNbinsX(); i++) {
       if (i>=minbin&&i<=maxbin)
       d1->SetBinContent(i + 1,source[i-minbin]);
       else
      d1->SetBinContent(i + 1,h->GetBinContent(i+1));

   }

   d1->SetLineColor(kRed);
   d1->Draw("SAME L");
   TH1F* d2=(TH1F*)h->Clone();
   d2->Add(d1,-1);
   d2->SetLineColor(kGreen);
   d2->Draw("same");

   s->Search(d2);
   Double_t* peaks=s->GetPositionX();


   //! sorting
   std::vector <Double_t> mappeak;
   std::vector <Double_t>::iterator mappeak_it;
   for (Int_t i=0;i<s->GetNPeaks();i++){
       mappeak.push_back(peaks[i]);
   }
   std::sort (mappeak.begin(), mappeak.begin()+4);
   k=0;
   for (mappeak_it=mappeak.begin();mappeak_it!=mappeak.end();mappeak_it++){
       peaks[k]=mappeak[k];
       cout<<peaks[k]<<endl;
       k++;
   }

   Double_t defaultrange=150;
   Double_t defaultheight=500;
   Double_t defaultsigma=20;

   Double_t peaklib[]={482,555,976,1049};


   TString fdef("gaus(0)+gaus(3)");
   TF1* ffinal=new TF1("ffinal",(char*)fdef.Data(),peaks[0]-defaultrange,peaks[1]+defaultrange);

   ffinal->SetParameter(0,defaultheight);// peak heigh
   ffinal->SetParameter(1,peaks[0]);// centroid
   ffinal->SetParameter(2,defaultsigma);// sigma
   ffinal->SetParameter(3,defaultheight);// peak heigh
   ffinal->SetParameter(4,peaks[1]);// centroid
   ffinal->SetParameter(5,defaultsigma);// sigma

   d2->Fit("ffinal","LQR","");
   for (Int_t i=0;i<6;i++) {
       cout<<i<<"\t"<<ffinal->GetParameter(i)<<endl;
   }

   TF1* ffinal2=new TF1("ffinal2",(char*)fdef.Data(),peaks[2]-defaultrange,peaks[3]+defaultrange);

   ffinal2->SetParameter(0,defaultheight);// peak heigh
   ffinal2->SetParameter(1,peaks[2]);// centroid
   ffinal2->SetParameter(2,defaultsigma);// sigma
   ffinal2->SetParameter(3,defaultheight);// peak heigh
   ffinal2->SetParameter(4,peaks[3]);// centroid
   ffinal2->SetParameter(5,defaultsigma);// sigma
   d2->Fit("ffinal2","LQR","");
   for (Int_t i=0;i<6;i++) {
       cout<<i<<"\t"<<ffinal->GetParameter(i)<<endl;
   }
   ffinal->Draw("same");

   c1->cd(2);
   TGraph *grcal=new TGraph(4,peaks,peaklib);
    grcal->Draw("AP*");
    grcal->Fit("pol1");
    Double_t offset=grcal->GetFunction("pol1")->GetParameter(0);
    Double_t gain=grcal->GetFunction("pol1")->GetParameter(1);

    cout<<"1\t"<<ffinal->GetParameter(2)*gain<<endl;
    cout<<"2\t"<<ffinal->GetParameter(5)*gain<<endl;
    cout<<"3\t"<<ffinal2->GetParameter(2)*gain<<endl;
    cout<<"4\t"<<ffinal2->GetParameter(5)*gain<<endl;

    //f->Close();


}
