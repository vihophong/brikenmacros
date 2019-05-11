#include "TChain.h"
#include "TLatex.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TString.h"
#include "TLine.h"
#include "TSpectrum.h"
#include "TMarker.h"
#include "TF1.h"
#include <fstream>

void docalib(char* infile)
{
    TFile* file1=TFile::Open(infile);


    //! get first group
    TH2F* hgroup[10];
    Int_t nbins[10];

    Double_t rlow[]={500000,300000,800000,800000,3000000};
    Double_t rhigh[]={1200000,1000000,2000000,2000000,9000000};


    TCanvas* c1=new TCanvas("c1","c1",900,700);
    c1->Divide(3,2);
    for (Int_t i=0;i<5;i++){
        hgroup[i]=(TH2F*)file1->Get(Form("h2dneugroup%d",i+1));

        c1->cd(i+1)->SetLogz();
        hgroup[i]->Draw("colz");
        nbins[i]= hgroup[i]->GetNbinsX();
        TH1F* hproj[50];
        for (Int_t j=0;j<nbins[i];j++){
            hproj[j]=(TH1F*) hgroup[i]->ProjectionY(Form("proj%d_%d",i,j),j+1,j+1);

            //! find maximum here
            Int_t binlow=hproj[j]->FindBin(rlow[i]);
            Int_t binhi=hproj[j]->FindBin(rhigh[i]);
            Int_t max=0;
            Double_t maxval=0;
            for (Int_t b=binlow;b<binhi;b++){
                if (hproj[j]->GetBinContent(b)>max) {
                    max=hproj[j]->GetBinContent(b);
                    maxval=b;
                }
            }
            maxval=hproj[j]->GetBinCenter(maxval);
            TMarker* mrk=new TMarker(hgroup[i]->GetXaxis()->GetBinCenter(j+1),maxval,20);
            mrk->SetMarkerColor(2);
            mrk->Draw();
            cout<<hgroup[i]->GetXaxis()->GetBinCenter(j+1)-0.5<<"\t"<<765./maxval<<endl;

        }

    }

}
