#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TCutG.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLine.h"

#include "fitexp.C"

void timebkgsub(char* infile,char* histname, Double_t projwidth, Double_t centroid, Double_t projbkg)
{
    TFile* f1=TFile::Open(infile);
    TH2F* h2=(TH2F*)f1->Get(histname);
    Int_t bincentroid=h2->GetYaxis()->FindBin(centroid);
    Int_t binminus=bincentroid-(Int_t)round(projwidth/h2->GetYaxis()->GetBinWidth(1)/2);
    Int_t binplus=bincentroid+(Int_t)round(projwidth/h2->GetYaxis()->GetBinWidth(1)/2);
    Int_t nbinproj=binplus-binminus+1;
    if (nbinproj%2==1) {
        nbinproj--;
        binminus++;
    }

    Int_t binbkgminus=h2->GetYaxis()->FindBin(centroid-projbkg);
    Int_t binbkgplus=h2->GetYaxis()->FindBin(centroid+projbkg);




    TH1F* hprojpeak=(TH1F*)h2->ProjectionX("hcen",binminus,binplus);
    TH1F* hprojleftwing=(TH1F*)h2->ProjectionX("hleft",binbkgminus-nbinproj/2+1,binbkgminus);
    TH1F* hprojrightwing=(TH1F*)h2->ProjectionX("hright",binbkgplus,binbkgplus+nbinproj/2-1);

    TCanvas *c1 =new TCanvas("c1","c1",900,700);
    c1->cd();
    hprojpeak->Draw();
    hprojleftwing->SetLineColor(2);
    //hprojleftwing->Draw("same");
    hprojrightwing->SetLineColor(3);
    //hprojrightwing->Draw("same");


    TH1F* hprojbkg=(TH1F*) hprojleftwing->Clone();
    hprojbkg->Add(hprojrightwing);
    hprojbkg->SetName("hbkg");
    hprojbkg->SetLineColor(2);
    hprojbkg->Draw("same");
    /*
    TH1F* hsubtracted=(TH1F*)hprojpeak->Clone();
            hsubtracted->Add(hprojbkg,-1);
            hsubtracted->SetLineColor(1);
    hsubtracted->Draw("same");
    fitexp(hsubtracted,0.2,5);
    */



    TCanvas *c2 =new TCanvas("c2","c2",900,700);
    c2->cd();
    h2->Draw("colz");
    cout<<nbinproj<<endl;
    cout<<binminus<<"\t"<<binplus<<endl;
    cout<<binbkgminus-nbinproj/2+1<<"\t"<<binbkgminus<<endl;
    cout<<binbkgplus<<"\t"<<binbkgplus+nbinproj/2-1<<endl;

    TLine * linecenminus=new TLine(h2->GetXaxis()->GetXmin(),h2->GetYaxis()->GetBinCenter(binminus),h2->GetXaxis()->GetXmax(),h2->GetYaxis()->GetBinCenter(binminus));
    linecenminus->Draw();
    TLine * linecenplus=new TLine(h2->GetXaxis()->GetXmin(),h2->GetYaxis()->GetBinCenter(binplus),h2->GetXaxis()->GetXmax(),h2->GetYaxis()->GetBinCenter(binplus));
    linecenplus->Draw();


    TLine * linebkgminusleft=new TLine(h2->GetXaxis()->GetXmin(),h2->GetYaxis()->GetBinCenter(binbkgminus-nbinproj/2+1),h2->GetXaxis()->GetXmax(),h2->GetYaxis()->GetBinCenter(binbkgminus-nbinproj/2+1));
    linebkgminusleft->SetLineColor(2);
    linebkgminusleft->Draw();
    TLine * linebkgminusright=new TLine(h2->GetXaxis()->GetXmin(),h2->GetYaxis()->GetBinCenter(binbkgminus),h2->GetXaxis()->GetXmax(),h2->GetYaxis()->GetBinCenter(binbkgminus));
    linebkgminusright->SetLineColor(2);
    linebkgminusright->Draw();

    TLine * linebkgplusleft=new TLine(h2->GetXaxis()->GetXmin(),h2->GetYaxis()->GetBinCenter(binbkgplus),h2->GetXaxis()->GetXmax(),h2->GetYaxis()->GetBinCenter(binbkgplus));
    linebkgplusleft->SetLineColor(2);
    linebkgplusleft->Draw();
    TLine * linebkgplusright=new TLine(h2->GetXaxis()->GetXmin(),h2->GetYaxis()->GetBinCenter(binbkgplus+nbinproj/2-1),h2->GetXaxis()->GetXmax(),h2->GetYaxis()->GetBinCenter(binbkgplus+nbinproj/2-1));
    linebkgplusright->SetLineColor(2);
    linebkgplusright->Draw();


    TCanvas *c3 =new TCanvas("c3","c3",900,700);
    c3->cd();
    TH1F* hsubtracted=(TH1F*)hprojpeak->Clone();
    //hsubtracted->Add(hprojbkg,-1);
    hsubtracted->Draw();
    fitexp(hsubtracted,2,80);
    


}
