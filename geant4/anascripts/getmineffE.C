#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <string>
#include <TLatex.h>
#include <TLine.h>
using namespace std;
void getmineffE(Double_t QbnMeV=1){

    TFile* halle=new TFile("runalle/effcurvealle.root");

    TGraph* grAllE=(TGraph*) halle->Get("grAllE");

    TCanvas * c1=new TCanvas("c1","c1",900,700);
    c1->cd();
    c1->SetLogx();
    grAllE->Draw("APL");


    Double_t* arrY=grAllE->GetY();
    Double_t maxY=0;
    for (Int_t i=0;i<grAllE->GetN();i++){
        if (arrY[i]>maxY) maxY=arrY[i];
    }
    TLine line;
    TLatex latex;

    line.DrawLine(grAllE->GetXaxis()->GetXmin(),maxY,grAllE->GetXaxis()->GetXmax(),maxY);
    latex.SetTextSize(0.03);
    latex.DrawLatex(grAllE->GetXaxis()->GetXmin()+0.1,maxY,Form("Eff_Max=%f %%",maxY));

    Double_t Eeval=grAllE->Eval(QbnMeV*1000);
    line.DrawLine(QbnMeV*1000,grAllE->GetYaxis()->GetXmin(),QbnMeV*1000,grAllE->GetYaxis()->GetXmax());

    line.DrawLine(grAllE->GetXaxis()->GetXmin(),Eeval,grAllE->GetXaxis()->GetXmax(),Eeval);
    latex.SetTextSize(0.03);
    latex.DrawLatex(grAllE->GetXaxis()->GetXmin()+100,Eeval,Form("Eff_Min=%f %%",Eeval));

    cout<<"Eff_min(%)\tEff_Max(%)"<<endl;
    cout<<Eeval<<"\t"<<maxY<<endl;
    cout<<"Effective_Eff_min(%)\tEffective_Eff_Max(%)"<<endl;
    cout<<(Eeval/100-Eeval/100*0.04)*100<<"\t"<<(maxY/100-maxY/100*0.04)*100<<endl;
}
