#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TCanvas.h>
#include <iostream>
#include <TGraph.h>
#include <TLine.h>
#include <TSpline.h>
#include <TMarker.h>

const int nmaxfiles=100;

void plotgeant4RvsE(char* dir,Int_t nfiles,Int_t trialsep=11)
{
    TFile *fileso[nmaxfiles];
    TTree *treeo[nmaxfiles];
    TH1F *hists[nmaxfiles];
    TH1F *hists2[nmaxfiles];
    Double_t arrsep[]={8.39295,8.58278,10.3195,10.6853,10.8032,12.6495,13.3655,14.0398,14.5,15.0877,16.7466,17.9478,18.2995,19.8366,20.9691,21.6478,22.7995,24.8015,26.2727,28.0596,29.0693,30.4814,31.2728,31.4496,34.2561,35.3504,35.8773,35.9917};
    Int_t narrsep=28;
    TGraph* gr[30];

    Double_t xarr[nmaxfiles];
    Double_t yarr[30][nmaxfiles];

    Double_t yskewarr[nmaxfiles];
    Double_t ymeanarr[nmaxfiles];
    Double_t ymeanandskewarr[nmaxfiles];

    //loop all files
    char tempname[500];
    for (Int_t i=0;i<nfiles;i++){
        //! get ratio parameter

        sprintf(tempname,"%stempHist%d.root",dir,i);
        fileso[i]=TFile::Open(tempname);
        treeo[i]=(TTree*)fileso[i]->Get("BRIKEN");

        for (Int_t j=0;j<narrsep;j++){
            Double_t ninner=(Double_t)treeo[i]->Draw("",Form("sqrt(x2*x2+y2*y2)<%f&&id<140",arrsep[j]),"goff");
            Double_t nouter=(Double_t)treeo[i]->Draw("",Form("sqrt(x2*x2+y2*y2)>=%f&&id<140",arrsep[j]),"goff");
            Double_t rparm=nouter/ninner;
            yarr[j][i]=rparm;
        }



        //! get original neutron energy
        treeo[i]->Draw(Form("partEnergy>>hpe%d(10000,0,10)",i),"","goff");
        sprintf(tempname,"hpe%i",i);
        hists[i]=(TH1F*)gDirectory->Get(tempname);
        Double_t partE=hists[i]->GetMean();
        xarr[i]=partE;

        //! get hit distribution histogram
        treeo[i]->Draw(Form("sqrt(x2*x2+y2*y2)>>hdist%d(200,0,50)",i),"","goff");
        sprintf(tempname,"hdist%i",i);
        hists2[i]=(TH1F*)gDirectory->Get(tempname);
        yskewarr[i]=hists2[i]->GetSkewness();
        ymeanarr[i]=hists2[i]->GetMean();
        ymeanandskewarr[i]=yskewarr[i]*ymeanarr[i];
    }

    for (Int_t i=0;i<narrsep;i++){
        gr[i]=new TGraph(nfiles,xarr,yarr[i]);
        gr[i]->SetName(Form("grr%d",i));
        gr[i]->SetMarkerStyle(20);
        gr[i]->SetLineColor(i+1);
    }

    TCanvas* c1=new TCanvas("c1","c1",900,700);
    treeo[3]->Draw("sqrt(x2*x2+y2*y2)>>h1(500,0,50)",Form("sqrt(x2*x2+y2*y2)<%f&&id<140",arrsep[trialsep]),"");
    treeo[3]->Draw("sqrt(x2*x2+y2*y2)>>h2(500,0,50)",Form("sqrt(x2*x2+y2*y2)>%f&&id<140",arrsep[trialsep]),"same");
    TH1F* h2=(TH1F*)gDirectory->Get("h2");
    h2->SetLineColor(2);
    h2->Draw("same");

    TCanvas* c2=new TCanvas("c2","c2",900,700);

    gr[0]->Draw("APL");
    gr[0]->GetYaxis()->SetRangeUser(0,26);
    for (Int_t i=1;i<narrsep;i++) gr[i]->Draw("PL same");


    TCanvas* c3=new TCanvas("c3","c3",900,700);
    TGraph* gr2[30];
    Double_t yarr2[30][nmaxfiles];


    //! normalized to first point
    for (Int_t i=0;i<narrsep;i++){
        yarr2[i][0]=0.;
        for (Int_t j=1;j<nfiles;j++){
            yarr2[i][j]=yarr[i][j]-yarr[i][0];
        }
    }

    for (Int_t i=0;i<narrsep;i++){
        gr2[i]=new TGraph(nfiles,xarr,yarr2[i]);
        gr2[i]->SetName(Form("grr2%d",i));
        gr2[i]->SetMarkerStyle(20);
        gr2[i]->SetLineColor(i+1);
    }

    gr2[0]->Draw("APL");
    gr2[0]->GetYaxis()->SetRangeUser(0,20);
    for (Int_t i=1;i<narrsep;i++) gr2[i]->Draw("PL same");

    TCanvas* c4=new TCanvas("c4","c4",900,700);
    c4->Divide(2,2);
    TGraph* grskew=new TGraph(nfiles,xarr,yskewarr);
    TGraph* grmean=new TGraph(nfiles,xarr,ymeanarr);
    TGraph* grmeanskew=new TGraph(nfiles,xarr,ymeanandskewarr);
    c4->cd(1);
    grskew->SetMarkerStyle(20);
    grskew->Draw("APL");
    c4->cd(2);
    grmean->SetMarkerStyle(20);
    grmean->Draw("APL");
    c4->cd(3);
    grmeanskew->SetMarkerStyle(20);
    grmeanskew->Draw("APL");
}

void calefficiency(char* dir,Int_t nfiles,Double_t inputr,Double_t sep=14.5,Double_t ntotal=1000000)
{
    TFile *fileso[nmaxfiles];
    TTree *treeo[nmaxfiles];
    TH1F *hists[nmaxfiles];

    Double_t xarr[nmaxfiles];
    Double_t yarr[nmaxfiles];
    Double_t effarr[nmaxfiles];

    //loop all files
    char tempname[500];
    for (Int_t i=0;i<nfiles;i++){
        //! get ratio parameter
        sprintf(tempname,"%stempHist%d.root",dir,i);
        fileso[i]=TFile::Open(tempname);
        treeo[i]=(TTree*)fileso[i]->Get("BRIKEN");

        Double_t ninner=(Double_t)treeo[i]->Draw("",Form("sqrt(x2*x2+y2*y2)<%f&&id<140",sep),"goff");
        Double_t nouter=(Double_t)treeo[i]->Draw("",Form("sqrt(x2*x2+y2*y2)>=%f&&id<140",sep),"goff");
        Double_t rparm=nouter/ninner;
        yarr[i]=rparm;

        //! get original neutron energy
        Double_t ndetected=(Double_t) treeo[i]->Draw(Form("partEnergy>>hpe%d(10000,0,10)",i),"id<140","goff");

        effarr[i]=ndetected/ntotal;
        sprintf(tempname,"hpe%i",i);
        hists[i]=(TH1F*)gDirectory->Get(tempname);
        Double_t partE=hists[i]->GetMean();
        xarr[i]=partE;
    }
    TGraph* gr=new TGraph(nfiles,xarr,yarr);
    TGraph* gr2=new TGraph(nfiles,yarr,xarr);
    gr->SetMarkerStyle(20);
    gr->SetLineColor(2);
    gr->Draw("APL");


    TMarker *mr=new TMarker(gr2->Eval(inputr),inputr,21);
    mr->SetMarkerColor(3);
    mr->Draw();
    TLine* l1=new TLine(xarr[0],inputr,xarr[nfiles-1],inputr);
    l1->SetLineColor(4);
    l1->Draw();
    TLine* l2=new TLine(gr2->Eval(inputr),yarr[0],gr2->Eval(inputr),yarr[nfiles-1]);
    l2->SetLineColor(4);
    l2->Draw();

    TGraph* greff=new TGraph(nfiles,xarr,effarr);
    greff->SetMarkerStyle(20);
    greff->SetLineColor(3);
    greff->Draw("PL same");

    TMarker *mr2=new TMarker(gr2->Eval(inputr),greff->Eval(gr2->Eval(inputr)),21);
    mr2->SetMarkerColor(5);
    mr2->Draw();

    //! summary
    cout<<"input_ratio\tmean_energy(MeV)\tefficiency(%)"<<endl;
    cout<<inputr<<"\t"<<gr2->Eval(inputr)<<"\t"<<greff->Eval(gr2->Eval(inputr))*100<<endl;

    TFile*outfile=new TFile("outfilereffvse.root","recreate");
    outfile->cd();
    gr->SetName("rvse");
    gr2->SetName("evsr");
    greff->SetName("effvsr");
    gr->Write();
    gr2->Write();
    greff->Write();
    outfile->Close();
}

void calefficiencyfast(Double_t inputr,Double_t inputrminus,Double_t inputrplus)
{
    TFile* file1=TFile::Open("outfilereffvse.root");
    TGraph* gr=(TGraph*)file1->Get("rvse");
    TGraph* gr2=(TGraph*)file1->Get("evsr");
    TGraph* greff=(TGraph*)file1->Get("effvsr");


    gr->Draw("APL");
    greff->Draw("PL same");

    TLine* l1=new TLine(gr->GetX()[0],inputr,gr->GetX()[gr->GetN()-1],inputr);
    l1->SetLineColor(2);
    l1->SetLineWidth(2);
    l1->Draw();
    TLine* l2=new TLine(gr2->Eval(inputr),gr->GetY()[0],gr2->Eval(inputr),gr->GetY()[gr->GetN()-1]);
    l2->SetLineColor(2);
    l2->SetLineWidth(2);
    l2->Draw();

    TMarker *mr=new TMarker(gr2->Eval(inputr),inputr,20);
    mr->SetMarkerColor(3);
    mr->SetMarkerSize(2);
    mr->Draw();
    TMarker *mr2=new TMarker(gr2->Eval(inputr),greff->Eval(gr2->Eval(inputr)),20);
    mr2->SetMarkerColor(5);
    mr2->SetMarkerSize(2);
    mr2->Draw();



    //! plus
    TLine* l1plus=new TLine(gr->GetX()[0],inputrplus,gr->GetX()[gr->GetN()-1],inputrplus);
    l1plus->SetLineColor(4);
    l1plus->Draw();

    TLine* l2plus=new TLine(gr2->Eval(inputrplus),gr->GetY()[0],gr2->Eval(inputrplus),gr->GetY()[gr->GetN()-1]);
    l2plus->SetLineColor(4);
    l2plus->Draw();

    TMarker *mrplus=new TMarker(gr2->Eval(inputrplus),inputrplus,21);
    mrplus->SetMarkerColor(3);
    mrplus->Draw();
    TMarker *mrplus2=new TMarker(gr2->Eval(inputrplus),greff->Eval(gr2->Eval(inputrplus)),21);
    mrplus2->SetMarkerColor(5);
    mrplus2->Draw();

    //! minus
    TLine* l1minus=new TLine(gr->GetX()[0],inputrminus,gr->GetX()[gr->GetN()-1],inputrminus);
    l1minus->SetLineColor(4);
    l1minus->Draw();

    TLine* l2minus=new TLine(gr2->Eval(inputrminus),gr->GetY()[0],gr2->Eval(inputrminus),gr->GetY()[gr->GetN()-1]);
    l2minus->SetLineColor(4);
    l2minus->Draw();

    TMarker *mrminus=new TMarker(gr2->Eval(inputrminus),inputrminus,21);
    mrminus->SetMarkerColor(3);
    mrminus->Draw();
    TMarker *mrminus2=new TMarker(gr2->Eval(inputrminus),greff->Eval(gr2->Eval(inputrminus)),21);
    mrminus2->SetMarkerColor(5);
    mrminus2->Draw();




    //! summary
    cout<<"input_ratio\tmean_energy(MeV)\tefficiency(%)"<<endl;
    cout<<inputr<<"\t"<<gr2->Eval(inputr)<<"\t"<<greff->Eval(gr2->Eval(inputr))*100<<endl;
    cout<<inputrplus<<"\t"<<gr2->Eval(inputrplus)<<"\t"<<greff->Eval(gr2->Eval(inputrplus))*100<<endl;
    cout<<inputrminus<<"\t"<<gr2->Eval(inputrminus)<<"\t"<<greff->Eval(gr2->Eval(inputrminus))*100<<endl;

}
