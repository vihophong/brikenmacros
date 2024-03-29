#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TCanvas.h>
#include <iostream>
#include <TGraph.h>


void ploteffvse()
{
    std::ifstream inpf("outGrid.txt");
    Double_t temp[500][12];
    Int_t nlines=0;
    while (inpf.good()){
        for (Int_t i=0;i<12;i++){
            inpf>>temp[nlines][i];
            cout<<temp[nlines][i]<<"\t";
        }
        nlines++;
        cout<<endl;
    }
    cout<<"\n\n"<<endl;
    cout<<"\n\n"<<endl;
    TCanvas* c2=new TCanvas("c2","c2",900,700);
    c2->cd();
    c2->SetLogx();
    c2->SetGrid();
    Double_t enearr[]={0.1,1,10,100,500,1000,2000,3000,4000,5000};
    Double_t effarr[10];
    for (Int_t i=0;i<12;i++){
        cout<<temp[nlines-2][i]<<"\t";
        if (i>1){
            effarr[i-2]=temp[nlines-2][i];
        }
    }
    TGraph* gr=new TGraph(10,enearr,effarr);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(2);
    gr->SetLineColor(2);
    gr->SetMarkerSize(1.2);
    gr->SetLineWidth(2);
    gr->GetYaxis()->SetRangeUser(0,90);
    gr->Draw("APL");
    cout<<endl;
}

void makeoutputfiles(){
    TString* legentryname[11];
    Int_t colorcode[]={1,2,3,5,6,7,8,9,11,12,13};
    legentryname[0]=new TString("100 eV");
    legentryname[1]=new TString("1 keV");
    legentryname[2]=new TString("10 keV");
    legentryname[3]=new TString("100 keV");
    legentryname[4]=new TString("500 keV");
    legentryname[5]=new TString("1000 keV");
    legentryname[6]=new TString("2000 keV");
    legentryname[7]=new TString("3000 keV");
    legentryname[8]=new TString("4000 keV");
    legentryname[9]=new TString("5000 keV");
    legentryname[10]=new TString("252Cfsimulation");
    TString* entryname[11];
    entryname[0]=new TString("0_1keV");
    entryname[1]=new TString("1keV");
    entryname[2]=new TString("10keV");
    entryname[3]=new TString("100keV");
    entryname[4]=new TString("500keV");
    entryname[5]=new TString("1000keV");
    entryname[6]=new TString("2000keV");
    entryname[7]=new TString("3000keV");
    entryname[8]=new TString("4000keV");
    entryname[9]=new TString("5000keV");
    entryname[10]=new TString("252Cfsimulation");

    TFile *fileso[11];
    TTree *treeo[11];
    TH1F *hists[11];
    TH1F *hists2[11];

    TCanvas* c1=new TCanvas("c1","c1",900,700);
    c1->Divide(1,2);
    c1->cd(1);
    TLegend* leg = new TLegend(0.2, 0.2, .8, .8);
    TLegend* leg2 = new TLegend(0.2, 0.2, .8, .8);

    for (Int_t i=0;i<11;i++){
        char tempname[500];
        sprintf(tempname,"tempHist%d.root",i);
        fileso[i]=TFile::Open(tempname);
        treeo[i]=(TTree*)fileso[i]->Get("BRIKEN");
        treeo[i]->Draw(Form("x*x+y*y>>hcont%s(200,0,2000)",(char*)entryname[i]->Data()),"","goff");

        sprintf(tempname,"hcont%s",(char*)entryname[i]->Data());
        hists[i]=(TH1F*)gDirectory->Get(tempname);
        treeo[i]->Draw(Form("x2*x2+y2*y2>>hdisc%s(200,0,2000)",(char*)entryname[i]->Data()),"","goff");
        sprintf(tempname,"hdisc%s",(char*)entryname[i]->Data());
        hists2[i]=(TH1F*)gDirectory->Get(tempname);

        hists[i]->SetLineColor(colorcode[i]);
        hists[i]->SetLineWidth(2);
        hists2[i]->SetLineColor(colorcode[i]);
        hists2[i]->SetLineWidth(2);

        leg->AddEntry(hists[i],(char*)legentryname[i]->Data());
        leg2->AddEntry(hists2[i],(char*)legentryname[i]->Data());
    }
    for (Int_t i=0;i<11;i++){
        if (i==0) hists[i]->Draw();
        else hists[i]->Draw("same");
    }
    leg->Draw();

    c1->cd(2);
    for (Int_t i=0;i<11;i++){
        if (i==0) hists2[i]->Draw();
        else hists2[i]->Draw("same");
    }
    leg2->Draw();

    std::ifstream inpf("outGrid.txt");
    Double_t temp[500][12];
    Int_t nlines=0;
    while (inpf.good()){
        for (Int_t i=0;i<12;i++){
            inpf>>temp[nlines][i];
            cout<<temp[nlines][i]<<"\t";
        }
        nlines++;
        cout<<endl;
    }
    cout<<"\n\n"<<endl;
    TCanvas* c2=new TCanvas("c2","c2",900,700);
    c2->cd();
    c2->SetLogx();
    c2->SetGrid();
    Double_t enearr[]={0.1,1,10,100,500,1000,2000,3000,4000,5000};
    Double_t effarr[10];
    for (Int_t i=0;i<12;i++){
        cout<<temp[nlines-2][i]<<"\t";
        if (i>1){
            effarr[i-2]=temp[nlines-2][i];
        }
    }
    TGraph* gr=new TGraph(10,enearr,effarr);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(2);
    gr->SetLineColor(2);
    gr->SetMarkerSize(1.2);
    gr->SetLineWidth(2);
    gr->GetYaxis()->SetRangeUser(0,90);
    gr->Draw("APL");
    cout<<endl;

    TFile* ofile[10];
    for (Int_t i=0;i<11;i++){
        char tempname[500];
        sprintf(tempname,"hist/neuhit%s.root",(char*)entryname[i]->Data());
        ofile[i]=new TFile(tempname,"recreate");
        hists2[i]->SetName("h1");
        hists2[i]->Write();
        ofile[i]->Close();
    }
}
void plotdistrib(char* outfilename){

    TString* legentryname[11];
    Int_t colorcode[]={1,2,3,5,6,7,8,9,11,12,13};
    legentryname[0]=new TString("100 eV");
    legentryname[1]=new TString("1 keV");
    legentryname[2]=new TString("10 keV");
    legentryname[3]=new TString("100 keV");
    legentryname[4]=new TString("500 keV");
    legentryname[5]=new TString("1000 keV");
    legentryname[6]=new TString("2000 keV");
    legentryname[7]=new TString("3000 keV");
    legentryname[8]=new TString("4000 keV");
    legentryname[9]=new TString("5000 keV");
    legentryname[10]=new TString("252Cf");
    TString* entryname[11];
    entryname[0]=new TString("0_1keV");
    entryname[1]=new TString("1keV");
    entryname[2]=new TString("10keV");
    entryname[3]=new TString("100keV");
    entryname[4]=new TString("500keV");
    entryname[5]=new TString("1000keV");
    entryname[6]=new TString("2000keV");
    entryname[7]=new TString("3000keV");
    entryname[8]=new TString("4000keV");
    entryname[9]=new TString("5000keV");
    entryname[10]=new TString("252Cfsimulation");

    Bool_t flag[100];
    flag[0]=false;
    flag[1]=false;
    flag[2]=false;
    flag[3]=true;
    flag[4]=false;
    flag[5]=true;
    flag[6]=false;
    flag[7]=true;
    flag[8]=false;
    flag[9]=true;
    flag[10]=false;

    TFile *fileso[11];
    TTree *treeo[11];
    TH1F *hists[11];
    TH1F *hists2[11];

    TCanvas* c1=new TCanvas("c1","c1",900,700);
    c1->Divide(1,2);
    c1->cd(1);
    TLegend* leg = new TLegend(0.2, 0.2, .8, .8);
    TLegend* leg2 = new TLegend(0.2, 0.2, .8, .8);

    for (Int_t i=0;i<11;i++){
        char tempname[500];
        sprintf(tempname,"tempHist%d.root",i);
        fileso[i]=TFile::Open(tempname);
        treeo[i]=(TTree*)fileso[i]->Get("BRIKEN");
        treeo[i]->Draw(Form("x*x+y*y>>hcont%s(200,0,2000)",(char*)entryname[i]->Data()),"","goff");

        sprintf(tempname,"hcont%s",(char*)entryname[i]->Data());
        hists[i]=(TH1F*)gDirectory->Get(tempname);
        treeo[i]->Draw(Form("x2*x2+y2*y2>>hdisc%s(200,0,2000)",(char*)entryname[i]->Data()),"","goff");
        sprintf(tempname,"hdisc%s",(char*)entryname[i]->Data());
        hists2[i]=(TH1F*)gDirectory->Get(tempname);

        hists[i]->SetLineColor(colorcode[i]);
        hists[i]->SetLineWidth(2);
        hists2[i]->SetLineColor(colorcode[i]);
        hists2[i]->SetLineWidth(2);        

    }
    for (Int_t i=0;i<11;i++){
        leg->AddEntry(hists[i],(char*)legentryname[i]->Data());
        if (i==0) hists[i]->Draw();
        else hists[i]->Draw("same");
    }
    leg->Draw();

    c1->cd(2);
    for (Int_t i=0;i<11;i++){
        if (flag[i]){
            leg2->AddEntry(hists2[i],(char*)legentryname[i]->Data());
            if (i==0) hists2[i]->Draw();
            else hists2[i]->Draw("same");
        }
    }
    leg2->Draw();

    std::ifstream inpf("outGrid.txt");
    Double_t temp[500][12];
    Int_t nlines=0;
    while (inpf.good()){
        for (Int_t i=0;i<12;i++){
            inpf>>temp[nlines][i];
            cout<<temp[nlines][i]<<"\t";
        }
        nlines++;
        cout<<endl;
    }
    cout<<"\n\n"<<endl;
    TCanvas* c2=new TCanvas("c2","c2",900,700);
    c2->cd();
    c2->SetLogx();
    c2->SetGrid();
    Double_t enearr[]={0.1,1,10,100,500,1000,2000,3000,4000,5000};
    Double_t effarr[10];
    for (Int_t i=0;i<12;i++){
        cout<<temp[nlines-2][i]<<"\t";
        if (i>1){
            effarr[i-2]=temp[nlines-2][i];
        }
    }
    TGraph* gr=new TGraph(10,enearr,effarr);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(2);
    gr->SetLineColor(2);
    gr->SetMarkerSize(1.2);
    gr->SetLineWidth(2);
    gr->GetYaxis()->SetRangeUser(0,90);
    gr->Draw("APL");
    cout<<endl;

    TFile* ofile=new TFile(outfilename,"recreate");
    ofile->cd();
    c1->Write();
    c2->Write();
    gr->Write();
    for (Int_t i=0;i<11;i++){
        hists[i]->Write();
        hists2[i]->Write();
    }
    ofile->Close();
}



void plotneudists(char* outfilename)
{
    Int_t colorcode[]={1,2,3,4,5,6,7,9,11,12,13,14,15,16,17,18,19,20,21,22,23};
    TString ri[100];
    Bool_t flag[100];
    Int_t nri=6;
    //ri[0]="Sn134";
    //ri[1]="Sn135";
    //ri[2]="Sn136";
    ri[0]="500keV";
    ri[1]="1000keV";
    ri[2]="2000keV";
    ri[3]="3000keV";
    ri[4]="252Cfsimulation";
    ri[5]="252Cfexperiment";


    TString riname[100];
    riname[0]="500keV simulation";
    riname[1]="1000keV simulation";
    riname[2]="2000keV simulation";
    riname[3]="3000keV simulation";
    riname[4]="252Cf simulation";
    riname[5]="252Cf experiment";

    flag[0]=true;
    flag[1]=true;
    flag[2]=true;
    flag[3]=true;
    flag[4]=true;
    flag[5]=true;
    TFile *fileso[nri];
    TH1F *hists[nri];
    TGraph *grhists[nri];

    TCanvas* c1=new TCanvas("c1","c1",900,700);
    TLegend* leg = new TLegend(0.2, 0.2, .8, .8);


    Double_t bin8countsnorm=10000.;
    Double_t bin8counts;


    for (Int_t i=0;i<nri;i++){
        Double_t gxx[200];
        Double_t gyy[200];
        Int_t gnpoints=0;

        char tempchar1[1000];
        sprintf(tempchar1,"hist/neuhit%s.root",(char*)ri[i].Data());
        fileso[i]=TFile::Open(tempchar1);
        hists[i]=(TH1F*)gDirectory->Get("h1");
        hists[i]->SetName((char*)ri[i].Data());
        hists[i]->SetMarkerStyle(kFullCircle);
        //hists[i]->Scale(norm, "width");
        bin8counts=(Double_t)hists[i]->GetBinContent(8);
        hists[i]->Scale(bin8countsnorm/bin8counts);

        //hists[i]->SetLineWidth(0);

        for (Int_t j=0;j<200;j++){
            if (hists[i]->GetBinContent(j+1)>0){
                gxx[gnpoints]=hists[i]->GetXaxis()->GetBinCenter(j+1);
                gyy[gnpoints]=hists[i]->GetBinContent(j+1);
                gnpoints++;
            }
        }
        grhists[i]=new TGraph(gnpoints,gxx,gyy);
        grhists[i]->SetMarkerStyle(kFullCircle);
        grhists[i]->SetFillColor(0);
        if (i==4||i==5) {
            grhists[i]->SetLineWidth(2);
            grhists[i]->SetMarkerSize(2);

            grhists[i]->SetMarkerStyle(kFullSquare);
        }
    }
    Int_t ncolor=0;
    for (Int_t i=0;i<nri;i++){
        if (flag[i]){
            hists[i]->SetMarkerColor(colorcode[ncolor]);
            hists[i]->SetLineColor(colorcode[ncolor]);
            grhists[i]->SetMarkerColor(colorcode[ncolor]);
            grhists[i]->SetLineColor(colorcode[ncolor]);
            leg->AddEntry(grhists[i],(char*)riname[i].Data());
            if (ncolor==0) grhists[i]->Draw("APL");
            else grhists[i]->Draw("PLSAME");
            ncolor++;
        }
    }
    leg->Draw();
    TFile* ofile=new TFile(outfilename,"recreate");
    ofile->cd();
    c1->Write();
    for (Int_t i=0;i<nri;i++){
        hists[i]->Write();
    }
    ofile->Close();
}



