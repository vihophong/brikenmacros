#include "stdio.h"
#include "string.h"
#include "TROOT.h"
#include "TString.h"
#include "TH1.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"

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

    for (Int_t i=0;i<10;i++){
        char tempname[500];
        sprintf(tempname,"tempHist%d.root",i);
        fileso[i]=TFile::Open(tempname);
        treeo[i]=(TTree*)fileso[i]->Get("BRIKEN");
        treeo[i]->Draw(Form("sqrt(x*x+y*y)>>hcont%s(200,0,2000)",(char*)entryname[i]->Data()),"","goff");

        sprintf(tempname,"hcont%s",(char*)entryname[i]->Data());
        hists[i]=(TH1F*)gDirectory->Get(tempname);
        treeo[i]->Draw(Form("sqrt(x2*x2+y2*y2)>>hdisc%s(200,0,2000)",(char*)entryname[i]->Data()),"","goff");
        sprintf(tempname,"hdisc%s",(char*)entryname[i]->Data());
        hists2[i]=(TH1F*)gDirectory->Get(tempname);

        hists[i]->SetLineColor(colorcode[i]);
        hists[i]->SetLineWidth(2);
        hists2[i]->SetLineColor(colorcode[i]);
        hists2[i]->SetLineWidth(2);

        leg->AddEntry(hists[i],(char*)legentryname[i]->Data());
        leg2->AddEntry(hists2[i],(char*)legentryname[i]->Data());
    }
    for (Int_t i=0;i<10;i++){
        if (i==0) hists[i]->Draw();
        else hists[i]->Draw("same");
    }
    leg->Draw();

    c1->cd(2);
    for (Int_t i=0;i<10;i++){
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
    for (Int_t i=0;i<10;i++){
        char tempname[500];
        sprintf(tempname,"hist/%shit.root",(char*)entryname[i]->Data());
        ofile[i]=new TFile(tempname,"recreate");
        hists2[i]->SetName("h1");
        hists2[i]->Write();
        ofile[i]->Close();
    }
}

void plotneudists(char* outfilename,Int_t isnorml=0)
{
    Int_t colorcode[]={1,2,3,4,5,6,7,9,11,12,13,14,15,16,17,18,19,20,21,22,23};
    TString ri[100];
    Bool_t flag[100];
    Int_t nri=21;
    ri[0]="Sn134";
    ri[1]="Sn135";
    ri[2]="Sn136";
    //ri[0]="1000keV";
    //ri[1]="2000keV";
    //ri[2]="3000keV";

    ri[3]="Sn137";
    ri[4]="Sn138";
    ri[5]="Sn138";
    ri[6]="In131";
    ri[7]="In132";
    ri[8]="In133";
    ri[9]="In134";
    ri[10]="In135";
    ri[11]="In135";
    ri[12]="Cd130";
    ri[13]="Cd130";
    ri[14]="Cd131";
    ri[15]="Cd132";
    ri[16]="Cd133";
    ri[17]="Cd133";
    ri[18]="Ag130";
    ri[19]="Ag131";
    //ri[20]="Ag131";
    ri[20]="Ag131";

    flag[0]=false;
    flag[1]=false;
    flag[2]=true;
    flag[3]=true;
    flag[4]=true;
    flag[5]=false;
    flag[6]=false;
    flag[7]=false;
    flag[8]=true;
    flag[9]=true;
    flag[10]=true;
    flag[11]=false;
    flag[12]=false;
    flag[13]=false;
    flag[14]=true;
    flag[15]=true;
    flag[16]=false;
    flag[17]=false;
    flag[18]=false;
    flag[19]=false;
    flag[20]=false;
    /*
    for (Int_t i=0;i<21;i++) {
        if (i<4)
        flag[i]=true;
        else flag[i]=false;
    }
    */

    TFile *fileso[nri];
    TH1F *hists[nri];
    TGraphErrors *grhists[nri];

    TCanvas* c1=new TCanvas("c1","c1",900,700);
    TLegend* leg = new TLegend(0.2, 0.2, .8, .8);


    Double_t bin8countsnorm=100.;
    Double_t bin8counts;


    char temppchar[500];
    sprintf(temppchar,"%s.txt",outfilename);
    std::ofstream ofs(temppchar);

    for (Int_t i=0;i<nri;i++){
        Double_t gxx[200];
        Double_t gyy[200];
        Double_t gyyerr[200];
        Int_t gnpoints=0;

        char tempchar1[1000];
        sprintf(tempchar1,"neuhit/%shit.root",(char*)ri[i].Data());
        fileso[i]=TFile::Open(tempchar1);
        hists[i]=(TH1F*)gDirectory->Get("h1");
        hists[i]->SetName((char*)ri[i].Data());
        hists[i]->SetMarkerStyle(kFullCircle);
        //hists[i]->Scale(norm, "width");
        bin8counts=(Double_t)hists[i]->GetBinContent(8);
        if (isnorml) hists[i]->Scale(bin8countsnorm/bin8counts);
        //hists[i]->SetLineWidth(0);

        Double_t seppts=13.;
        Int_t innercnt=0;
        Int_t outercnt=0;
        for (Int_t j=0;j<200;j++){
            if (hists[i]->GetBinContent(j+1)>0){
                gxx[gnpoints]=sqrt(hists[i]->GetXaxis()->GetBinCenter(j+1));
                gyy[gnpoints]=hists[i]->GetBinContent(j+1);
                gyyerr[gnpoints]=hists[i]->GetBinError(j+1);
                if (gxx[gnpoints]<seppts) innercnt+=gyy[gnpoints];
                else outercnt+=gyy[gnpoints];
                gnpoints++;
            }
        }
        /*
        ofs<<ri[i]<<"\t";
        for (Int_t k=0;k<gnpoints;k++) ofs<<gxx[k]<<"\t";
        ofs<<endl;
        */
        ofs<<ri[i]<<"\t";
        //for (Int_t k=0;k<gnpoints;k++) ofs<<gyy[k]<<"\t";
        ofs<<innercnt<<"\t"<<outercnt;
        ofs<<endl;

        grhists[i]=new TGraphErrors(gnpoints,gxx,gyy,0,gyyerr);
        grhists[i]->SetMarkerStyle(kFullCircle);
        grhists[i]->SetFillColor(0);
    }
    Int_t ncolor=0;
    for (Int_t i=0;i<nri;i++){
        if (flag[i]){
            hists[i]->SetMarkerColor(colorcode[ncolor]);
            hists[i]->SetLineColor(colorcode[ncolor]);
            grhists[i]->SetMarkerColor(colorcode[ncolor]);
            grhists[i]->SetLineColor(colorcode[ncolor]);
            leg->AddEntry(grhists[i],(char*)ri[i].Data());
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

