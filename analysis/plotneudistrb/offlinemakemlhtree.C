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
    TGraph *grhists[nri];

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

        grhists[i]=new TGraph(gnpoints,gxx,gyy);
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

