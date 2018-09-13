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
void offlinemakemlhtree(char* input)
{
    TString ri[100];
    Int_t nri=21;
    ri[0]="Sn134";
    ri[1]="Sn135";
    ri[2]="Sn136";
    ri[3]="Sn137";
    ri[4]="Sn138";
    ri[5]="Sn139";
    ri[6]="In131";
    ri[7]="In132";
    ri[8]="In133";
    ri[9]="In134";
    ri[10]="In135";
    ri[11]="In136";
    ri[12]="Cd129";
    ri[13]="Cd130";
    ri[14]="Cd131";
    ri[15]="Cd132";
    ri[16]="Cd133";
    ri[17]="Cd134";
    ri[18]="Ag130";
    ri[19]="Ag131";
    ri[20]="Ag132";

    gROOT->ProcessLine(".L makemlhtree.C");
    for (Int_t i=0;i<nri;i++){
        char tempchar1[1000];
        sprintf(tempchar1,"outree%s.root",(char*)ri[i].Data());
        cout<<Form("makemlhtree o%d(\"%s\",\"%s\");",i,input,(char*)ri[i].Data())<<endl;
        gROOT->ProcessLine(Form("makemlhtree o%d(\"%s\",\"%s\");",i,input,(char*)ri[i].Data()));
        gROOT->ProcessLine(Form("o%d.Loop(\"fittrees/%s\");",i,tempchar1));
    }

}
void offlineplotneutrondist(char* input)
{
    TString ri[100];
    Int_t nri=21;
    ri[0]="Sn134";
    ri[1]="Sn135";
    ri[2]="Sn136";
    ri[3]="Sn137";
    ri[4]="Sn138";
    ri[5]="Sn139";
    ri[6]="In131";
    ri[7]="In132";
    ri[8]="In133";
    ri[9]="In134";
    ri[10]="In135";
    ri[11]="In136";
    ri[12]="Cd129";
    ri[13]="Cd130";
    ri[14]="Cd131";
    ri[15]="Cd132";
    ri[16]="Cd133";
    ri[17]="Cd134";
    ri[18]="Ag130";
    ri[19]="Ag131";
    ri[20]="Ag132";

    gROOT->ProcessLine(".L makemlhtree.C");
    for (Int_t i=0;i<nri;i++){
        char tempchar1[1000];
        sprintf(tempchar1,"neuhit%s.root",(char*)ri[i].Data());
        cout<<Form("makemlhtree o%d(\"%s\",\"%s\");",i,input,(char*)ri[i].Data())<<endl;
        gROOT->ProcessLine(Form("makemlhtree o%d(\"%s\",\"%s\");",i,input,(char*)ri[i].Data()));
        gROOT->ProcessLine(Form("o%d.PlotNeuHitPattern(\"neutron_distr/%s\");",i,tempchar1));
    }

}

void plotneudists(char* outfilename,Int_t isnorml=0)
{
    Int_t colorcode[]={1,2,3,4,5,6,7,9,11,12,13,14,15,16,17,18,19,20,21,22,23};
    TString ri[100];
    Bool_t flag[100];
    Int_t nri=21;
    //ri[0]="Sn134";
    //ri[1]="Sn135";
    //ri[2]="Sn136";
    ri[0]="1000keV";
    ri[1]="2000keV";
    ri[2]="3000keV";

    ri[3]="Sn137";
    ri[4]="Sn138";
    ri[5]="Sn139";
    ri[6]="In131";
    ri[7]="In132";
    ri[8]="In133";
    ri[9]="In134";
    ri[10]="In135";
    ri[11]="In136";
    ri[12]="Cd129";
    ri[13]="Cd130";
    ri[14]="Cd131";
    ri[15]="Cd132";
    ri[16]="Cd133";
    ri[17]="Cd134";
    ri[18]="Ag130";
    ri[19]="Ag131";
    ri[20]="Ga83";
    //ri[20]="Ag132";

    flag[0]=true;
    flag[1]=true;
    flag[2]=true;
    flag[3]=false;
    flag[4]=false;
    flag[5]=false;
    flag[6]=false;
    flag[7]=false;
    flag[8]=false;
    flag[9]=false;
    flag[10]=false;
    flag[11]=false;
    flag[12]=false;
    flag[13]=false;
    flag[14]=false;
    flag[15]=true;
    flag[16]=false;
    flag[17]=false;
    flag[18]=false;
    flag[19]=false;
    flag[20]=true;
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
        sprintf(tempchar1,"neutron_distr/neuhit%s.root",(char*)ri[i].Data());
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


void testtscorr()
{
    unsigned long long tsarr[]={6616109261992	,
    6616109262192	,
    6616109262392	,
    6616109262592	,
    6616109262792	,
    6616109262992	,
    6616109263192	,
    6616109263392	,
    6616109263592	,
    6616109263792	,
    6616109263992	,
    6616109264192	,
    6616109264392	,
    6616109264592	,
    6616109264792	,
    6616109264992
                               };
    unsigned long long fprev_ASICS_ts=0;
    unsigned long long ffirst_ASICS_ts=0;
    unsigned short fprev_ASICS_cnt=0;

    for (Int_t i=0;i<16;i++){
        if ((tsarr[i]-fprev_ASICS_ts)==200){
            if(ffirst_ASICS_ts==0){
                ffirst_ASICS_ts=fprev_ASICS_ts;
            }
            cout<<tsarr[i]<<"-"<<ffirst_ASICS_ts<<endl;
            fprev_ASICS_cnt++;
        }else{
            fprev_ASICS_cnt=0;
        }
        fprev_ASICS_ts=tsarr[i];
    }
}
