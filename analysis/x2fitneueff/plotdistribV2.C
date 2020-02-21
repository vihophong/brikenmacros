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
#include <stdio.h>
#include <string>


const int nmaxfiles=100;

void makeg4histos(char* dir,Int_t nfiles)
{
    TFile *fileso[nmaxfiles];
    TTree *treeo[nmaxfiles];
    TH1F *hists[nmaxfiles];
    TH1F *hists2[nmaxfiles];


    //loop all files
    char tempname[500];
    std::ofstream ofs("geant4hitsall.txt");
    for (Int_t i=0;i<nfiles;i++){
        //! get ratio parameter

        sprintf(tempname,"%stempHist%d.root",dir,i);
        fileso[i]=TFile::Open(tempname);
        treeo[i]=(TTree*)fileso[i]->Get("BRIKEN");

        //! get original neutron energy
        treeo[i]->Draw(Form("partEnergy>>hpe%d(10000,0,10)",i),"","goff");
        sprintf(tempname,"hpe%i",i);
        hists[i]=(TH1F*)gDirectory->Get(tempname);
        Double_t partE=hists[i]->GetMean()*1e6;
        cout<<partE<<endl;

        //! get hit distribution histogram
        treeo[i]->Draw(Form("sqrt(x2*x2+y2*y2)>>hdist%d(200,0,50)",i),"id<140","goff");
        sprintf(tempname,"hdist%i",i);
        hists2[i]=(TH1F*)gDirectory->Get(tempname);
        sprintf(tempname,"hdist%i_%i",i,(Int_t)round(partE));
        hists2[i]->SetName(tempname);
        ofs<<tempname<<endl;
    }
    TFile* f0=new TFile("geant4hitsall.root","recreate");
    for (Int_t i=0;i<nfiles;i++){
        hists2[i]->Write();
    }
    f0->Close();
}
void plotdistr(char* inputfilename)
{
    TFile* fg4in=TFile::Open(inputfilename);
    TH1F* h1=(TH1F*)fg4in->Get("hdistr");
    TH1F* h4=(TH1F*)fg4in->Get("hdistrbwib");
    TH1F* h5=(TH1F*)fg4in->Get("hdistrbwbn");
    TH1F* h6=(TH1F*)fg4in->Get("hdistrbwbnbwib");
    TH1F* h1c=(TH1F*)h1->Clone();
    h1c->Add(h4,-1);
    h1c->Add(h5,-1);
    h1c->Add(h6,-1);
    //fg4in->Close();

    std::ifstream ifs("geant4hitsall.txt");
    std::string filelist[1000];

    Int_t nfiles=0;
    while (!ifs.eof()){
        ifs>>filelist[nfiles];
        cout<<filelist[nfiles]<<endl;
        nfiles++;
    }
    nfiles=nfiles-1;
    cout<<"There are "<<nfiles<<" histograms in total!"<<endl;

    TFile* fin=TFile::Open("geant4hitsall.root");

    TH1F *histsG4[100];
    Double_t partE[100];//in MeV
    Double_t KolmVal[100];



    for (Int_t i=0;i<nfiles;i++){
        int partEi,ii;
        sscanf((char*)filelist[i].c_str(),"hdist%d_%d",&ii,&partEi);
        partE[i]=(Double_t)partEi/1000.;
        //cout<<partE[i]<<endl;
        histsG4[i]=(TH1F*) fin->Get((char*)filelist[i].c_str());
        KolmVal[i]=histsG4[i]->KolmogorovTest(h1c,"UO");
    }


    //! find maximum
    Double_t maxKolval=0;
    Double_t maxKolvalE=0;
    Int_t maxKolvali=0;

    for (Int_t i=0;i<nfiles;i++){
        if (i==0) {
            maxKolval=KolmVal[i];
        }else{
            if (KolmVal[i]>maxKolval){
                maxKolval=KolmVal[i];
                maxKolvali=i;
                maxKolvalE=partE[i];
            }
        }

    }

    TGraph* gr=new TGraph(nfiles,partE,KolmVal);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(2);
    gr->SetLineColor(2);
    gr->SetMarkerSize(1.2);
    gr->SetLineWidth(2);


    TCanvas* c1=new TCanvas("c1","c1",900,700);
    c1->Draw();
    c1->Divide(2,2);
    c1->cd(1);
    gr->Draw("APL");

    //! sort all
    TH1F *histsall=(TH1F*)histsG4[0]->Clone();
    Int_t binwcontent[200];
    Int_t nbinwcontent=0;
    for (Int_t i=1;i<nfiles;i++) histsall->Add(histsG4[i]);
    for (Int_t i=1;i<201;i++){
        if (histsall->GetBinContent(i)>0){
            binwcontent[nbinwcontent]=i;
            nbinwcontent++;
        }
    }
    //c1->cd(3);
    //histsall->Draw();
    cout<<nbinwcontent<<endl;

    TH1F* h1csorted=new TH1F("h1csorted","h1csorted",nbinwcontent,0,nbinwcontent);
    for (Int_t j=1;j<201;j++){
        for (Int_t k=0;k<nbinwcontent;k++){
            if (j==binwcontent[k]){
                h1csorted->SetBinContent(k,h1c->GetBinContent(binwcontent[k]));
            }
        }
    }

    TH1F *histsG4sort[100];
    for (Int_t i=0;i<nfiles;i++) {
        histsG4sort[i]=new TH1F(Form("hneudist_sort%d",i),Form("hneudist_sort%d",i),nbinwcontent,0,nbinwcontent);
        for (Int_t j=1;j<201;j++){
            for (Int_t k=0;k<nbinwcontent;k++){
                if (j==binwcontent[k]){
                    histsG4sort[i]->SetBinContent(k,histsG4[i]->GetBinContent(binwcontent[k]));
                }
            }
        }
    }

    Double_t KolmValSort[100];

    for (Int_t i=0;i<nfiles;i++) {
        KolmValSort[i]=histsG4sort[i]->KolmogorovTest(h1csorted,"UO");
    }

    //! find maximum
    maxKolval=0;
    maxKolvalE=0;
    maxKolvali=0;

    for (Int_t i=0;i<nfiles;i++){
        if (i==0) {
            maxKolval=KolmValSort[i];
        }else{
            if (KolmValSort[i]>maxKolval){
                maxKolval=KolmValSort[i];
                maxKolvali=i;
                maxKolvalE=partE[i];
            }
        }
    }

    c1->cd(3);
    histsG4sort[maxKolvali]->Draw("hist");
    Double_t maxbinexp=h1csorted->GetBinContent(h1csorted->GetMaximumBin());
    Double_t maxbing4=histsG4sort[maxKolvali]->GetBinContent(histsG4sort[maxKolvali]->GetMaximumBin());
    h1csorted->Scale(maxbing4/maxbinexp);
    h1csorted->SetLineColor(2);
    h1csorted->Draw("hist same");

    c1->cd(4);
    TGraph* gr2=new TGraph(nfiles,partE,KolmValSort);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(2);
    gr2->SetLineColor(2);
    gr2->SetMarkerSize(1.2);
    gr2->SetLineWidth(2);
    gr2->Draw("APL");

    c1->cd(2);
    histsG4[maxKolvali]->Draw("hist");

    maxbinexp=h1c->GetBinContent(h1c->GetMaximumBin());
    maxbing4=histsG4[maxKolvali]->GetBinContent(histsG4[maxKolvali]->GetMaximumBin());
    h1c->Scale(maxbing4/maxbinexp);
    h1c->SetLineColor(2);
    h1c->Draw("hist same");


    //c1.SaveAs("test.root");
}
