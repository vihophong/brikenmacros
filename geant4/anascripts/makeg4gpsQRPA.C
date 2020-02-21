#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TCanvas.h>
#include <iostream>
#include <TGraph.h>

#define kmaxlines 1000
#define knormalization_de 0.1
void makeg4gpsQRPA(char* inputfile, char* outputfile,Int_t rebin=1)
{
    std::ifstream inpf(inputfile);
    Double_t Eedges[kmaxlines];

    Double_t Emin[kmaxlines];
    Double_t Emax[kmaxlines];
    Double_t neuspec[kmaxlines];

    Double_t temp;


    Int_t nlines=0,ientry=0;
    std::string line;


    while (std::getline(inpf, line))
    {
        std::istringstream iss(line);
        if (nlines>2){
            iss>>Emin[ientry]>>Emax[ientry]>>temp>>neuspec[ientry]>>temp>>temp;
            Eedges[ientry]=Emin[ientry];
            std::cout<<Emin[ientry]<<"\t"<<Emax[ientry]<<"\t"<<neuspec[ientry]<<std::endl;
            ientry++;
        }
        nlines++;
    }
    Eedges[ientry]=Emax[ientry-1];


    TH1F* hSpec=new TH1F("hSpec","hSpec",ientry,Eedges);

    for (Int_t i=0;i<ientry;i++){
        hSpec->SetBinContent(i+1,neuspec[i]);
    }
    hSpec->SetLineColor(2);


    for (Int_t i=0;i<ientry;i++){
        std::cout<<hSpec->GetBinCenter(i+1)<<"\t"<<hSpec->GetBinContent(i+1)<<"\t"<<hSpec->GetBinWidth(i+1)<<std::endl;
    }
    TH1F* hSpecRebin=(TH1F*)hSpec->Clone();
    hSpecRebin->SetName("hSpecRebin");
    hSpecRebin->Rebin(rebin);
    hSpecRebin->SetLineColor(3);
    hSpecRebin->Draw("hist");
    hSpec->Draw("same");

    std::ofstream oupf(outputfile);
    nlines=0;
    oupf<<"/gps/particle neutron\n/gps/pos/type Point\n/gps/pos/centre 0 0 5 mm\n/gps/ang/type iso\n/gps/ene/type Arb\n/gps/hist/type arb\n"<<std::endl;
    for (Int_t i=0;i<hSpecRebin->GetNbinsX();i++){
        if (hSpecRebin->GetBinContent(i+1)>0){
            std::cout<<hSpecRebin->GetBinCenter(i+1)<<"\t"<<hSpecRebin->GetBinContent(i+1)<<"\t"<<hSpecRebin->GetBinWidth(i+1)<<std::endl;
            oupf<<"/gps/hist/point\t"<<hSpecRebin->GetBinCenter(i+1)<<"\t"<<hSpecRebin->GetBinContent(i+1)/hSpecRebin->GetBinWidth(i+1)<<std::endl;
            nlines++;
        }
    }
    oupf<<"\n/gps/hist/inter Spline"<<std::endl;
    //oupf<<"\n/gps/hist/inter Lin"<<std::endl;
    std::cout<<nlines<<" lines"<<std::endl;

    char tempchar[1000];
    sprintf(tempchar,"%s.root",outputfile);
    hSpec->SaveAs(tempchar);
}
