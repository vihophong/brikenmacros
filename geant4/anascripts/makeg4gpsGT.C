#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TCanvas.h>
#include <iostream>
#include <TGraph.h>

#define kMAXE 10

#define kmaxlines 10000
#define knormalization_de 0.1
void makeg4gpsGT(char* inputfile, char* outputfile,Int_t rebin=8)
{

    std::ifstream inpf(inputfile);
    Double_t Eedges[kmaxlines];

    Double_t neuspec[kmaxlines];


    Int_t nlines=0,ientry=0;
    std::string line;


    while (std::getline(inpf, line))
    {
        std::istringstream iss(line);
        //if (nlines>2){
            iss>>Eedges[ientry]>>neuspec[ientry];
            std::cout<<Eedges[ientry]<<"\t"<<neuspec[ientry]<<std::endl;
            ientry++;
        //}
        nlines++;
    }
    std::cout<<ientry<<std::endl;
    TH1F* hSpec=new TH1F("hSpec","hSpec",ientry-1,Eedges);


    for (Int_t i=0;i<ientry-1;i++){
        hSpec->SetBinContent(i+1,neuspec[i]);
    }
    hSpec->SetLineColor(2);


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
            if (hSpecRebin->GetBinCenter(i+1)<kMAXE){
                std::cout<<hSpecRebin->GetBinCenter(i+1)<<"\t"<<hSpecRebin->GetBinContent(i+1)<<"\t"<<hSpecRebin->GetBinWidth(i+1)<<std::endl;
                oupf<<"/gps/hist/point\t"<<hSpecRebin->GetBinCenter(i+1)<<"\t"<<hSpecRebin->GetBinContent(i+1)/hSpecRebin->GetBinWidth(i+1)<<std::endl;
                nlines++;
            }

        }
    }
    oupf<<"\n/gps/hist/inter Spline"<<std::endl;
    //oupf<<"\n/gps/hist/inter Lin"<<std::endl;
    std::cout<<nlines<<" lines"<<std::endl;

    char tempchar[1000];
    sprintf(tempchar,"%s.root",outputfile);
    hSpec->SaveAs(tempchar);


}
