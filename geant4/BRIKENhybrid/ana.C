#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TCanvas.h>
#include <iostream>
#include <TGraph.h>
using namespace std;

void ana(Int_t nthreads=20){
    std::ifstream inpf("outGrid.txt");
    std::ofstream ofs("outGrid_joint.txt");
    while (inpf.good()){
        Double_t x,y,z,E;
        Int_t ndetect=0,nall=0;
        for (int i=0;i<nthreads;i++){
            Int_t ndetecti,nalli;
            inpf>>x>>y>>z>>E>>ndetecti>>nalli;
            cout<<x<<"\t"<<y<<"\t"<<z<<"\t"<<E<<"\t"<<ndetecti<<"\t"<<nalli<<endl;
            ndetect+=ndetecti;
            nall+=nalli;
        }
        cout<<nall<<"\t"<<ndetect<<endl;
        ofs<<x<<"\t"<<y<<"\t"<<z<<"\t"<<E<<"\t"<<(Double_t)ndetect/(Double_t)nall<<"\t"<<sqrt((Double_t)ndetect)/(Double_t)nall<<endl;
        cout<<"------------------------"<<endl;
    }
}
