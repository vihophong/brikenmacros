#include "TGraph.h"
#include "TSpline.h"
#include "TTree.h"
#include "TCanvas.h"
void calattenuation()
{
    std::ifstream ifs("xraytable_Al.txt");
    Double_t xdots[1000];
    Double_t ydots[1000];


    Int_t nlines=0;

    while (!ifs.eof()){
        ifs>>xdots[nlines]>>ydots[nlines];
        cout<<xdots[nlines]<<"\t"<<ydots[nlines]<<endl;
        nlines++;
    }
    nlines--;

    TCanvas *c1=new TCanvas("c1","c1",900,700);
    TGraph* gr=new TGraph(nlines,xdots,ydots);
    gr->SetMarkerStyle(20);
    gr->Draw("AP");
    gPad->SetLogy();
    gPad->SetLogx();


    TSpline3* sp=new TSpline3("grs",gr);
    sp->SetLineColor(2);
    sp->Draw("same");

    TCanvas *c2=new TCanvas("c2","c2",900,700);
    cout<<sp->Eval(0.0567)<<endl;

    Double_t thickness[]={1,2,3,4,5};
    Double_t density=2.33;
    Double_t ene[2000];


    TGraph* grabp[5];
    for (Int_t i=0;i<5;i++){
        Double_t arr[2000];
        for (Int_t j=1;j<1001;j++){
            ene[j-1]=(Double_t)j;
            Double_t eval=sp->Eval(ene[j-1]/1000.);
            arr[j-1]=exp(-eval*density*thickness[i]/10)*100;
        }
        grabp[i]=new TGraph(1000,ene,arr);
    }


    TLegend* leg=new TLegend(0.1,0.7,0.5,0.9);
    grabp[0]->Draw("AL");
    grabp[0]->SetFillColor(0);
    grabp[0]->SetLineWidth(2);
    leg->AddEntry(grabp[0],"1 mm Al");
    for (Int_t i=1;i<5;i++){
        grabp[i]->SetFillColor(0);
        grabp[i]->SetLineWidth(2);
        grabp[i]->SetLineColor(i+1);
        grabp[i]->Draw("L SAME");
        leg->AddEntry(grabp[i],Form("%d mm Al",i+1));
    }
    leg->Draw();


}
