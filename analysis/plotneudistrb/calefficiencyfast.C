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

void calefficiencyfast(char* ri,Double_t inputr,Double_t inputrminus,Double_t inputrplus)
{
    TFile* file1=TFile::Open("outfilereffvse.root");
    TGraph* gr=(TGraph*)file1->Get("rvse");
    TGraph* gr2=(TGraph*)file1->Get("evsr");
    TGraph* greff=(TGraph*)file1->Get("effvsr");


    TCanvas* c1=new TCanvas("c1","c1",900,700);
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



    TFile* fileout=new TFile(Form("effcal%s.root",ri),"recreate");
    fileout->cd();
    c1->Write();

    fileout->Close();

    //! summary
    //!
    std::ofstream str("outimpandratio.txt",std::ios::app);
    str<<ri<<"\t"<<inputr<<"\t"<<gr2->Eval(inputr)<<"\t"<<greff->Eval(gr2->Eval(inputr))*100<<"\t"
      <<inputrplus<<"\t"<<gr2->Eval(inputrplus)<<"\t"<<greff->Eval(gr2->Eval(inputrplus))*100
      <<"\t"<<inputrminus<<"\t"<<gr2->Eval(inputrminus)<<"\t"<<greff->Eval(gr2->Eval(inputrminus))*100<<endl;
    cout<<"input_ratio\tmean_energy(MeV)\tefficiency(%)"<<endl;
    cout<<inputr<<"\t"<<gr2->Eval(inputr)<<"\t"<<greff->Eval(gr2->Eval(inputr))*100<<endl;
    cout<<inputrplus<<"\t"<<gr2->Eval(inputrplus)<<"\t"<<greff->Eval(gr2->Eval(inputrplus))*100<<endl;
    cout<<inputrminus<<"\t"<<gr2->Eval(inputrminus)<<"\t"<<greff->Eval(gr2->Eval(inputrminus))*100<<endl;

}
