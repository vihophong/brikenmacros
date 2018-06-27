#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TGraphErrors.h>

#include <TMultiGraph.h>

#include <TLatex.h>

#include <string.h>

void plotresults()
{
    string s[]={"In132",
                "In133",
                "In134",
                "In135",
                "Sn134",
                "Sn135",
                "Sn136",
                "Sn137",
                "Sn138",
                "Cd131",
                "Cd132",
                "Cd133",
                "Ag130",
                "Ag131"};
    cout<<s[0]<<endl;
    TGraphErrors* gr1hl=new TGraphErrors("grtxt/fitresults1.txt","%*s %*lg %lg %lg %lg %*lg %*lg %*lg %*lg");
    gr1hl->SetMarkerStyle(21);
    gr1hl->SetMarkerColor(1);

    for (Int_t i=0;i<gr1hl->GetN();i++) {
        TLatex *latex = new TLatex(gr1hl->GetX()[i], gr1hl->GetY()[i],s[i].c_str());
        latex->SetTextSize(0.04);
        gr1hl->GetListOfFunctions()->Add(latex);;
    }

    TGraphErrors* gr2hl=new TGraphErrors("grtxt/fitresults2.txt","%*s %*lg %lg %lg %lg %*lg %*lg %*lg %*lg");
    gr2hl->SetMarkerStyle(22);
    gr2hl->SetMarkerColor(2);
    TGraphErrors* gr3hl=new TGraphErrors("grtxt/fitresults3.txt","%*s %*lg %lg %lg %lg %*lg %*lg %*lg %*lg");
    gr3hl->SetMarkerStyle(23);
    gr3hl->SetMarkerColor(3);
    TGraphErrors* gr4hl=new TGraphErrors("grtxt/fitresults4.txt","%*s %*lg %lg %lg %lg %*lg %*lg %*lg %*lg");
    gr4hl->SetMarkerStyle(24);
    gr4hl->SetMarkerColor(4);
    TGraphErrors* gr5hl=new TGraphErrors("grtxt/fitresults5.txt","%*s %*lg %lg %lg %lg %*lg %*lg %*lg %*lg");
    gr5hl->SetMarkerStyle(25);
    gr5hl->SetMarkerColor(5);
    TGraphErrors* gr6hl=new TGraphErrors("grtxt/fitresults6.txt","%*s %*lg %lg %lg %lg %*lg %*lg %*lg %*lg");
    gr6hl->SetMarkerStyle(26);
    gr6hl->SetMarkerColor(6);
    TGraphErrors* gr7hl=new TGraphErrors("grtxt/fitresults7.txt","%*s %*lg %lg %lg %lg %*lg %*lg %*lg %*lg");
    gr7hl->SetMarkerStyle(27);
    gr7hl->SetMarkerColor(7);
    TGraphErrors* gr8hl=new TGraphErrors("grtxt/fitresults8.txt","%*s %*lg %lg %lg %lg %*lg %*lg %*lg %*lg");
    gr8hl->SetMarkerStyle(28);
    gr8hl->SetMarkerColor(8);

    TCanvas* cc1=new TCanvas("cc1","cc1",900,700);
    cc1->Divide(2,2);
    cc1->cd(1);
    TMultiGraph* mghl=new TMultiGraph();
    gr1hl->SetTitle("100 keV");
    gr2hl->SetTitle("120 keV");
    gr3hl->SetTitle("140 keV");
    gr4hl->SetTitle("160 keV");
    gr5hl->SetTitle("180 keV");
    gr6hl->SetTitle("200 keV");
    gr7hl->SetTitle("220 keV");
    gr8hl->SetTitle("240 keV");

    mghl->Add(gr1hl);
    mghl->Add(gr2hl);
    mghl->Add(gr3hl);
    mghl->Add(gr4hl);
    mghl->Add(gr5hl);
    mghl->Add(gr6hl);
    mghl->Add(gr7hl);
    mghl->Add(gr8hl);
    mghl->Draw("ap");
    cc1->cd(1)->BuildLegend();


    TGraphErrors* gr1p1n=new TGraphErrors("grtxt/fitresults1.txt","%*s %*lg %lg %*lg %*lg %lg %lg %*lg %*lg");
        gr1p1n->SetMarkerStyle(21);
        gr1p1n->SetMarkerColor(1);

     for (Int_t i=0;i<gr1p1n->GetN();i++) {
         TLatex *latex = new TLatex(gr1p1n->GetX()[i], gr1p1n->GetY()[i],s[i].c_str());
         latex->SetTextSize(0.04);
         gr1p1n->GetListOfFunctions()->Add(latex);;
     }

        TGraphErrors* gr2p1n=new TGraphErrors("grtxt/fitresults2.txt","%*s %*lg %lg %*lg %*lg %lg %lg %*lg %*lg");
        gr2p1n->SetMarkerStyle(22);
        gr2p1n->SetMarkerColor(2);
        TGraphErrors* gr3p1n=new TGraphErrors("grtxt/fitresults3.txt","%*s %*lg %lg %*lg %*lg %lg %lg %*lg %*lg");
        gr3p1n->SetMarkerStyle(23);
        gr3p1n->SetMarkerColor(3);
        TGraphErrors* gr4p1n=new TGraphErrors("grtxt/fitresults4.txt","%*s %*lg %lg %*lg %*lg %lg %lg %*lg %*lg");
        gr4p1n->SetMarkerStyle(24);
        gr4p1n->SetMarkerColor(4);
        TGraphErrors* gr5p1n=new TGraphErrors("grtxt/fitresults5.txt","%*s %*lg %lg %*lg %*lg %lg %lg %*lg %*lg");
        gr5p1n->SetMarkerStyle(25);
        gr5p1n->SetMarkerColor(5);
        TGraphErrors* gr6p1n=new TGraphErrors("grtxt/fitresults6.txt","%*s %*lg %lg %*lg %*lg %lg %lg %*lg %*lg");
        gr6p1n->SetMarkerStyle(26);
        gr6p1n->SetMarkerColor(6);
        TGraphErrors* gr7p1n=new TGraphErrors("grtxt/fitresults7.txt","%*s %*lg %lg %*lg %*lg %lg %lg %*lg %*lg");
        gr7p1n->SetMarkerStyle(27);
        gr7p1n->SetMarkerColor(7);
        TGraphErrors* gr8p1n=new TGraphErrors("grtxt/fitresults8.txt","%*s %*lg %lg %*lg %*lg %lg %lg %*lg %*lg");
        gr8p1n->SetMarkerStyle(28);
        gr8p1n->SetMarkerColor(8);



        cc1->cd(2);
        TMultiGraph* mgp1n=new TMultiGraph();
        gr1p1n->SetTitle("100 keV");
        gr2p1n->SetTitle("120 keV");
        gr3p1n->SetTitle("140 keV");
        gr4p1n->SetTitle("160 keV");
        gr5p1n->SetTitle("180 keV");
        gr6p1n->SetTitle("200 keV");
        gr7p1n->SetTitle("220 keV");
        gr8p1n->SetTitle("240 keV");

        mgp1n->Add(gr1p1n);
        mgp1n->Add(gr2p1n);
        mgp1n->Add(gr3p1n);
        mgp1n->Add(gr4p1n);
        mgp1n->Add(gr5p1n);
        mgp1n->Add(gr6p1n);
        mgp1n->Add(gr7p1n);
        mgp1n->Add(gr8p1n);
        mgp1n->Draw("ap");
        cc1->cd(2)->BuildLegend();

        TGraphErrors* gr1p2n=new TGraphErrors("grtxt/fitresults1.txt","%*s %*lg %lg %*lg %*lg %*lg %*lg %lg %lg");
            gr1p2n->SetMarkerStyle(21);
            gr1p2n->SetMarkerColor(1);


            for (Int_t i=0;i<gr1p2n->GetN();i++) {
                TLatex *latex = new TLatex(gr1p2n->GetX()[i], gr1p2n->GetY()[i],s[i].c_str());
                latex->SetTextSize(0.04);
                gr1p2n->GetListOfFunctions()->Add(latex);;
            }

            TGraphErrors* gr2p2n=new TGraphErrors("grtxt/fitresults2.txt","%*s %*lg %lg %*lg %*lg %*lg %*lg %lg %lg");
            gr2p2n->SetMarkerStyle(22);
            gr2p2n->SetMarkerColor(2);
            TGraphErrors* gr3p2n=new TGraphErrors("grtxt/fitresults3.txt","%*s %*lg %lg %*lg %*lg %*lg %*lg %lg %lg");
            gr3p2n->SetMarkerStyle(23);
            gr3p2n->SetMarkerColor(3);
            TGraphErrors* gr4p2n=new TGraphErrors("grtxt/fitresults4.txt","%*s %*lg %lg %*lg %*lg %*lg %*lg %lg %lg");
            gr4p2n->SetMarkerStyle(24);
            gr4p2n->SetMarkerColor(4);
            TGraphErrors* gr5p2n=new TGraphErrors("grtxt/fitresults5.txt","%*s %*lg %lg %*lg %*lg %*lg %*lg %lg %lg");
            gr5p2n->SetMarkerStyle(25);
            gr5p2n->SetMarkerColor(5);
            TGraphErrors* gr6p2n=new TGraphErrors("grtxt/fitresults6.txt","%*s %*lg %lg %*lg %*lg %*lg %*lg %lg %lg");
            gr6p2n->SetMarkerStyle(26);
            gr6p2n->SetMarkerColor(6);
            TGraphErrors* gr7p2n=new TGraphErrors("grtxt/fitresults7.txt","%*s %*lg %lg %*lg %*lg %*lg %*lg %lg %lg");
            gr7p2n->SetMarkerStyle(27);
            gr7p2n->SetMarkerColor(7);
            TGraphErrors* gr8p2n=new TGraphErrors("grtxt/fitresults8.txt","%*s %*lg %lg %*lg %*lg %*lg %*lg %lg %lg");
            gr8p2n->SetMarkerStyle(28);
            gr8p2n->SetMarkerColor(8);

            cc1->cd(3);
            TMultiGraph* mgp2n=new TMultiGraph();
            gr1p2n->SetTitle("100 keV");
            gr2p2n->SetTitle("120 keV");
            gr3p2n->SetTitle("140 keV");
            gr4p2n->SetTitle("160 keV");
            gr5p2n->SetTitle("180 keV");
            gr6p2n->SetTitle("200 keV");
            gr7p2n->SetTitle("220 keV");
            gr8p2n->SetTitle("240 keV");

            mgp2n->Add(gr1p2n);
            mgp2n->Add(gr2p2n);
            mgp2n->Add(gr3p2n);
            mgp2n->Add(gr4p2n);
            mgp2n->Add(gr5p2n);
            mgp2n->Add(gr6p2n);
            mgp2n->Add(gr7p2n);
            mgp2n->Add(gr8p2n);
            mgp2n->Draw("ap");
            cc1->cd(3)->BuildLegend();

}
void plotresults2()
{
    TCanvas* cc=new TCanvas("cc","cc",900,700);

    TGraphErrors* grIn133hl=new TGraphErrors("grtxt/fitIn133.txt","%lg %lg %lg %*lg %*lg %*lg %*lg");
    grIn133hl->SetMarkerStyle(21);
    grIn133hl->SetMarkerColor(2);

    TGraphErrors* grIn133p1n=new TGraphErrors("grtxt/fitIn133.txt","%lg %*lg %*lg %lg %lg %*lg %*lg");
    grIn133p1n->SetMarkerStyle(21);
    grIn133p1n->SetMarkerColor(2);

    TGraphErrors* grIn133p2n=new TGraphErrors("grtxt/fitIn133.txt","%lg %*lg %*lg %*lg %*lg %lg %lg");
    grIn133p2n->SetMarkerStyle(21);
    grIn133p2n->SetMarkerColor(2);


    TGraphErrors* grSn136hl=new TGraphErrors("grtxt/fitSn136.txt","%lg %lg %lg %*lg %*lg %*lg %*lg");
    grSn136hl->SetMarkerStyle(21);
    grSn136hl->SetMarkerColor(2);

    TGraphErrors* grSn136p1n=new TGraphErrors("grtxt/fitSn136.txt","%lg %*lg %*lg %lg %lg %*lg %*lg");
    grSn136p1n->SetMarkerStyle(21);
    grSn136p1n->SetMarkerColor(2);

    TGraphErrors* grSn136p2n=new TGraphErrors("grtxt/fitSn136.txt","%lg %*lg %*lg %*lg %*lg %lg %lg");
    grSn136p2n->SetMarkerStyle(21);
    grSn136p2n->SetMarkerColor(2);

    TGraphErrors* grCd132hl=new TGraphErrors("grtxt/fitCd132.txt","%lg %lg %lg %*lg %*lg %*lg %*lg");
    grCd132hl->SetMarkerStyle(21);
    grCd132hl->SetMarkerColor(2);

    TGraphErrors* grCd132p1n=new TGraphErrors("grtxt/fitCd132.txt","%lg %*lg %*lg %lg %lg %*lg %*lg");
    grCd132p1n->SetMarkerStyle(21);
    grCd132p1n->SetMarkerColor(2);

    TGraphErrors* grCd132p2n=new TGraphErrors("grtxt/fitCd132.txt","%lg %*lg %*lg %*lg %*lg %lg %lg");
    grCd132p2n->SetMarkerStyle(21);
    grCd132p2n->SetMarkerColor(2);

    TGraphErrors* grAg130hl=new TGraphErrors("grtxt/fitAg130.txt","%lg %lg %lg %*lg %*lg %*lg %*lg");
    grAg130hl->SetMarkerStyle(21);
    grAg130hl->SetMarkerColor(2);

    TGraphErrors* grAg130p1n=new TGraphErrors("grtxt/fitAg130.txt","%lg %*lg %*lg %lg %lg %*lg %*lg");
    grAg130p1n->SetMarkerStyle(21);
    grAg130p1n->SetMarkerColor(2);

    TGraphErrors* grAg130p2n=new TGraphErrors("grtxt/fitAg130.txt","%lg %*lg %*lg %*lg %*lg %lg %lg");
    grAg130p2n->SetMarkerStyle(21);
    grAg130p2n->SetMarkerColor(2);



    grSn136hl->SetTitle("Half-life of Sn136");
    grSn136hl->GetXaxis()->SetTitle("Threshold (keV)");
    grSn136hl->GetYaxis()->SetTitle("Half-life (s)");

    grIn133hl->SetTitle("Half-life of In133");
    grIn133hl->GetXaxis()->SetTitle("Threshold (keV)");
    grIn133hl->GetYaxis()->SetTitle("Half-life (s)");

    grCd132hl->SetTitle("Half-life of Cd132");
    grCd132hl->GetXaxis()->SetTitle("Threshold (keV)");
    grCd132hl->GetYaxis()->SetTitle("Half-life (s)");

    grAg130hl->SetTitle("Half-life of Ag130");
    grAg130hl->GetXaxis()->SetTitle("Threshold (keV)");
    grAg130hl->GetYaxis()->SetTitle("Half-life (s)");


    grSn136p1n->SetTitle("P1n of Sn136");
        grSn136p1n->GetXaxis()->SetTitle("Threshold (keV)");
        grSn136p1n->GetYaxis()->SetTitle("P1n (/100 %)");

        grIn133p1n->SetTitle("P1n of In133");
        grIn133p1n->GetXaxis()->SetTitle("Threshold (keV)");
        grIn133p1n->GetYaxis()->SetTitle("P1n (/100 %)");

        grCd132p1n->SetTitle("P1n of Cd132");
        grCd132p1n->GetXaxis()->SetTitle("Threshold (keV)");
        grCd132p1n->GetYaxis()->SetTitle("P1n (/100 %)");

        grAg130p1n->SetTitle("P1n of Ag130");
        grAg130p1n->GetXaxis()->SetTitle("Threshold (keV)");
        grAg130p1n->GetYaxis()->SetTitle("P1n (/100 %)");

        grSn136p2n->SetTitle("P2n of Sn136");
            grSn136p2n->GetXaxis()->SetTitle("Threshold (keV)");
            grSn136p2n->GetYaxis()->SetTitle("P2n (/100 %)");

            grIn133p2n->SetTitle("P2n of In133");
            grIn133p2n->GetXaxis()->SetTitle("Threshold (keV)");
            grIn133p2n->GetYaxis()->SetTitle("P2n (/100 %)");

            grCd132p2n->SetTitle("P2n of Cd132");
            grCd132p2n->GetXaxis()->SetTitle("Threshold (keV)");
            grCd132p2n->GetYaxis()->SetTitle("P2n (/100 %)");

            grAg130p2n->SetTitle("P2n of Ag130");
            grAg130p2n->GetXaxis()->SetTitle("Threshold (keV)");
            grAg130p2n->GetYaxis()->SetTitle("P2n (/100 %)");




    cc->Divide(3,4);
    cc->cd(1);
    grSn136hl->Draw();
    cc->cd(2);
    grSn136p1n->Draw();
    cc->cd(3);
    grSn136p2n->Draw();

    cc->cd(4);
    grIn133hl->Draw();
    cc->cd(5);
    grIn133p1n->Draw();
    cc->cd(6);
    grIn133p2n->Draw();

    cc->cd(7);
    grCd132hl->Draw();
    cc->cd(8);
    grCd132p1n->Draw();
    cc->cd(9);
    grCd132p2n->Draw();

    cc->cd(10);
    grAg130hl->Draw();
    cc->cd(11);
    grAg130p1n->Draw();
    cc->cd(12);
    grAg130p2n->Draw();

}
