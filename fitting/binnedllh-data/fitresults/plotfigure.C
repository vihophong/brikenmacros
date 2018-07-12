/// \file
/// \ingroup tutorial_fit
/// \notebook
/// Combined (simultaneous) fit of two histogram with separate functions
/// and some common parameters
///
/// See http://root.cern.ch/phpBB3//viewtopic.php?f=3&t=11740#p50908
/// for a modified version working with Fumili or GSLMultiFit
///
/// N.B. this macro must be compiled with ACliC
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Lorenzo Moneta

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/PoissonLikelihoodFCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TFile.h"
#include <fstream>


void plotfigure(char* infile){

    TFile *file0 = TFile::Open(infile);

   TH1F*hdecay = (TH1F*)file0->Get("hdecay");
   TH1F*hdecay1n = (TH1F*)file0->Get("hdecay1n");
   TH1F*hdecay2n = (TH1F*)file0->Get("hdecay2n");

   TH1F*hdecaynuc1 = (TH1F*)file0->Get("hdecaynuc1");
   TH1F*hdecaynuc2 = (TH1F*)file0->Get("hdecaynuc2");
   TH1F*hdecaynuc3 = (TH1F*)file0->Get("hdecaynuc3");
   TH1F*hdecaynuc4 = (TH1F*)file0->Get("hdecaynuc4");
   TH1F*hdecaynuc5 = (TH1F*)file0->Get("hdecaynuc5");
   TH1F*hdecaynuc6 = (TH1F*)file0->Get("hdecaynuc6");
   TH1F*hdecaynuc9 = (TH1F*)file0->Get("hdecaynuc9");


   TH1F*h1ncomp1 = (TH1F*)file0->Get("h1ncomp1");
   TH1F*h1ncomp2 = (TH1F*)file0->Get("h1ncomp2");
   TH1F*h1ncomp3 = (TH1F*)file0->Get("h1ncomp3");

   TH1F*h2ncomp1 = (TH1F*)file0->Get("h2ncomp1");
   TH1F*h2ncomp2 = (TH1F*)file0->Get("h2ncomp2");
   TH1F*h2ncomp3 = (TH1F*)file0->Get("h2ncomp3");
   TH1F*h2ncomp4 = (TH1F*)file0->Get("h2ncomp4");



   TF1* fbBnuc1=(TF1*)file0->Get("fbBnuc1");
   TF1* fbBnuc2=(TF1*)file0->Get("fbBnuc2");
   TF1* fbBnuc3=(TF1*)file0->Get("fbBnuc3");
   TF1* fbBnuc4=(TF1*)file0->Get("fbBnuc4");
   TF1* fbBnuc5=(TF1*)file0->Get("fbBnuc5");
   TF1* fbBnuc6=(TF1*)file0->Get("fbBnuc6");
   TF1* fbBnuc9=(TF1*)file0->Get("fbBnuc9");

   TF1* fSBc1=(TF1*)file0->Get("fSBc1");
   TF1* fSBc2=(TF1*)file0->Get("fSBc2");
   TF1* fSBc3=(TF1*)file0->Get("fSBc3");


   TF1* fSB2c1=(TF1*)file0->Get("fSB2c1");
   TF1* fSB2c2=(TF1*)file0->Get("fSB2c2");
   TF1* fSB2c3=(TF1*)file0->Get("fSB2c3");
   TF1* fSB2c4=(TF1*)file0->Get("fSB2c4");


   TF1* fB=(TF1*)file0->Get("fB");
   TF1* fSB=(TF1*)file0->Get("fSB");
   TF1* fSB2=(TF1*)file0->Get("fSB2");


   gStyle->SetOptStat(0);
   TCanvas * c3 = new TCanvas("Simfit3","Simultaneous fit of 3 histograms",
                              10,10,800,600);

   c3->Divide(1,3);
   c3->cd(1);
   c3->cd(1)->SetLogy();
   c3->cd(1)->SetGrid();
   /*
   gPad->SetBottomMargin(0.001);
   gPad->SetTopMargin(0.01);
   gPad->SetRightMargin(0.01);
   */

   hdecay->Draw("PL");
   hdecay->SetMarkerStyle(20);
   hdecay->SetMarkerColor(1);
   hdecay->SetMarkerSize(0.7);
   hdecay->SetLineColor(1);

   hdecay->GetXaxis()->SetRangeUser(-1,10);
   hdecay->GetYaxis()->SetRangeUser(1,1e6);

   hdecaynuc1->Draw("PLsame");
   hdecaynuc2->Draw("PLsame");
   hdecaynuc3->Draw("PLsame");
   hdecaynuc4->Draw("PLsame");
   hdecaynuc5->Draw("PLsame");
   hdecaynuc6->Draw("PLsame");
   hdecaynuc9->Draw("PLsame");

   hdecaynuc1->SetMarkerStyle(20);
   hdecaynuc2->SetMarkerStyle(20);
   hdecaynuc3->SetMarkerStyle(20);
   hdecaynuc4->SetMarkerStyle(20);
   hdecaynuc5->SetMarkerStyle(20);
   hdecaynuc6->SetMarkerStyle(20);
   hdecaynuc9->SetMarkerStyle(20);

   hdecaynuc1->SetMarkerColor(1);
   hdecaynuc2->SetMarkerColor(1);
   hdecaynuc3->SetMarkerColor(1);
   hdecaynuc4->SetMarkerColor(1);
   hdecaynuc5->SetMarkerColor(1);
   hdecaynuc6->SetMarkerColor(1);
   hdecaynuc9->SetMarkerColor(1);

   hdecaynuc1->SetMarkerSize(0.5);
   hdecaynuc2->SetMarkerSize(0.5);
   hdecaynuc3->SetMarkerSize(0.5);
   hdecaynuc4->SetMarkerSize(0.5);
   hdecaynuc5->SetMarkerSize(0.5);
   hdecaynuc6->SetMarkerSize(0.5);
   hdecaynuc9->SetMarkerSize(0.5);

   hdecaynuc1->SetLineColor(2);
   hdecaynuc2->SetLineColor(2);
   hdecaynuc3->SetLineColor(2);
   hdecaynuc4->SetLineColor(2);
   hdecaynuc5->SetLineColor(2);
   hdecaynuc6->SetLineColor(2);
   hdecaynuc9->SetLineColor(2);

   fB->SetLineColor(3);
   fB->Draw("same");

   fbBnuc1->SetLineColor(4);
   fbBnuc2->SetLineColor(5);
   fbBnuc3->SetLineColor(6);
   fbBnuc4->SetLineColor(7);
   fbBnuc5->SetLineColor(9);
   fbBnuc6->SetLineColor(11);
   fbBnuc9->SetLineColor(12);

   fbBnuc1->Draw("same");
   fbBnuc2->Draw("same");
   fbBnuc3->Draw("same");
   fbBnuc4->Draw("same");
   fbBnuc5->Draw("same");
   fbBnuc6->Draw("same");
   fbBnuc9->Draw("same");


   c3->cd(2);
   c3->cd(2)->SetLogy();
   c3->cd(2)->SetGrid();
   /*
   gPad->SetBottomMargin(0.001);
   gPad->SetTopMargin(0.001);
   gPad->SetRightMargin(0.01);
   */

   hdecay1n->Draw("PL");
   hdecay1n->SetMarkerStyle(20);
   hdecay1n->SetMarkerColor(1);
   hdecay1n->SetMarkerSize(0.7);
   hdecay1n->SetLineColor(2);

   hdecay1n->GetXaxis()->SetRangeUser(-1,10);
   hdecay1n->GetYaxis()->SetRangeUser(1,5e5);

   h1ncomp1->Draw("PLsame");
   h1ncomp2->Draw("PLsame");
   h1ncomp3->Draw("PLsame");

   h1ncomp1->SetMarkerStyle(20);
   h1ncomp2->SetMarkerStyle(20);
   h1ncomp3->SetMarkerStyle(20);

   h1ncomp1->SetMarkerColor(1);
   h1ncomp2->SetMarkerColor(1);
   h1ncomp3->SetMarkerColor(1);

   h1ncomp1->SetMarkerSize(0.5);
   h1ncomp2->SetMarkerSize(0.5);
   h1ncomp3->SetMarkerSize(0.5);

   h1ncomp1->SetLineColor(2);
   h1ncomp2->SetLineColor(2);
   h1ncomp3->SetLineColor(2);

   fSB->SetLineColor(3);
   fSB->Draw("same");

   fSBc1->SetLineColor(4);
   fSBc2->SetLineColor(5);
   fSBc3->SetLineColor(6);
   fSBc1->Draw("same");
   fSBc2->Draw("same");
   fSBc3->Draw("same");

   c3->cd(3);
   c3->cd(3)->SetLogy();
   c3->cd(3)->SetGrid();
   /*
   gPad->SetBottomMargin(0.1);
   gPad->SetTopMargin(0.001);
   gPad->SetRightMargin(0.01);
   */

   hdecay2n->Draw("PL");
   hdecay2n->SetMarkerStyle(20);
   hdecay2n->SetMarkerColor(1);
   hdecay2n->SetMarkerSize(0.7);
   hdecay2n->SetLineColor(1);
   hdecay2n->GetXaxis()->SetRangeUser(-1,10);
   hdecay2n->GetYaxis()->SetRangeUser(1,1e5);

   h2ncomp1->Draw("PLsame");
   h2ncomp2->Draw("PLsame");
   h2ncomp3->Draw("PLsame");
   h2ncomp4->Draw("PLsame");

   h2ncomp1->SetMarkerStyle(20);
   h2ncomp2->SetMarkerStyle(20);
   h2ncomp3->SetMarkerStyle(20);
   h2ncomp4->SetMarkerStyle(20);

   h2ncomp1->SetMarkerColor(1);
   h2ncomp2->SetMarkerColor(1);
   h2ncomp3->SetMarkerColor(1);
   h2ncomp4->SetMarkerColor(1);

   h2ncomp1->SetMarkerSize(0.5);
   h2ncomp2->SetMarkerSize(0.5);
   h2ncomp3->SetMarkerSize(0.5);
   h2ncomp4->SetMarkerSize(0.5);

   h2ncomp1->SetLineColor(2);
   h2ncomp2->SetLineColor(2);
   h2ncomp3->SetLineColor(2);
   h2ncomp4->SetLineColor(2);

   fSB2->SetLineColor(3);
   fSB2->Draw("same");

   fSB2c1->SetLineColor(4);
   fSB2c2->SetLineColor(5);
   fSB2c3->SetLineColor(6);
   fSB2c4->SetLineColor(7);
   fSB2c1->Draw("same");
   fSB2c2->Draw("same");
   fSB2c3->Draw("same");
   fSB2c4->Draw("same");

}

