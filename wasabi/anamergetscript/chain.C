#include "TChain.h"
#include "TLatex.h"
#include "TString.h"
#include "TCutG.h"


TCutG* cutg1[1000];

void chainwithCut(char* listfile,char* pidfile,Int_t file_beg=0,Int_t file_end=38){
  char tempchar1[1000];
  sprintf(tempchar1,"tree");
  TChain* ch = new TChain(tempchar1);
  std::ifstream ifs(listfile);
  string filelist[1000];

  Int_t nfiles=0;
  while (!ifs.eof()){
      ifs>>filelist[nfiles];
      cout<<filelist[nfiles]<<endl;
      nfiles++;
  }
  nfiles=nfiles-1;
  cout<<"There are "<<nfiles<<" files in total!"<<endl;

  for (Int_t i=0;i<nfiles;i++){
      char tempchar2[1000];
      sprintf(tempchar2,"%s/tree",filelist[i].c_str());
      ch->Add(tempchar2);
  }

  std::ifstream ifspid(pidfile);
  Int_t nri=0;
  Int_t nbinszet,nbinsaoq;
  Double_t zetrange[2];
  Double_t aoqrange[2];
  ifspid>>nbinsaoq>>aoqrange[0]>>aoqrange[1]>>nbinszet>>zetrange[0]>>zetrange[1];
  ifspid>>nri;

  Int_t enablepid[nri];
  Int_t enablepid2[nri];

  TString nameri[nri];
  TString latexnametri[nri];
  Double_t parmsri[nri][7];
  TString tempriname,tempria;
  TLatex* pidtag[nri];

  Double_t halflife[nri];

  for (Int_t i=0;i<nri;i++){
      ifspid>>enablepid[i]>>enablepid2[i]>>tempria>>tempriname>>halflife[i];
      for(Int_t j=0;j<7;j++) ifspid>>parmsri[i][j];
      nameri[i]=tempriname+tempria;
      latexnametri[i]=TString("^{")+tempria+TString("}"+tempriname);
      cout<<nameri[i]<<"\t"<<latexnametri[i]<<"\t"<<halflife[i];
      for(Int_t j=0;j<7;j++) cout<<"\t"<<parmsri[i][j];
      cout<<endl;

      pidtag[i]=new TLatex(parmsri[i][0],parmsri[i][1]+0.2,latexnametri[i]);
      pidtag[i]->SetTextSize(0.025);
      pidtag[i]->SetTextColor(2);
  }

  Int_t ncutpts=20;// number of cut points
  for (Int_t i=0;i<nri;i++){
      cutg1[i]=new TCutG(nameri[i],ncutpts);
      cutg1[i]->SetTitle(nameri[i]);
      cutg1[i]->SetVarX("aoq");
      cutg1[i]->SetVarY("zet");

      Double_t theta=0;
      Double_t x1=parmsri[i][0];
      Double_t y1=parmsri[i][1];
      Double_t r1=parmsri[i][2];
      Double_t r2=parmsri[i][3];

      Double_t phi1 = TMath::Min(0,360);
      Double_t phi2 = TMath::Max(0,360);
      Double_t kPI = 3.14159265358979323846;
      Double_t angle,dx,dy;

      Double_t x[ncutpts], y[ncutpts];
      Double_t dphi = (phi2-phi1)*kPI/(180*ncutpts);
      Double_t ct   = TMath::Cos(kPI*theta/180);
      Double_t st   = TMath::Sin(kPI*theta/180);

      for (Int_t j=0;j<ncutpts;j++) {
         angle = phi1*kPI/180 + Double_t(j)*dphi;
         dx    = r1*TMath::Cos(angle);
         dy    = r2*TMath::Sin(angle);
         x[j]  = x1 + dx*ct - dy*st;
         y[j]  = y1 + dx*st + dy*ct;
         //cout<<x[j]<<"<"<<y[j]<<endl;
         cutg1[i]->SetPoint(j,x[j],y[j]);
      }
      cutg1[i]->SetPoint(ncutpts,x[0],y[0]);
  }
  TCanvas* c1=new TCanvas("pid","pid",900,700);
  c1->cd();

  TH2F* pidall;
  TH2F* pidi[nfiles];
  for (Int_t i=file_beg;i<file_end+1;i++){
      TFile *f1 = TFile::Open(filelist[i].c_str());
      pidi[i]=(TH2F*) f1->Get("pidall");
/*
      if (i==0){
              for (Int_t j=0;j<nri;j++){
                 cutg[j]=(TCutG* )f1->Get(nameri[j]);
                 cout<<nameri[j]<<endl;
              }
      }
*/
      if (i==0) pidall=(TH2F*) pidi[i]->Clone();
      else pidall->Add(pidi[i]);

      cout<<i<<endl;
  }
  pidall->Draw("colz");
  for (Int_t i=0;i<nri;i++){
       cutg1[i]->Draw("same");
       pidtag[i]->Draw("same");
  }
  TCanvas* c2=new TCanvas("plot","plot",900,700);
  c2->cd();

}

void PlotPid(char* listfile,char* pidfile){
    std::ifstream ifspid(pidfile);
    Int_t nri=0;
    Int_t nbinszet,nbinsaoq;
    Double_t zetrange[2];
    Double_t aoqrange[2];
    ifspid>>nbinsaoq>>aoqrange[0]>>aoqrange[1]>>nbinszet>>zetrange[0]>>zetrange[1];
    ifspid>>nri;

    Int_t enablepid[nri];
    Int_t enablepid2[nri];

    TString nameri[nri];
    TString latexnametri[nri];
    Double_t parmsri[nri][7];
    TString tempriname,tempria;
    TCutG* cutg[nri];
    TLatex* pidtag[nri];

    Double_t halflife[nri];

    for (Int_t i=0;i<nri;i++){
        ifspid>>enablepid[i]>>enablepid2[i]>>tempria>>tempriname>>halflife[i];
        for(Int_t j=0;j<7;j++) ifspid>>parmsri[i][j];
        nameri[i]=tempriname+tempria;
        latexnametri[i]=TString("^{")+tempria+TString("}"+tempriname);
        cout<<nameri[i]<<"\t"<<latexnametri[i]<<"\t"<<halflife[i];
        for(Int_t j=0;j<7;j++) cout<<"\t"<<parmsri[i][j];
        cout<<endl;

        pidtag[i]=new TLatex(parmsri[i][0],parmsri[i][1]+0.2,latexnametri[i]);
        pidtag[i]->SetTextSize(0.025);
        pidtag[i]->SetTextColor(2);
    }

  Int_t ncutpts=20;// number of cut points
  for (Int_t i=0;i<nri;i++){      
      cutg[i]=new TCutG(nameri[i],ncutpts);      
      cutg[i]->SetTitle(nameri[i]);
      cutg[i]->SetVarX("decayaoq");
      cutg[i]->SetVarY("decay.zet");

      Double_t theta=0;
      Double_t x1=parmsri[i][0];
      Double_t y1=parmsri[i][1];
      Double_t r1=parmsri[i][2];
      Double_t r2=parmsri[i][3];

      Double_t phi1 = TMath::Min(0,360);
      Double_t phi2 = TMath::Max(0,360);
      Double_t kPI = 3.14159265358979323846;
      Double_t angle,dx,dy;

      Double_t x[ncutpts], y[ncutpts];
      Double_t dphi = (phi2-phi1)*kPI/(180*ncutpts);
      Double_t ct   = TMath::Cos(kPI*theta/180);
      Double_t st   = TMath::Sin(kPI*theta/180);

      for (Int_t j=0;j<ncutpts;j++) {
         angle = phi1*kPI/180 + Double_t(j)*dphi;
         dx    = r1*TMath::Cos(angle);
         dy    = r2*TMath::Sin(angle);
         x[j]  = x1 + dx*ct - dy*st;
         y[j]  = y1 + dx*st + dy*ct;
         //cout<<x[j]<<"<"<<y[j]<<endl;
         cutg[i]->SetPoint(j,x[j],y[j]);
      }
      cutg[i]->SetPoint(ncutpts,x[0],y[0]);
  }


    std::ifstream ifs(listfile);
    string filelist[1000];
    Int_t nfiles=0;
    while (!ifs.eof()){
        ifs>>filelist[nfiles];
        cout<<filelist[nfiles]<<endl;
        nfiles++;
    }
    nfiles=nfiles-1;
    cout<<"There are "<<nfiles<<" files in total!"<<endl;
    TH2F* pidall;
    TH2F* pidi[nfiles];
    for (Int_t i=0;i<nfiles;i++){
	TFile *f1 = TFile::Open(filelist[i].c_str());        
        pidi[i]=(TH2F*) f1->Get("pidall");
/*
	if (i==0){
        	for (Int_t j=0;j<nri;j++){
		   cutg[j]=(TCutG* )f1->Get(nameri[j]);
		   cout<<nameri[j]<<endl;
		} 
	}   
*/
        if (i==0) pidall=(TH2F*) pidi[i]->Clone();
        else pidall->Add(pidi[i]);

	cout<<i<<endl;
    }
    TCanvas* c1=new TCanvas("hprofile","hprofile",900,700);
    pidall->Draw("colz");
    for (Int_t i=0;i<nri;i++){
	if (enablepid2[i]){ 
        cutg[i]->Draw("same");
	pidtag[i]->Draw("same");
	}
    }
}


void Internalchaineach(char* listfile,char* pid){
  char tempchar1[1000];
  sprintf(tempchar1,"tree%s",pid);
  TChain* ch = new TChain(tempchar1);
  std::ifstream ifs(listfile);
  string filelist[1000];

  Int_t nfiles=0;
  while (!ifs.eof()){
      ifs>>filelist[nfiles];
      cout<<filelist[nfiles]<<endl;
      nfiles++;
  }
  nfiles=nfiles-1;
  cout<<"There are "<<nfiles<<" files in total!"<<endl;

  for (Int_t i=0;i<nfiles;i++){
      char tempchar2[1000];
      sprintf(tempchar2,"%s/tree%s",filelist[i].c_str(),pid);
      ch->Add(tempchar2);
  }
}

void Internalchain(char* listfile,char* pidfile){
    std::ifstream ifspid(pidfile);
    Int_t nri=0;
    Int_t nbinszet,nbinsaoq;
    Double_t zetrange[2];
    Double_t aoqrange[2];
    ifspid>>nbinsaoq>>aoqrange[0]>>aoqrange[1]>>nbinszet>>zetrange[0]>>zetrange[1];
    ifspid>>nri;

    Int_t enablepid[nri];
    Int_t enablepid2[nri];

    TString nameri[nri];
    TString latexnametri[nri];
    Double_t parmsri[nri][7];
    TString tempriname,tempria;
    TCutG* cutg[nri];
    TLatex* pidtag[nri];

    Double_t halflife[nri];

    for (Int_t i=0;i<nri;i++){
        ifspid>>enablepid[i]>>enablepid2[i]>>tempria>>tempriname>>halflife[i];
        for(Int_t j=0;j<7;j++) ifspid>>parmsri[i][j];
        nameri[i]=tempriname+tempria;
        latexnametri[i]=TString("^{")+tempria+TString("}"+tempriname);
        //cout<<nameri[i]<<"\t"<<latexnametri[i]<<"\t"<<halflife[i];
        for(Int_t j=0;j<7;j++) cout<<"\t"<<parmsri[i][j];
        //cout<<endl;
        pidtag[i]=new TLatex(parmsri[i][0],parmsri[i][1]+0.2,latexnametri[i]);
        pidtag[i]->SetTextSize(0.025);
        pidtag[i]->SetTextColor(2);
        Internalchaineach(listfile,(char*)nameri[i].Data());
    }
    Internalchaineach(listfile,(char*)TString("").Data());
}


void InternalplotPIDbetaion(char* listfile,char* pidfile){;
    std::ifstream ifspid(pidfile);
    Int_t nri=0;
    Int_t nbinszet,nbinsaoq;
    Double_t zetrange[2];
    Double_t aoqrange[2];
    ifspid>>nbinsaoq>>aoqrange[0]>>aoqrange[1]>>nbinszet>>zetrange[0]>>zetrange[1];
    ifspid>>nri;

    Int_t enablepid[nri];
    Int_t enablepid2[nri];

    TString nameri[nri];
    TString latexnametri[nri];
    Double_t parmsri[nri][7];
    TString tempriname,tempria;
    TCutG* cutg[nri];
    TLatex* pidtag[nri];

    Double_t halflife[nri];

    for (Int_t i=0;i<nri;i++){
        ifspid>>enablepid[i]>>enablepid2[i]>>tempria>>tempriname>>halflife[i];
        for(Int_t j=0;j<7;j++) ifspid>>parmsri[i][j];
        nameri[i]=tempriname+tempria;
        latexnametri[i]=TString("^{")+tempria+TString("}"+tempriname);
        cout<<nameri[i]<<"\t"<<latexnametri[i]<<"\t"<<halflife[i];
        for(Int_t j=0;j<7;j++) cout<<"\t"<<parmsri[i][j];
        cout<<endl;

        pidtag[i]=new TLatex(parmsri[i][0],parmsri[i][1]+0.2,latexnametri[i]);
        pidtag[i]->SetTextSize(0.025);
        pidtag[i]->SetTextColor(2);
    }
    std::ifstream ifs(listfile);
    string filelist[1000];
    Int_t nfiles=0;
    while (!ifs.eof()){
        ifs>>filelist[nfiles];
        cout<<filelist[nfiles]<<endl;
        nfiles++;
    }
    TFile* f0=new TFile(filelist[0].c_str());
    for (Int_t i=0;i<nri;i++){
         cutg[i]=(TCutG* )f0->Get(nameri[i]);
    }
    TChain* ch = new TChain("tree");
    nfiles=nfiles-1;
    cout<<"There are "<<nfiles<<" files in total!"<<endl;

    for (Int_t i=0;i<nfiles;i++){
        char tempchar2[1000];
        sprintf(tempchar2,"%s/tree",filelist[i].c_str());
        ch->Add(tempchar2);
    }
    TCanvas* c1=new TCanvas("pid","pid",900,700);
    ch->Draw(Form("zet:aoq>>h1(%d,%f,%f,%d,%f,%f)",nbinsaoq,aoqrange[0],aoqrange[1],nbinszet,zetrange[0],zetrange[1]),"zet>0&&aoq>0&&ion.fz!=5","colz");
    for (Int_t i=0;i<nri;i++){
      //cutg[i]->Draw("same");
         pidtag[i]->Draw("same");
    }
}

void InternalPlotProfile(char* listfile,char* riname){
    std::ifstream ifs(listfile);
    string filelist[1000];
    Int_t nfiles=0;
    while (!ifs.eof()){
        ifs>>filelist[nfiles];
        cout<<filelist[nfiles]<<endl;
        nfiles++;
    }
    nfiles=nfiles-1;
    cout<<"There are "<<nfiles<<" files in total!"<<endl;
    char tempchar2[1000];
    sprintf(tempchar2,"%s_depth",riname);
    TH1F* hprofileall;
    TH1F* hprofile[nfiles];
    for (Int_t i=0;i<nfiles;i++){
        TFile *f1 = TFile::Open(filelist[i].c_str());
        hprofile[i]=(TH1F*)f1->Get(tempchar2);
        //f1->Close();
        if (i==0) hprofileall=(TH1F*) hprofile[i]->Clone();
        else hprofileall->Add(hprofile[i]);
    }
    TCanvas* c1=new TCanvas("hprofile","hprofile",900,700);
/*
    c1->Divide(4,1);
    c1->cd(1);
    hprofile[0]->Draw();
    c1->cd(2);
    hprofile[1]->Draw();
    c1->cd(3);
    hprofile[2]->Draw();
    c1->cd(4);
*/
    hprofileall->Draw();
    Int_t nimpWasabireal02=0;
    Int_t nimpWasabireal=0;
    for (int i=0;i<4;i++) {
        cout<<"dssd "<<i+1<<" : nimplant"<<riname<<" = "<<hprofileall->GetBinContent(i+1)<<endl;
        if (i<3) nimpWasabireal02+=hprofileall->GetBinContent(i+1);
        nimpWasabireal+=hprofileall->GetBinContent(i+1);
    }
    cout<<"YSOimplant = "<<hprofileall->GetBinContent(5)<<endl;
    cout<<"Nimplant in Wasabi (layer 1-3) = "<<nimpWasabireal02<<endl;
    cout<<"Nimplant in Wasabi (layer 1-4) = "<<nimpWasabireal-hprofileall->GetBinContent(5)<<endl;
    cout<<"Ntotal implantation"<<nimpWasabireal<<endl;
}


void Internalgetnimplant(Int_t imparr[],char* listfile,char* riname){
    std::ifstream ifs(listfile);
    string filelist[1000];
    Int_t nfiles=0;
    while (!ifs.eof()){
        ifs>>filelist[nfiles];
        //cout<<filelist[nfiles]<<endl;
        nfiles++;
    }
    nfiles=nfiles-1;
    //cout<<"There are "<<nfiles<<" files in total!"<<endl;
    char tempchar2[1000];
    sprintf(tempchar2,"%s_depth",riname);
    TH1F* hprofileall;
    TH1F* hprofile[nfiles];
    TH1::AddDirectory(kFALSE);
     TFile *f1[nfiles];
    for (Int_t i=0;i<nfiles;i++){
        f1[i]= TFile::Open(filelist[i].c_str());
        hprofile[i]=(TH1F*)f1[i]->Get(tempchar2);
        hprofile[i]->AddDirectory(kFALSE);
        //f1->Close();
        if (i==0) hprofileall=(TH1F*) hprofile[i]->Clone();
        else hprofileall->Add(hprofile[i]);
    }
    for (int i=0;i<5;i++) {
        imparr[i]=hprofileall->GetBinContent(i+1);
    }
    for (Int_t i=0;i<nfiles;i++){
        f1[i]->Close();
        delete f1[i];
    }

}

void InternalMakeImplantTable(char* listfile,char* outtextfile, char* pidfile){
    std::ifstream ifspid(pidfile);
    Int_t nri=0;
    Int_t nbinszet,nbinsaoq;
    Double_t zetrange[2];
    Double_t aoqrange[2];
    ifspid>>nbinsaoq>>aoqrange[0]>>aoqrange[1]>>nbinszet>>zetrange[0]>>zetrange[1];
    ifspid>>nri;

    Int_t enablepid[nri];
    Int_t enablepid2[nri];

    TString nameri[nri];
    TString latexnametri[nri];
    Double_t parmsri[nri][7];
    TString tempriname,tempria;

    Double_t halflife[nri];


    for (Int_t i=0;i<nri;i++){
        ifspid>>enablepid[i]>>enablepid2[i]>>tempria>>tempriname>>halflife[i];
        for(Int_t j=0;j<7;j++) ifspid>>parmsri[i][j];
        nameri[i]=tempriname+tempria;
        latexnametri[i]=TString("^{")+tempria+TString("}"+tempriname);
        cout<<nameri[i]<<"\t"<<latexnametri[i]<<"\t"<<halflife[i];
        for(Int_t j=0;j<7;j++) cout<<"\t"<<parmsri[i][j];
        cout<<endl;
        /*
        pidtag[i]=new TLatex(parmsri[i][0],parmsri[i][1]+0.2,latexnametri[i]);
        pidtag[i]->SetTextSize(0.025);
        pidtag[i]->SetTextColor(2);
        */
    }
    std::ofstream ofs(outtextfile);
    for (Int_t i=0;i<nri;i++){
        if (enablepid2[i]){
            Int_t implantarr[5];
            Internalgetnimplant(implantarr,listfile,(char*)nameri[i].Data());
            ofs<<nameri[i]<<"\t"<<implantarr[0]<<"\t"<<implantarr[1]<<"\t"<<implantarr[2]<<"\t"<<implantarr[3]-implantarr[4]<<"\t"<<implantarr[4]<<"\t"<<implantarr[0]+implantarr[1]+implantarr[2]+implantarr[3]<<endl;
            cout<<nameri[i]<<"\t"<<implantarr[0]<<"\t"<<implantarr[1]<<"\t"<<implantarr[2]<<"\t"<<implantarr[3]-implantarr[4]<<"\t"<<implantarr[4]<<"\t"<<implantarr[0]+implantarr[1]+implantarr[2]+implantarr[3]<<endl;
        }
    }
}

