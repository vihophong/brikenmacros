#include "TChain.h"
#include "TLatex.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
void plotPid(char* listfile, char* pidfile, char* outfile){

  //! get pid cut
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
           cutg[i]->SetPoint(j,x[j],y[j]);
        }
        cutg[i]->SetPoint(ncutpts,x[0],y[0]);
    }
    //! get pid

    char tempchar1[1000];
    sprintf(tempchar1,"treeimpall");
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
        sprintf(tempchar2,"%s/treeimpall",filelist[i].c_str());
        ch->Add(tempchar2);
    }

    TCanvas* c1=new TCanvas("pid","pid",900,700);
    ch->Draw(Form("implant.beam.zet:implant.beam.aoq>>h1(%d,%f,%f,%d,%f,%f)",nbinsaoq,aoqrange[0],aoqrange[1],nbinszet,zetrange[0],zetrange[1]),"implant.beam.zet>0&&implant.beam.aoq>0","colz");
    TH2F* pid=(TH2F*) gDirectory->Get("h1");
    for (Int_t i=0;i<nri;i++){
         cutg[i]->Draw("same");
         pidtag[i]->Draw("same");
    }

    TFile* fout=new TFile(outfile,"recreate");
    pid->Write();
    fout->Close();
}

void slicepid(char* pidfile, char* infile,char* fitparms,Double_t Zcenter,Double_t width,Double_t low=2.65, Double_t high=2.81, Int_t idpidlow=0,Int_t idpidhigh=0)
{
    //! get pid cut
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
          cout<<i<<"\t"<<nameri[i]<<"\t"<<latexnametri[i]<<"\t"<<halflife[i];
          for(Int_t j=0;j<7;j++) cout<<"\t"<<parmsri[i][j];
          cout<<endl;

          pidtag[i]=new TLatex(parmsri[i][0],parmsri[i][1]+0.2,latexnametri[i]);
          pidtag[i]->SetTextSize(0.045);
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
             cutg[i]->SetPoint(j,x[j],y[j]);
          }
          cutg[i]->SetPoint(ncutpts,x[0],y[0]);
      }

    TFile* fin=TFile::Open(infile);
    TH2F* pid=(TH2F*) fin->Get("h1");
    Int_t startbin=pid->GetYaxis()->FindBin(Zcenter-width/2);
    Int_t stopbin=pid->GetYaxis()->FindBin(Zcenter+width/2);
    cout<<startbin<<"-"<<stopbin<<endl;
    TH1F* hproj=(TH1F*) pid->ProjectionX("prj",startbin,stopbin);
    TCanvas* c1=new TCanvas("c1","c1",900,1200);
    c1->Divide(1,2);
    c1->cd(1)->SetLogz();
    TLine* l1=new TLine(pid->GetXaxis()->GetXmin(), Zcenter-width/2,pid->GetXaxis()->GetXmax(),Zcenter-width/2);
    TLine* l2=new TLine(pid->GetXaxis()->GetXmin(), Zcenter+width/2,pid->GetXaxis()->GetXmax(),Zcenter+width/2);
    l1->SetLineColor(2);
    l2->SetLineColor(2);
    pid->Draw("colz");
    l1->Draw("same");
    l2->Draw("same");


    if (idpidlow==0&&idpidhigh==0)
        idpidhigh=nri;

    for (Int_t i=0;i<nri;i++){
        if (i>=idpidlow&&i<idpidhigh){
            cutg[i]->Draw("same");
            pidtag[i]->Draw("same");
        }
    }

    c1->cd(2)->SetLogy();
    //fitting

    std::ifstream ifs(fitparms);
    Double_t meanval[1000];
    Double_t meanvalvar[1000];
    Double_t sigmavar[1000];

    Int_t npars=0;
    while (!ifs.eof()){
        ifs>>meanval[npars]>>meanvalvar[npars]>>sigmavar[npars];
        cout<<meanval[npars]<<endl;
        npars++;
    }
    npars--;
    cout<<npars<<" lines"<<endl;

    Double_t parms[1000];
    Double_t parmsmin[1000];
    Double_t parmsmax[1000];
    Double_t cconst=1e4;
    Double_t cconstmin=0;
    Double_t cconstmax=1e9;
    Double_t sigma=1.34352e-03;

    TString fdef("");

    for (Int_t i=0;i<npars;i++){
        parms[i*3]=cconst;
        parms[i*3+1]=meanval[i];
        parms[i*3+2]=sigma;
        parmsmin[i*3]=cconstmin;
        parmsmin[i*3+1]=meanval[i]-meanvalvar[i];
        parmsmin[i*3+2]=sigma-sigmavar[i];
        parmsmax[i*3]=cconstmax;
        parmsmax[i*3+1]=meanval[i]+meanvalvar[i];
        parmsmax[i*3+2]=sigma+sigmavar[i];
        fdef=fdef+TString(Form("gaus(%d)+",3*i));
    }
    fdef=fdef+TString("0");

    TF1* f1=new TF1("f1",(char*)fdef.Data(),low,high);
    for (Int_t i=0;i<npars*3;i++) f1->SetParLimits(i,parmsmin[i],parmsmax[i]);

    f1->SetParameters(parms);

    f1->SetNpx(2000);
    f1->SetLineColor(0);
    f1->SetLineWidth(0);
    //f1->Draw();
    hproj->Fit("f1","R");

    f1->FixParameter(36,8.39175e+00);
    f1->FixParameter(37,2.77539e+00);
    f1->FixParameter(38,-1.75946e-03 );
    f1->FixParameter(39,2.65119e+00);
    f1->FixParameter(40,2.78183e+00 );
    f1->FixParameter(41,1.57427e-03);
    TF1* fdecom[npars];

    Int_t itemp=0;
    Int_t tempp[]={6,8,10,12};
    for (Int_t i=0;i<npars;i++){
        fdecom[i]=new TF1(Form("fi%d",i),"gaus",low,high);
        fdecom[i]->SetParameter(0,f1->GetParameter(i*3));
        fdecom[i]->SetParameter(1,f1->GetParameter(i*3+1));
        fdecom[i]->SetParameter(2,f1->GetParameter(i*3+2));
        fdecom[i]->SetLineColor(3);

        if (i==6||i==8||i==10||i==12) {
            fdecom[i]->SetLineColor(2);
            cout<<nameri[86+itemp]<<"\t"<<fdecom[tempp[itemp]]->Integral(parmsri[86+itemp][0]-parmsri[86+itemp][2],high)<<"\t"<<parmsri[86+itemp][0]-parmsri[86+itemp][2]<<"\t"<<high<<endl;
            itemp++;
        }
    }

    for (Int_t i=0;i<npars;i++) fdecom[i]->Draw("same");


    //! Draw line
    TLine* sepll[500];
    TLine* seplr[500];

    Double_t maxy=hproj->GetMaximum();
    Int_t ni=0;
    for (Int_t i=idpidlow;i<idpidhigh;i++){
        sepll[ni]=new TLine(parmsri[i][0]-parmsri[i][2],0,parmsri[i][0]-parmsri[i][2],maxy);
        sepll[ni]->Draw("same");
        seplr[ni]=new TLine(parmsri[i][0]+parmsri[i][2],0,parmsri[i][0]+parmsri[i][2],maxy);
        seplr[ni]->Draw("same");
        ni++;
    }

    cout<<"\n\n"<<endl;
    for (Int_t i=0;i<npars;i++) cout<<i<<"\t"<<meanval[i]<<"\t"<<f1->GetParameter(i*3+0)<<"\t"<<f1->GetParameter(i*3+1)<<"\t"<<f1->GetParameter(i*3+2)<<"\t"<<f1->GetParameter(i*3+0)*f1->GetParameter(i*3+2)*sqrt(2*TMath::Pi())<<"\t"<<fdecom[i]->Integral(low,high)<<endl;

}


