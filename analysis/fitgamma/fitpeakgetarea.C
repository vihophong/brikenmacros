#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TSpectrum.h"
#include <iostream>
using namespace std;


void fitpeaksingle(TH1F* hist,Double_t low, Double_t high,Double_t bkglvl,Double_t sigma=2.)
{
    TString fdef("pol1(0)+gaus(2)");
    TF1* ffinal=new TF1("ffinal",(char*)fdef.Data(),low,high);
    TF1* fgaus=new TF1("fgaus","gaus",low,high);

    ffinal->SetParameter(0,bkglvl);
    ffinal->SetParameter(1,0.);

    //hist->SetAxisRange(low,high);
    TSpectrum* sp=new TSpectrum();
    sp->Search(hist);
    Double_t * xpeaks=sp->GetPositionX();

    Double_t peakcenter=xpeaks[0];
    Double_t peakheight=hist->GetBinContent(hist->GetXaxis()->FindBin(xpeaks[0]))-bkglvl;

    cout<<xpeaks[0]<<endl;

    ffinal->SetParameter(2,peakheight);// peak heigh
    ffinal->SetParameter(3,peakcenter);// centroid
    ffinal->SetParameter(4,sigma);// sigma

    hist->Fit("ffinal","LQR");
    cout<<"final fit result"<<endl;
    for (Int_t i=0;i<5;i++) {
        cout<<i<<"\t"<<ffinal->GetParameter(i)<<"\t"<<ffinal->GetParError(i)<<endl;
    }
    Double_t peakarea=ffinal->GetParameter(2)*ffinal->GetParameter(4)*sqrt(2*TMath::Pi());
    Double_t peakareaerr=peakarea*sqrt(pow((ffinal->GetParError(2)/ffinal->GetParameter(2)),2.) + pow((ffinal->GetParError(4)/ffinal->GetParameter(4)),2));

    fgaus->SetParameter(0,ffinal->GetParameter(2));
    fgaus->SetParameter(1,ffinal->GetParameter(3));
    fgaus->SetParameter(2,ffinal->GetParameter(4));
    cout<<"Chi-square/NDF = "<<fgaus->GetChisquare()/fgaus->GetNDF()<<endl;
    Double_t peakarea2=fgaus->Integral(low,high);
    cout<<"peak area="<<peakarea<<"\t error="<<peakareaerr<<" rel. err="<<peakareaerr/peakarea<<endl;

    cout<<peakarea<<"\t"<<peakareaerr<<"\t"<<peakareaerr/peakarea<<endl;
    cout<<peakarea<<"\t"<<peakareaerr<<"\t"<<peakareaerr/peakarea<<endl;

}

void fitpeakdouble(TH1F* hist,Double_t low, Double_t high,Double_t bkglvl,Double_t sigma=2.)
{
    TString fdef("pol1(0)+gaus(2)+gaus(5)");
    TF1* ffinal=new TF1("ffinal",(char*)fdef.Data(),low,high);

    ffinal->SetParameter(0,bkglvl);
    ffinal->SetParameter(1,0.);

    //hist->SetAxisRange(low,high);
    TSpectrum* sp=new TSpectrum();
    sp->Search(hist);
    Double_t * xpeaks=sp->GetPositionX();

    Double_t peakcenter[2]={xpeaks[0],xpeaks[1]};

    Double_t peakheight[2]={hist->GetBinContent(hist->GetXaxis()->FindBin(xpeaks[0]))-bkglvl,hist->GetBinContent(hist->GetXaxis()->FindBin(xpeaks[1]))-bkglvl};

    ffinal->SetParameter(2,peakheight[0]);// peak heigh
    ffinal->SetParameter(3,peakcenter[0]);// centroid
    ffinal->SetParameter(4,sigma);// sigma
    ffinal->SetParameter(5,peakheight[1]);// peak heigh
    ffinal->SetParameter(6,peakcenter[1]);// centroid
    ffinal->SetParameter(7,sigma);// sigma

    hist->Fit("ffinal","LQR");
    cout<<"final fit result"<<endl;
    for (Int_t i=0;i<8;i++) {
        cout<<i<<"\t"<<ffinal->GetParameter(i)<<endl;
    }

    Double_t peakarea1=ffinal->GetParameter(2)*ffinal->GetParameter(4)*sqrt(2*TMath::Pi());
    Double_t peakareaerr1=peakarea1*sqrt(pow((ffinal->GetParError(2)/ffinal->GetParameter(2)),2.) + pow((ffinal->GetParError(4)/ffinal->GetParameter(4)),2));
    Double_t peakarea2=ffinal->GetParameter(5)*ffinal->GetParameter(7)*sqrt(2*TMath::Pi());
    Double_t peakareaerr2=peakarea2*sqrt(pow((ffinal->GetParError(5)/ffinal->GetParameter(5)),2.) + pow((ffinal->GetParError(7)/ffinal->GetParameter(7)),2));

    cout<<"peak at "<<ffinal->GetParameter(3)<<" area= "<<peakarea1<<" error= "<<peakareaerr1<<" rel. err="<<peakareaerr1/peakarea1<<endl;
    cout<<"peak at "<<ffinal->GetParameter(6)<<" area= "<<peakarea2<<" error= "<<peakareaerr2<<" rel. err="<<peakareaerr2/peakarea2<<endl;

    cout<<ffinal->GetParameter(3)<<"\t"<<peakarea1<<"\t"<<peakareaerr1<<"\t"<<peakareaerr1/peakarea1<<endl;
    cout<<ffinal->GetParameter(6)<<"\t"<<peakarea2<<"\t"<<peakareaerr2<<"\t"<<peakareaerr2/peakarea2<<endl;
}


void fitpeaksingle(TH1F* hist,Double_t low, Double_t high,Double_t bkglvl,Double_t peak1,Double_t sigma=2.)
{
    TString fdef("pol1(0)+gaus(2)");
    TF1* ffinal=new TF1("ffinal",(char*)fdef.Data(),low,high);
    TF1* fgaus=new TF1("fgaus","gaus",low,high);

    ffinal->SetParameter(0,bkglvl);
    ffinal->SetParameter(1,0.);

    Double_t peakcenter=peak1;
    Double_t peakheight=hist->GetBinContent(hist->GetXaxis()->FindBin(peak1))-bkglvl;

    ffinal->SetParameter(2,peakheight);// peak heigh
    ffinal->SetParameter(3,peakcenter);// centroid
    ffinal->SetParameter(4,sigma);// sigma

    hist->Fit("ffinal","LQR");
    cout<<"final fit result"<<endl;
    for (Int_t i=0;i<5;i++) {
        cout<<i<<"\t"<<ffinal->GetParameter(i)<<"\t"<<ffinal->GetParError(i)<<endl;
    }
    Double_t peakarea=ffinal->GetParameter(2)*ffinal->GetParameter(4)*sqrt(2*TMath::Pi());
    Double_t peakareaerr=peakarea*sqrt(pow((ffinal->GetParError(2)/ffinal->GetParameter(2)),2.) + pow((ffinal->GetParError(4)/ffinal->GetParameter(4)),2));

    fgaus->SetParameter(0,ffinal->GetParameter(2));
    fgaus->SetParameter(1,ffinal->GetParameter(3));
    fgaus->SetParameter(2,ffinal->GetParameter(4));
    Double_t peakarea2=fgaus->Integral(low,high);
    cout<<"peak area="<<peakarea<<"\t error="<<peakareaerr<<" rel. err="<<peakareaerr/peakarea<<endl;

    cout<<ffinal->GetParameter(3)<<"\t"<<peakarea<<"\t"<<peakareaerr<<"\t"<<peakareaerr/peakarea<<endl;
    cout<<ffinal->GetParameter(3)<<"\t"<<peakarea<<"\t"<<peakareaerr<<"\t"<<peakareaerr/peakarea<<endl;

}

void fitpeakdouble(TH1F* hist,Double_t low, Double_t high,Double_t bkglvl,Double_t peak1, Double_t peak2,Double_t sigma=2.)
{
    TString fdef("pol1(0)+gaus(2)+gaus(5)");
    TF1* ffinal=new TF1("ffinal",(char*)fdef.Data(),low,high);

    ffinal->SetParameter(0,bkglvl);
    ffinal->SetParameter(1,0.);

    Double_t peakcenter[2]={peak1,peak2};

    Double_t peakheight[2]={hist->GetBinContent(hist->GetXaxis()->FindBin(peak2))-bkglvl,hist->GetBinContent(hist->GetXaxis()->FindBin(peak2))-bkglvl};

    ffinal->SetParameter(2,peakheight[0]);// peak heigh
    ffinal->SetParameter(3,peakcenter[0]);// centroid
    ffinal->SetParameter(4,sigma);// sigma
    ffinal->SetParameter(5,peakheight[1]);// peak heigh
    ffinal->SetParameter(6,peakcenter[1]);// centroid
    ffinal->SetParameter(7,sigma);// sigma

    hist->Fit("ffinal","LQR");
    cout<<"final fit result"<<endl;
    for (Int_t i=0;i<8;i++) {
        cout<<i<<"\t"<<ffinal->GetParameter(i)<<endl;
    }

    Double_t peakarea1=ffinal->GetParameter(2)*ffinal->GetParameter(4)*sqrt(2*TMath::Pi());
    Double_t peakareaerr1=peakarea1*sqrt(pow((ffinal->GetParError(2)/ffinal->GetParameter(2)),2.) + pow((ffinal->GetParError(4)/ffinal->GetParameter(4)),2));
    Double_t peakarea2=ffinal->GetParameter(5)*ffinal->GetParameter(7)*sqrt(2*TMath::Pi());
    Double_t peakareaerr2=peakarea2*sqrt(pow((ffinal->GetParError(5)/ffinal->GetParameter(5)),2.) + pow((ffinal->GetParError(7)/ffinal->GetParameter(7)),2));

    cout<<"peak at "<<ffinal->GetParameter(3)<<" area= "<<peakarea1<<" error= "<<peakareaerr1<<" rel. err="<<peakareaerr1/peakarea1<<endl;
    cout<<"peak at "<<ffinal->GetParameter(6)<<" area= "<<peakarea2<<" error= "<<peakareaerr2<<" rel. err="<<peakareaerr2/peakarea2<<endl;

    cout<<ffinal->GetParameter(3)<<"\t"<<peakarea1<<"\t"<<peakareaerr1<<"\t"<<peakareaerr1/peakarea1<<endl;
    cout<<ffinal->GetParameter(6)<<"\t"<<peakarea2<<"\t"<<peakareaerr2<<"\t"<<peakareaerr2/peakarea2<<endl;
}





