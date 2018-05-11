#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <unistd.h>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "TH1.h"
#include "TH2.h"

#include "TCutG.h"
#include "TLatex.h"
#include "TMath.h"
#include "TCanvas.h"

#include "TGraph.h"
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/GSLSimAnMinimizer.h"
#include "Math/Functor.h"


Double_t ioffset[128][128];
Double_t islope[128][128];
Double_t ioffseterr[128][128];
Double_t islopeerr[128][128];

Double_t iadcoffsetx[128];
Double_t iadcoffsety[128];

void ReadPulserCalibTable(char *infile,Int_t dssdno)
{
    Double_t dssd_cal[6][256][2];
    //clean up
    for (Int_t i=0;i<6;i++){
        for (Int_t j=0;j<256;j++){
            dssd_cal[i][j][0]=0.;
            dssd_cal[i][j][1]=1.;
        }
    }

    ifstream inpf(infile);
    if (inpf.fail()){
        cout<<"No Calibration table is given"<<endl;
        return;
    }

    cout<<"Start reading calibration table: "<<infile<<endl;
    Int_t dssd_index,strip_index;
    Double_t cal1,cal2;
    Int_t mm=0;

    while (inpf.good()){
    //for (Int_t i=0;i<100;i++){
        inpf>>dssd_index>>strip_index>>cal1>>cal2;
        dssd_cal[dssd_index][strip_index][0]=cal1;
        dssd_cal[dssd_index][strip_index][1]=cal2;
        mm++;
    }

    for (int i=0;i<256;i++){
        //if (i<127) iadcoffsetx[i]=-dssd_cal[dssdno][i][0]/dssd_cal[dssdno][i][1];
        //else iadcoffsety[i-128]=-dssd_cal[dssdno][i-128][0]/dssd_cal[dssdno][i-128][1];
        if (i<127) iadcoffsetx[i]=dssd_cal[dssdno][i][0];
        else iadcoffsety[i-128]=dssd_cal[dssdno][i-128][0];

        if (i<127) cout<<iadcoffsetx[i]<<endl;
        else cout<<iadcoffsety[i-128]<<endl;
    }
    cout<<"Read "<<mm<<" line"<<endl;
    inpf.close();
}


void correctemptypixel(char* infile,char* outfile)
{

    Int_t nearbypixelsR=1;
    Double_t badpixelserrorfactor=1;
    Double_t verybadpixelserrorfactor=1;

    TH2F* hpixels=new TH2F("hpixels","hpixels",128,0,128,128,0,128);
    std::ifstream ifs(infile);
    Int_t nlines=0;
    Double_t offset[16384];
    Double_t slope[16384];
    Double_t offseterr[16384];
    Double_t slopeerr[16384];
    Int_t np[16384];
    Double_t chi2ndf[16384];

    Bool_t chFlag[16384];

    Double_t aveoffset=0;
    Double_t aveslope=0;
    Double_t aveoffseterr=0;
    Double_t aveslopeerr=0;
    Int_t nave=0;

    for (int i=0;i<16384;i++){
        //assign defaut value
        chFlag[i]=false;
        offset[i]=-99999;
        slope[i]=1;
        offseterr[i]=1;
        slopeerr[i]=1;
        np[i]=1;
        chi2ndf[i]=1;
    }

    Int_t ich;
    while (!ifs.eof()){
        ifs>>ich;
        ifs>>offset[ich]>>slope[ich]>>offseterr[ich]>>slopeerr[ich]>>np[ich]>>chi2ndf[ich];
        //cout<<ich<<"\t"<<offset[ich]<<"\t"<<slope[ich]<<"\t"<<offseterr[ich]<<"\t"<<slopeerr[ich]<<"\t"<<temp<<"\t"<<temp<<endl;;
        if (offset[ich]!=-99999) chFlag[ich]=true;
        if (np[ich]>3){
            aveoffset=aveoffset+offset[ich];
            aveslope=aveslope+slope[ich];
            aveoffseterr=aveoffseterr+offseterr[ich];
            aveslopeerr=aveslopeerr+slopeerr[ich];
            nave++;
        }

        Int_t x=ich/128;
        Int_t y=ich%128;
        if (chFlag[ich]) hpixels->Fill(x,y);
        nlines++;
    }
    //!calculated average values for "good" pixels
    aveoffset=aveoffset/(Double_t)nave;
    aveoffseterr=aveoffseterr/(Double_t)nave;
    aveslope=aveslope/(Double_t)nave;
    aveslopeerr=aveslopeerr/(Double_t)nave;

    cout<<aveoffset<<"-"<<aveoffseterr<<"-"<<aveslope<<"-"<<aveslopeerr<<"-"<<endl;

    Int_t nbadpixels=0;
    Int_t nverybadpixels=0;

    for (int i=0;i<16384;i++){
        //!assign average value, or neighboring pixel value for "empty" pixels + 5 time more uncertainty
        if (chFlag[i]) continue;

        offset[i]=aveoffset;
        slope[i]=aveslope;
        offseterr[i]=aveoffseterr+aveoffseterr*verybadpixelserrorfactor;
        slopeerr[i]=aveslopeerr+aveslopeerr*verybadpixelserrorfactor;
        Int_t x=i/128;
        Int_t y=i%128;
        Int_t maxnp=0;
        Int_t maxnpj=0;
        for (int j=0;j<16384;j++){
            if (!chFlag[j]) continue; //ommit bad nearby pixels
            Int_t xx=j/128;
            Int_t yy=j%128;
            Double_t r=sqrt((Double_t)((x-xx)*(x-xx)+(y-yy)*(y-yy)));
            if (r<nearbypixelsR+1){
                if (np[j]>maxnp){
                    maxnp=np[j];
                    maxnpj=j;
                }
            }
        }
        if (maxnp>0){
            if(y==0) cout<<x<<"-"<<y<<"-"<<maxnpj/128<<"-"<<maxnpj%128<<"-"<<maxnpj<<endl;
            offset[i]=offset[maxnpj];
            slope[i]=slope[maxnpj];
            offseterr[i]=offseterr[maxnpj]+offseterr[maxnpj]*badpixelserrorfactor;
            slopeerr[i]=slopeerr[maxnpj]+slopeerr[maxnpj]*badpixelserrorfactor;
            nbadpixels++;
        }else{
            nverybadpixels++;
        }
    }
    cout<<"Number of bad pixels = "<<nbadpixels<<" = "<<(Double_t)nbadpixels/16384*100<<" %"<<endl;
    cout<<"Number of very bad pixels = "<<nverybadpixels<<" = "<<(Double_t)nverybadpixels/16384*100<<" %"<<endl;

    std::ofstream ofs(outfile);
    for (int ich=0;ich<16384;ich++){
        ofs<<ich<<"\t"<<offset[ich]<<"\t"<<slope[ich]<<"\t"<<offseterr[ich]<<"\t"<<slopeerr[ich]<<"\t"<<np[ich]<<"\t"<<chi2ndf[ich]<<endl;
    }
    hpixels->Draw("colz");
}

/*
double RosenBrock(const double *xx )
{
  Double_t sumchi2=0.;
  for (int i=0;i<128;i++){
      if (i<127){
          for (int j=0;j<128;j++){
              sumchi2=sumchi2+(islope[i][j]-xx[256+i]/xx[j])*(islope[i][j]-xx[256+i]/xx[j])/islopeerr[i][j]/islopeerr[i][j]+
                      (ioffset[i][j]-((-xx[256+i]*iadcoffsetx[i]-xx[j]*iadcoffsety[j])/xx[j]))*(ioffset[i][j]-((-xx[256+i]*iadcoffsetx[i]-xx[j]*iadcoffsety[j])/xx[j]))/ioffseterr[i][j]/ioffseterr[i][j];
          }
      }else{
          for (int j=0;j<128;j++){
              sumchi2=sumchi2+(islope[i][j]-xx[510]/xx[j])*(islope[i][j]-xx[510]/xx[j])/islopeerr[i][j]/islopeerr[i][j]+
                      (ioffset[i][j]-((-xx[510]*iadcoffsetx[i]+xx[j]*iadcoffsety[j])/xx[j]))*(ioffset[i][j]-((xx[510]*iadcoffsetx[i]-xx[j]*iadcoffsety[j])/xx[j]))/ioffseterr[i][j]/ioffseterr[i][j];
          }
      }
  }
  return sumchi2;
}
*/

/*
double RosenBrock(const double *xx )
{
  Double_t sumchi2=0.;
  for (int i=0;i<128;i++){
      if (i<127){
          for (int j=0;j<128;j++){
              sumchi2=sumchi2+(islope[i][j]-xx[256+i]/xx[j])*(islope[i][j]-xx[256+i]/xx[j])/islopeerr[i][j]/islopeerr[i][j]+
                      (ioffset[i][j]-((-xx[256+i]*iadcoffsetx[i]+xx[j]*iadcoffsety[j])/xx[j]))*(ioffset[i][j]-((-xx[256+i]*iadcoffsetx[i]+xx[j]*iadcoffsety[j])/xx[j]))/ioffseterr[i][j]/ioffseterr[i][j];
          }
      }else{
          for (int j=0;j<128;j++){
              sumchi2=sumchi2+(islope[i][j]-1./xx[j])*(islope[i][j]-1./xx[j])/islopeerr[i][j]/islopeerr[i][j]+
                      (ioffset[i][j]-((1.+xx[j]*iadcoffsety[j])/xx[j]))*(ioffset[i][j]-((1.+xx[j]*iadcoffsety[j])/xx[j]))/ioffseterr[i][j]/ioffseterr[i][j];
          }
      }
  }
  return sumchi2;
}

*/


double RosenBrock(const double *xx )
{
  Double_t sumchi2=0.;
  for (int i=0;i<128;i++){
      //if (i<127){
          for (int j=0;j<128;j++){
              sumchi2=sumchi2+(islope[i][j]-xx[256+i]/xx[j])*(islope[i][j]-xx[256+i]/xx[j])/islopeerr[i][j]/islopeerr[i][j]+
                      (ioffset[i][j]-((xx[383+i]-xx[128+j])/xx[j]))*(ioffset[i][j]-((xx[383+i]-xx[128+j])/xx[j]))/ioffseterr[i][j]/ioffseterr[i][j];
          }
      //}else{
      //    for (int j=0;j<128;j++){
      //        sumchi2=sumchi2+(islope[i][j]-1./xx[j])*(islope[i][j]-1./xx[j])/islopeerr[i][j]/islopeerr[i][j]+
      //                (ioffset[i][j]-((1.-xx[128+j])/xx[j]))*(ioffset[i][j]-((1.-xx[128+j])/xx[j]))/ioffseterr[i][j]/ioffseterr[i][j];
      //    }
      //}
  }
  return sumchi2;
}



void dominfull(char* infile, char* outfile, int dssd)
{

    //ReadPulserCalibTable("cal_R28_29_may2017.txt",dssd);
    //ReadPulserCalibTable("cal_fake.txt",dssd);
    std::ifstream ifs(infile);
    Int_t nlines=0;
    Double_t offset[16384];
    Double_t slope[16384];
    Double_t offseterr[16384];
    Double_t slopeerr[16384];


    for (int i=0;i<16384;i++){
        offset[i]=-241.142;
        slope[i]=1.02413;
        offseterr[i]=612.515*3;
        slopeerr[i]=0.0557309*3;
        //assign average value
    }

    Double_t temp;
    Int_t ich;
    while (!ifs.eof()){
        ifs>>ich;
        ifs>>offset[ich]>>slope[ich]>>offseterr[ich]>>slopeerr[ich]>>temp>>temp;
        //cout<<ich<<"\t"<<offset[ich]<<"\t"<<slope[ich]<<"\t"<<offseterr[ich]<<"\t"<<slopeerr[ich]<<"\t"<<temp<<"\t"<<temp<<endl;;
        nlines++;
    }

    for (int i=0;i<128;i++){ //dimension: 128x128
        for (int j=0;j<128;j++){
            int idxx=i*128+j;
            ioffset[i][j]=offset[idxx];
            ioffseterr[i][j]=offseterr[idxx];
            islope[i][j]=slope[idxx];
            islopeerr[i][j]=slopeerr[idxx];
        }
    }
    ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
    //ROOT::Math::GSLSimAnMinimizer min;

    //min.SetMaxFunctionCalls(1000000);
    //min.SetMaxIterations(100000);
    min.SetMaxFunctionCalls(1000000000);
    min.SetMaxIterations(100000000);

    min.SetTolerance(0.001);
    ROOT::Math::Functor f(&RosenBrock,511);
    min.SetFunction(f);
      // Set the free variables to be minimized!
    for (int i=0;i<511;i++){
        if (i<128||(i>=256&&i<383)) min.SetVariable(i,Form("x%d",i),0.6103515625, 0.001);
        else min.SetVariable(i,Form("x%d",i),0., 0.001);//offset
    }
    min.SetFixedVariable(336,"x336",0.6103515625);//gain of X strip 80 as ref
    min.SetFixedVariable(463,"x463",0.);//offset of X strip 80 (from a pulser run) as ref

    cout<<"Be patient! We are busy feeding our birds..."<<endl;
    //min.SetPrintLevel(1);
    double xxx[511];
    for (int i=0;i<511;i++){
        if (i<128||(i>=256&&i<383)) xxx[i]=1.;
        else xxx[i]=0.;
    }
    cout<<"start chisquare ="<<RosenBrock(xxx)<<endl;
    min.Minimize();
    const double *xs = min.X();
    cout<<"minimized chisquare ="<<RosenBrock(xs)<<endl;
    double fslope[6][256];
    double foffset[6][256];
    //! set default value
    for (int i=0;i<6;i++){
        for (int j=0;j<256;j++){
            fslope[i][j]=1.;
            foffset[i][j]=0.;
        }
    }
    for (int i=0;i<511;i++){
        //cout<<i<<"\t"<<xs[i]<<endl;
        if (i>=256&&i<383) {//x gain
            fslope[dssd][i-256]=xs[i];
        }else if (i<128){//y gain
            fslope[dssd][i+128]=xs[i];
        }else if (i>=383){//x offset
            foffset[dssd][i-383]=xs[i];
        }else{//y offset
            foffset[dssd][i-128+128]=xs[i];
        }
    }
    std::ofstream ofs(outfile,std::ofstream::out | std::ofstream::app);
    for (int i=0;i<6;i++){
        if (i!=dssd) continue;
        for (int j=0;j<256;j++){
            //if (j<128) foffset[i][j]=-fslope[i][j]*iadcoffsetx[j];
            //else foffset[i][j]=-fslope[i][j]*iadcoffsety[j-128];
            cout<<"dssd = "<<i<<"\t"<<j<<"\t"<<foffset[i][j]<<"\t"<<fslope[i][j]<<endl;
            ofs<<i<<"\t"<<j<<"\t"<<foffset[i][j]<<"\t"<<fslope[i][j]<<endl;
        }
    }
}

