#include <math.h>
#include <stdlib.h>
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <map>
#include <iostream>
#include <fstream>



Double_t betaeff=50.; //percentage of an isotope
Double_t neutroneff=68*0.95; //neutron detection efficiency in percentage

Double_t deltaxy=2.5;
Double_t dxbetamean=0;Double_t dxbetasigma=deltaxy/4.;
Double_t dybetamean=0;Double_t dybetasigma=deltaxy/4.;
Double_t betaneutronmodtime=0.000021; //21 us moderationtime
Double_t beamneutronmodtime=0.000027; //27 us moderationtime

Double_t ximpmean=64;Double_t ximpsigma=32;
Double_t yimpmean=64;Double_t yimpsigma=32;
Double_t xbetabkgmean=64;Double_t xbetabkgsigma=32;
Double_t ybetabkgmean=64;Double_t ybetabkgsigma=32;

Double_t xmin=0;Double_t xmax=128;
Double_t ymin=0;Double_t ymax=128;
Double_t tsoffset=3600; //! there is a bug on this offset (don't now why -> to be fixed later)
Double_t neuwbeamperctg=40.;

Int_t isSpatialDistFromHist=0;
TH1F* hhdx=NULL;
TH1F* hhdy=NULL;


TH2F* hbetadxy;



typedef struct {
    double T; 	 // Calibrated time
    double Tcorr; //correlated time
    double x,y,z;// number of pixel for AIDA, or number of tube for BELEN
    int type;
    int type2;//neutron index
    int evt;
} datatype;

void copydata(datatype* datain,datatype &dataout){
    dataout.T=datain->T;
    dataout.Tcorr=datain->Tcorr;
    dataout.x=datain->x;
    dataout.y=datain->y;
    dataout.z=datain->z;
    dataout.type=datain->type;
    dataout.Tcorr=datain->Tcorr;
    dataout.evt=datain->evt;
}

bool gendecay(TRandom3 *rseed,Double_t ximp,Double_t yimp,Double_t tsin,Double_t halflife, Double_t p1n,Double_t p2n, datatype &sim,Int_t &bflag,Int_t &nflag){
    //! Reset value (just to stay in the safe side)
    bflag=-1;
    nflag=-1;
    sim.T=-1.;
    sim.Tcorr=-1.;
    sim.x=-1.;
    sim.y=-1.;
    sim.z=0.;
    sim.type=-1.;
    sim.type2=-1.;
    sim.evt=-1;

    //!gen decay start here!
    Double_t dxbeta=rseed->Gaus(dxbetamean,dxbetasigma);
    Double_t dybeta=rseed->Gaus(dybetamean,dybetasigma);

    if (isSpatialDistFromHist!=0){
        dxbeta=hhdx->GetRandom();
        dybeta=hhdy->GetRandom();
    }

    hbetadxy->Fill(dxbeta*0.56,dybeta*0.56);

    Double_t decayt=rseed->Exp(halflife/TMath::Log(2));
    sim.Tcorr=decayt;
    sim.T=tsin+decayt;
    Double_t xbeta=dxbeta+ximp;
    Double_t ybeta=dybeta+yimp;
    //! take into account boundary effect
    while(!((dxbeta*dxbeta<=deltaxy*deltaxy&&dybeta*dybeta<=deltaxy*deltaxy))||xbeta>=xmax||xbeta<xmin||ybeta>=ymax||ybeta<ymin){
        dxbeta=rseed->Gaus(dxbetamean,dxbetasigma);
        dybeta=rseed->Gaus(dybetamean,dybetasigma);
        xbeta=dxbeta+ximp;
        ybeta=dybeta+yimp;
    }
    sim.x=xbeta;
    sim.y=ybeta;
    sim.z=0;
    Double_t pbeta=rseed->Rndm()*100;
    Double_t ppn=rseed->Rndm()*100;

    if (pbeta<betaeff) bflag=1; else bflag=0;

    if (ppn>=(p1n+p2n)){//isobaric decay
        nflag=0;
    }else if (ppn>=p2n&&ppn<p1n+p2n){//decay with 1 delayed neutron
        nflag=1;
    }else{//decay with 2 delayed neutron
        nflag=2;
    }
    return true;
}

bool genneutron(TRandom3 *rseed,Double_t tsin,datatype &sim, Int_t &nflag){
    //! Reset value (just to stay in the safe side)
    nflag=-1;
    sim.T=-1.;
    sim.Tcorr=0;
    sim.x=-1.;
    sim.y=-1.;
    sim.z=0.;
    sim.type=-1.;
    sim.type2=-1.;
    sim.evt=-1;

    Double_t pneu=rseed->Rndm()*100;
    Double_t dtneu=rseed->Exp(betaneutronmodtime/TMath::Log(2));
    sim.T=tsin+dtneu;
    sim.Tcorr=dtneu;
    sim.x=0;
    sim.y=0;
    sim.z=0;
    if (pneu<=neutroneff) nflag=1; else nflag=0;
    return true;
}

bool genneutronwithbkg(Double_t prob,TRandom3 *rseed,Double_t tsin,datatype &sim, Int_t &nflag){
    //! Reset value (just to stay in the safe side)
    nflag=-1;
    sim.T=-1.;
    sim.Tcorr=0;
    sim.x=0.;
    sim.y=0.;
    sim.z=0.;
    sim.type=-1.;
    sim.type2=-1.;
    sim.evt=-1;

    Double_t pneu=rseed->Rndm()*100;
    Double_t dtneu=rseed->Exp(betaneutronmodtime/TMath::Log(2));
    sim.T=tsin+dtneu;
    sim.Tcorr=dtneu;

    if (pneu<=prob) nflag=1; else nflag=0;
    return true;
}
bool gen2neutronwithbkg(Double_t prob,TRandom3 *rseed,Double_t tsin,datatype &simneu1,datatype &simneu2, Int_t &nflag){
    //! Reset value (just to stay in the safe side)
    nflag=-1;
    simneu1.T=-1.;
    simneu1.Tcorr=0;
    simneu1.x=0.;
    simneu1.y=0.;
    simneu1.z=0.;
    simneu1.type=-1.;
    simneu1.type2=-1.;
    simneu1.evt=-1;

    simneu2.T=-1.;
    simneu2.Tcorr=0;
    simneu1.x=0.;
    simneu1.y=0.;
    simneu1.z=0.;
    simneu2.type=-1.;
    simneu2.type2=-1.;
    simneu2.evt=-1;

    Double_t pneu=rseed->Rndm()*100;


    //! 1st neutron
    Double_t dtneu=rseed->Exp(betaneutronmodtime/TMath::Log(2));
    simneu1.T=tsin+dtneu;
    simneu1.Tcorr=dtneu;

    //! 2nd neutron
    dtneu=rseed->Exp(betaneutronmodtime/TMath::Log(2));
    simneu2.T=tsin+dtneu;
    simneu2.Tcorr=dtneu;

    if (pneu<=prob) nflag=1; else nflag=0;
    return true;
}

void simabinitio()
{

    Int_t isFixImplantPosition = 0;
    Int_t isSpatialDistFromHist = 0;

    //! some default rate values
    Double_t BeamTime=3600;
    Double_t rate=141.;
    Double_t isoperctg=60.;
    Double_t betabkgrateg=1.;
    Double_t betabkgrateu=1.;
    Double_t neurndbkgrate=1.;//random single neutron background  per second
    Double_t r2neurndbkgrate=1.;//random 2 neutron background per second

    Double_t percentage_bkgnb=7.;
    Double_t percentage_bkg2nb=1.;


    //! reading beam condition files
    std::ifstream ifscond("simparmsex.txt");
    std::string line;

    std::string line_head;
    Double_t  line_val;

    while (!ifscond.eof()){
        std::getline(ifscond,line);
        std::stringstream ss(line);
        ss>>line_head;
        if (line_head.at(0)=='#') continue;
        ss>>line_val;
        //cout<<line_head<<"-"<<line_val<<endl;

        if (line_head=="isfiximppos") isFixImplantPosition=(Int_t)line_val;
        if (line_head=="ishistgausbkg") isSpatialDistFromHist=(Int_t)line_val;

        if (line_head=="beamtime") BeamTime=line_val;
        if (line_head=="beamrate") rate=line_val;
        if (line_head=="isoperctg") isoperctg=line_val;
        if (line_head=="betabkgrateg") betabkgrateg=line_val;
        if (line_head=="betabkgrateu") betabkgrateu=line_val;
        if (line_head=="neubkgrate") neurndbkgrate=line_val;
        if (line_head=="r2neubkgrate") r2neurndbkgrate=line_val;
        if (line_head=="randbetaneuperctg") percentage_bkgnb=line_val;
        if (line_head=="randbeta2neuperctg") percentage_bkg2nb=line_val;
        if (line_head=="betaeff") betaeff=line_val;
        if (line_head=="neueff") neutroneff=line_val;
        if (line_head=="betaneutronmodtime") betaneutronmodtime=line_val;
        if (line_head=="beamneutronmodtime") beamneutronmodtime=line_val;

        if (line_head=="xmin") xmin=line_val;
        if (line_head=="xmax") xmax=line_val;
        if (line_head=="ymin") ymin=line_val;
        if (line_head=="ymax") ymax=line_val;

        if (line_head=="ximpmean") ximpmean=line_val;
        if (line_head=="ximpsigma") ximpsigma=line_val;
        if (line_head=="yimpmean") yimpmean=line_val;
        if (line_head=="yimpsigma") yimpsigma=line_val;

        if (line_head=="deltaxylimit") deltaxy=line_val;
        if (line_head=="dxbetamean") dxbetamean=line_val;
        if (line_head=="dxbetasigma") dxbetasigma=line_val;
        if (line_head=="dybetamean") dxbetamean=line_val;
        if (line_head=="dybetasigma") dxbetasigma=line_val;


        if (line_head=="xbetabkgmean") xbetabkgmean=line_val;
        if (line_head=="xbetabkgsigma") xbetabkgsigma=line_val;
        if (line_head=="ybetabkgmean") ybetabkgmean=line_val;
        if (line_head=="ybetabkgsigma") ybetabkgsigma=line_val;


        if (line_head=="tsoffset") tsoffset=line_val;
        if (line_head=="neuwbeamperctg") neuwbeamperctg=line_val;
    }

    cout<<"*****************\nSimulation parameters:\n"<<endl;
    cout<<"Fix implant position?  "<<isFixImplantPosition<<endl;
    cout<<"Beam time = "<<BeamTime<<endl;
    cout<<"Beam rate = "<<rate<<endl;
    cout<<"Isotope percentage = "<<isoperctg<<endl;
    cout<<"Beta Gaussian background = "<<betabkgrateg<<endl;
    cout<<"Beta Uniform background = "<<betabkgrateu<<endl;
    cout<<"Neutron background = "<<neurndbkgrate<<endl;
    cout<<"2 Neutron background = "<<r2neurndbkgrate<<endl;
    cout<<"Percentage of neutron correlated with random beta background = "<<percentage_bkgnb<<endl;
    cout<<"Percentage of 2 neutron correlated with random beta background = "<<percentage_bkg2nb<<endl;
    cout<<"Beta Efficiency = "<<betaeff<<endl;
    cout<<"Neutron Efficiency = "<<neutroneff<<endl;
    cout<<"Beta neutron moderation time = "<<betaneutronmodtime<<endl;
    cout<<"Beam neutron moderation time = "<<beamneutronmodtime<<endl;

    cout<<"xmin = "<<xmin<<endl;
    cout<<"xmax = "<<xmax<<endl;
    cout<<"ymin = "<<ymin<<endl;
    cout<<"ymax = "<<ymax<<endl;

    cout<<"ximpmean = "<<ximpmean<<endl;
    cout<<"ximpsigma = "<<ximpsigma<<endl;
    cout<<"yimpmean = "<<yimpmean<<endl;
    cout<<"yimpsigma = "<<yimpsigma<<endl;

    cout<<"deltaxylimit = "<<deltaxy<<endl;
    cout<<"dxbetamean = "<<dxbetamean<<endl;
    cout<<"dxbetasigma = "<<dxbetasigma<<endl;
    cout<<"dybetamean = "<<dybetamean<<endl;
    cout<<"dybetasigma = "<<dybetasigma<<endl;


    cout<<"xbetabkgmean = "<<ximpmean<<endl;
    cout<<"xbetabkgsigma = "<<ximpsigma<<endl;
    cout<<"ybetabkgmean = "<<ybetabkgmean<<endl;
    cout<<"ybetabkgsigma = "<<ybetabkgsigma<<endl;

    cout<<"ybetabkgmean = "<<ybetabkgmean<<endl;
    cout<<"ybetabkgsigma = "<<ybetabkgsigma<<endl;

    cout<<"Offset time stamp = "<<tsoffset<<endl;
    cout<<"Percentage of neutron correlated with beam = "<<neuwbeamperctg<<endl;

    cout<<"*****************\n"<<endl;


    TH1F* hhx=NULL;
    TH1F* hhy=NULL;
    TH1F* hhimpx=NULL;
    TH1F* hhimpy=NULL;


    if (isSpatialDistFromHist!=0){
        cout<<"reading histograms hy.root, hx.root, himpx.root, himpy.root, hdx.root and hdy.root for beta background distribution"<<endl;
        TFile* f1x = TFile::Open("hx.root");
        hhx=(TH1F*)f1x->Get("hx");
        TFile* f1y = TFile::Open("hy.root");
        hhy=(TH1F*)f1y->Get("hy");
        TFile* f1impx = TFile::Open("himpx.root");
        hhimpx=(TH1F*)f1impx->Get("himpx");
        TFile* f1impy = TFile::Open("himpy.root");
        hhimpy=(TH1F*)f1impy->Get("himpy");
        TFile* f1hdx = TFile::Open("hdx.root");
        hhdx=(TH1F*)f1hdx->Get("hdx");
        TFile* f1hdy = TFile::Open("hdy.root");
        hhdy=(TH1F*)f1hdy->Get("hdy");
    }
    //return;

    //! declare default parameters
    Double_t halflife=0.128;//half life parent in second
    Double_t halflife2=0.468;//half life daugter in second
    Double_t halflife3=0.637;//half life delayed 1 neutron in second
    Double_t halflife4=2.8;//half life granddaugter second
    Double_t halflife5=5.7;//half life granddaugter in delayed 1 neutron branch in second
    Double_t halflife6=1.222;//half life daugter in delayed 2 neutron bacnh in second
    Double_t halflife7=13.2;//half life grand-granddaugter in second
    Double_t halflife8=32.6;//half life grand-granddaugter in delayed 1 neutron bacnh in second
    Double_t halflife9=10.2;//half life granddaugter in delayed 2 neutron bacnh in second

    Double_t p1n=30.;//pn value in %
    Double_t p2n=0.;//pn value in %
    Double_t p1n2=30.3;//pn value of daugter in %
    Double_t p1n3=7.2;//pn value of daugter in 1 delayed neutron branch in %
    Double_t p1n4=0.00001;//pn value of grand daugter in %
    Double_t p1n5=0.00001;//pn value of grand-daugter in 1 delayed neutron branch %
    Double_t p1n6=0.0;//pn value of grand-daugter in 1 delayed neutron branch %
    Double_t p1n9=0.0;//pn value of grand-daugter in 1 delayed neutron branch %


    //! read from file (copy from fitting script)
    std::ifstream ifs("parmsex.txt");
    Int_t rino;
    Double_t temp;
    Int_t knri=9;
    Bool_t flagfix[knri][3];
    Double_t decayparms[knri][3];
    Double_t decayparms_p[knri][3];
    Double_t decayparms_m[knri][3];
    for (int i=0;i<knri;i++){
        ifs>>rino;
        for (int j=0;j<3;j++){
            ifs>>temp;
            if (temp>=0){
                flagfix[i][j]=true;
                if (j==0) {
                    decayparms[i][j]=temp;
                    ifs>>temp;decayparms_m[i][j]=temp;decayparms_p[i][j]=temp;
                }else{
                    decayparms[i][j]=temp;
                    ifs>>temp;decayparms_p[i][j]=temp;decayparms_m[i][j]=temp;
                }
            }else{
                flagfix[i][j]=false;
                if (j==0) {
                    decayparms[i][j]=(-temp);
                    ifs>>temp;decayparms_p[i][j]=temp;
                    decayparms_m[i][j]=temp;
                }else{
                    decayparms[i][j]=-temp;
                    ifs>>temp;decayparms_p[i][j]=temp+decayparms[i][j];
                    decayparms_m[i][j]=decayparms[i][j]-temp;
                }

            }
        }
    }

    halflife=decayparms[0][0];
    halflife2=decayparms[1][0];
    halflife3=decayparms[2][0];
    halflife4=decayparms[3][0];
    halflife5=decayparms[4][0];
    halflife6=decayparms[5][0];
    halflife7=decayparms[6][0];
    halflife8=decayparms[7][0];
    halflife9=decayparms[8][0];
    p1n=decayparms[0][1]*100;//pn value in %
    p2n=decayparms[0][2]*100;;//pn value in %
    p1n2=decayparms[1][1]*100;//pn value of daugter in %
    p1n3=decayparms[2][1]*100;//pn value of daugter in 1 delayed neutron branch in %

    p1n4=decayparms[3][1]*100;//pn value of grand daugter in %
    p1n5=decayparms[4][1]*100;//pn value of grand-daugter in 1 delayed neutron branch %
    p1n6=decayparms[5][1]*100;//pn value of grand-daugter in 1 delayed neutron branch %
    p1n9=decayparms[8][1]*100;//pn value of grand-daugter in 1 delayed neutron branch %

    cout<<"P2n="<<p2n<<endl;
   Bool_t enableFlag[17]={true,//1
                         true,//2
                         true,//3
                         true,//4
                         true,//5
                         true,//6
                         true,//7
                         true,//8
                         true,//9
                         true,//10
                         true,//11
                         true,//12
                         true,//13
                         true,//14
                         true,//15
                         true,//16
                         true};//17

   /*
   for (Int_t i=0;i<17;i++) enableFlag[i]=false;
  enableFlag[0]=true;
  enableFlag[3]=true;
  enableFlag[5]=true;
*/

  Bool_t neutronEnableFlag[24];for (int i=0;i<24;i++) neutronEnableFlag[i]=false;
  neutronEnableFlag[4]=true;
  neutronEnableFlag[6]=true;

  neutronEnableFlag[3]=true;
  neutronEnableFlag[8]=true;neutronEnableFlag[10]=true;
  neutronEnableFlag[14]=true;
  neutronEnableFlag[16]=true;neutronEnableFlag[17]=true;
  neutronEnableFlag[0]=false;//beam induced neutron


  neutronEnableFlag[20]=true;//constant neutron background
  neutronEnableFlag[21]=true;//correlated neutron background with beta background
  neutronEnableFlag[22]=true;//corrlated 2neutron background with beta background
  neutronEnableFlag[23]=true;//corrlated 2neutron background within unknown source



  std::multimap < double, datatype > betaMap;
  std::multimap < double, datatype >::iterator betaMap_it;
  std::multimap < double, datatype > neuMap;
  std::multimap < double, datatype >::iterator neuMap_it;


  betaMap.clear();
  neuMap.clear();
  TRandom3 *r = new TRandom3(0);

  TFile* ofile = new TFile("outtree.root","recreate");
  ofile->cd();
  datatype ion;
  ion.T=0;
  ion.Tcorr=0;
  ion.x=0;
  ion.y=0;
  ion.z=0;
  ion.type=-1;
  ion.type2=-1;
  ion.evt=-1;
  datatype beta;
  beta.T=0;
  beta.Tcorr=0;
  beta.x=0;
  beta.y=0;
  beta.z=0;
  beta.type=-1;
  beta.type2=-1;
  beta.evt=-1;
  datatype neu;
  neu.T=0;
  neu.Tcorr=0;
  neu.x=0;
  neu.y=0;
  neu.z=0;
  neu.type=-1;
  neu.type2=-1;
  neu.evt=-1;

  datatype veto;
  veto.T=0;
  veto.Tcorr=0;
  veto.x=0;
  veto.y=0;
  veto.z=0;
  veto.type=-1;
  veto.type2=-1;
  veto.evt=-1;


  TTree* treeion=new TTree("ion","aida tree ion");
  TTree* treebeta=new TTree("beta","aida tree beta");
  TTree* treeneu=new TTree("neutron","neutron");
  treeion->Branch("ion",&ion,"T/D:Tcorr/D:x/D:y/D:z/D:type/I:type2/I:evt/I");
  treebeta->Branch("beta",&beta,"T/D:Tcorr/D:x/D:y/D:z/D:type/I:type2/I:evt/I");
  treeneu->Branch("neutron",&neu,"T/D:Tcorr/D:x/D:y/D:z/D:type/I:type2/I:evt/I");
  treeneu->Branch("veto",&veto,"T/D:Tcorr/D:x/D:y/D:z/D:type/I:type2/I:evt/I");


  TH1F* hdtimp=new TH1F("hdtimp","hdtimp",2000,0,0.1);
  TH2F* hxyimp=new TH2F("hxyimp","hxyimp",148,-10*0.56,138*0.56,148,-10*0.56,138*0.56);
  TH1F* hpimp=new TH1F("hpimp","hpimp",100,0,100);

  hbetadxy=new TH2F("hbetadxy","hbetadxy",200,-6*0.56,6*0.56,200,-6*0.56,6*0.56);
  TH1F* hpbeta=new TH1F("hpbeta","hpbeta",100,0,100);
  TH2F* hxybeta=new TH2F("hxybeta","hxybeta",148,-10,138,148,-10,138);
  //TH1F* hdtbetaimp=new TH1F("hdtbetaimp","hdtbetaimp",2000,0,halflife*20);
  TH1F* hdtbetaimp=new TH1F("hdtbetaimp","hdtbetaimp",2000,-11,11);

  TH2F* hxybetabkgg=new TH2F("hxybetabkgg","hxybetabkgg",148,-10,138,148,-10,138);
  TH2F* hxybetabkgu=new TH2F("hxybetabkgu","hxybetabkgu",148,-10,138,148,-10,138);

  TH1F* hdbemneu=new TH1F("hdbemneu","hdbemneu",2000,0,beamneutronmodtime*20);
  TH1F* hdbetaneu3=new TH1F("hdbetaneu3","hdbetaneu3",2000,0,betaneutronmodtime*20);
  TH1F* hdbetaneu4=new TH1F("hdbetaneu4","hdbetaneu4",2000,0,betaneutronmodtime*20);
  TH1F* hdbetaneu61=new TH1F("hdbetaneu61","hdbetaneu61",2000,0,betaneutronmodtime*20);
  TH1F* hdbetaneu62=new TH1F("hdbetaneu62","hdbetaneu62",2000,0,betaneutronmodtime*20);

  TH1F* hdbetaimpwneu4=new TH1F("hdbetaimpwneu4","hdbetaimpwneu4",2000,-10,10);

  TH1F* hdecayall=new TH1F("hdecayall","hdecayall",2000,-10e9,10e9);
  TH1F* hdecayw1neu=new TH1F("hdecayw1neu","hdecayw1neu",2000,-10e9,10e9);
  TH1F* hdecayw1neu4=new TH1F("hdecayw1neu4","hdecayw1neu4",2000,-10e9,10e9);

  TH1F* hdecayw2neu=new TH1F("hdecayw2neu","hdecayw2neu",2000,-10e9,10e9);
  TH1F* hdecayw1neu2neu=new TH1F("hdecayw1neu2neu","hdecayw1neu2neu",2000,-10e9,10e9);


  Int_t nimplant=0;
  Double_t tsimp=tsoffset;
  Double_t tsbeta=0;
  Double_t tsneu=0;


  while (tsimp<BeamTime+tsoffset){
      //! beam outside AIDA
      Double_t ximp=r->Gaus(ximpmean,ximpsigma);
      Double_t yimp=r->Gaus(yimpmean,yimpsigma);

      if (isFixImplantPosition!=0){
          ximp=ximpmean;
          yimp=yimpmean;
      }

      if (isSpatialDistFromHist!=0) {
          ximp=hhimpx->GetRandom();
          yimp=hhimpy->GetRandom();
      }

      if (ximp<xmax&&yimp<ymax&&ximp>=xmin&&yimp>=ymin){
          Double_t dtimp=r->Exp(1/rate);
          tsimp=dtimp+tsimp;//time stamp implantation

          //! neutron associated with beam
          Double_t p1neui=r->Rndm()*100;
          Double_t dtp1neui=r->Exp(beamneutronmodtime/TMath::Log(2));
          tsneu=dtp1neui+tsimp;
          if (p1neui<neuwbeamperctg){//note this neuwbeamperctg includes the efficiency
              Double_t tsneunsi=tsneu;
              datatype aida;
              aida.Tcorr=dtp1neui;
              aida.x=0;
              aida.y=0;
              aida.z=0;
              aida.T=tsneunsi;
              aida.type=0;
              aida.type2=-1;
              aida.evt=nimplant;
              neuMap.insert(make_pair(tsneunsi,aida));
          }

          //! beam inside aida
          Double_t pimp=r->Rndm()*100;
          if (pimp<isoperctg){
              hdbemneu->Fill(dtp1neui);
              //!isotope after PID gate
              hxyimp->Fill(ximp*0.56,yimp*0.56);
              hdtimp->Fill(dtimp);
              hpimp->Fill(pimp);
              ion.T=tsimp;
              ion.x=ximp;
              ion.y=yimp;
              ion.type=0;
              ion.evt=nimplant;
              treeion->Fill();

              //! Parent decay
              datatype ri1;Int_t bflag1;Int_t nflag1;
              gendecay(r,ximp,yimp,ion.T,halflife,p1n,p2n,ri1,bflag1,nflag1);
              if (nflag1==0){//isobaric decay
                  //! register beta
                  if (bflag1>0){
                      datatype betadecay;
                      copydata(&ri1,betadecay);
                      betadecay.type=1;
                      betadecay.evt=nimplant;
                      betaMap.insert(make_pair(betadecay.T,betadecay));
                  }
                  if (bflag1>0) hdecayall->Fill(ri1.T*1e9-ion.T*1e9);

                  //! Daugter decay
                  datatype ri2;Int_t bflag2;Int_t nflag2;
                  gendecay(r,ximp,yimp,ri1.T,halflife2,p1n2,0.,ri2,bflag2,nflag2);
                  if (nflag2==0){//isobaric decay
                      //!register beta
                      if (bflag2>0){
                          datatype betadecay;
                          copydata(&ri2,betadecay);
                          betadecay.type=2;
                          betadecay.evt=nimplant;
                          betaMap.insert(make_pair(betadecay.T,betadecay));
                      }
                      if (bflag2>0) hdecayall->Fill(ri2.T*1e9-ion.T*1e9);

                      //! Granddaugter decay
                      datatype ri4;Int_t bflag4;Int_t nflag4;
                      gendecay(r,ximp,yimp,ri2.T,halflife4,p1n4,0.,ri4,bflag4,nflag4);
                      if (nflag4==0){//isobaric decay
                          //!register beta
                          if (bflag4>0){
                              datatype betadecay;
                              copydata(&ri4,betadecay);
                              betadecay.type=7;
                              betadecay.evt=nimplant;
                              betaMap.insert(make_pair(betadecay.T,betadecay));
                          }
                          if (bflag4>0) hdecayall->Fill(ri4.T*1e9-ion.T*1e9);

                          //! Grandgranddaugter decay (assume no n-emission)
                          datatype ri7;Int_t bflag7;Int_t nflag7;
                          gendecay(r,ximp,yimp,ri4.T,halflife7,0.,0.,ri7,bflag7,nflag7);
                          //!register beta
                          if (bflag7>0){
                              datatype betadecay;
                              copydata(&ri7,betadecay);
                              betadecay.type=12;
                              betadecay.evt=nimplant;
                              betaMap.insert(make_pair(betadecay.T,betadecay));
                          }
                          if (bflag7>0) hdecayall->Fill(ri7.T*1e9-ion.T*1e9);

                      }else if (nflag4==1){//delayed 1 neutron
                          //!register beta
                          if (bflag4>0){
                              datatype betadecay;
                              copydata(&ri4,betadecay);
                              betadecay.type=8;
                              betadecay.evt=nimplant;
                              betaMap.insert(make_pair(betadecay.T,betadecay));
                          }
                          if (bflag4>0) hdecayall->Fill(ri4.T*1e9-ion.T*1e9);

                          //!register neutron
                          datatype bd1n4;Int_t bd1nflag4;
                          genneutron(r,ri4.T,bd1n4,bd1nflag4);
                          if (bd1nflag4>0){
                              datatype neu;
                              copydata(&bd1n4,neu);
                              neu.type=8;
                              neu.type2=-1;
                              neu.evt=nimplant;
                              neuMap.insert(make_pair(neu.T,neu));
                          }
                          if (bflag4>0&&bd1nflag4>0) hdecayw1neu->Fill(ri4.T*1e9-ion.T*1e9);

                          //! decay in 1n branch (assume no n-emission)
                          datatype ri8;Int_t bflag8;Int_t nflag8;
                          gendecay(r,ximp,yimp,ri4.T,halflife8,0.,0.,ri8,bflag8,nflag8);
                          //!register beta
                          if (bflag8>0){
                              datatype betadecay;
                              copydata(&ri8,betadecay);
                              betadecay.type=13;
                              betadecay.evt=nimplant;
                              betaMap.insert(make_pair(betadecay.T,betadecay));
                          }
                          if (bflag8>0) hdecayall->Fill(ri8.T*1e9-ion.T*1e9);
                      }else{ cout<<"sth wrong"<<endl;}// I will leave it here and don't say anything

                  }else if (nflag2==1){//delayed 1 neutron
                      //!register beta
                      if (bflag2>0){
                          datatype betadecay;
                          copydata(&ri2,betadecay);
                          betadecay.type=3;
                          betadecay.evt=nimplant;
                          betaMap.insert(make_pair(betadecay.T,betadecay));
                      }
                      if (bflag2>0) hdecayall->Fill(ri2.T*1e9-ion.T*1e9);

                      //!register neutron
                      datatype bd1n2;Int_t bd1nflag2;
                      genneutron(r,ri2.T,bd1n2,bd1nflag2);
                      if (bd1nflag2>0){
                          datatype neu;
                          copydata(&bd1n2,neu);
                          neu.type=3;
                          neu.type2=-1;
                          neu.evt=nimplant;
                          neuMap.insert(make_pair(neu.T,neu));
                      }
                      if (bflag2>0&&bd1nflag2>0) hdecayw1neu->Fill(ri2.T*1e9-ion.T*1e9);

                      //! decay in 1 neutron branch
                      datatype ri5;Int_t bflag5;Int_t nflag5;
                      gendecay(r,ximp,yimp,ri2.T,halflife5,p1n5,0.,ri5,bflag5,nflag5);
                      if (nflag5==0){//isobaric decay
                          //!register beta
                          if (bflag5>0){
                              datatype betadecay;
                              copydata(&ri5,betadecay);
                              betadecay.type=9;
                              betadecay.evt=nimplant;
                              betaMap.insert(make_pair(betadecay.T,betadecay));
                          }
                          if (bflag5>0) hdecayall->Fill(ri5.T*1e9-ion.T*1e9);

                          //! continue decay from 5 (no neutron emission)
                          datatype ri8;Int_t bflag8;Int_t nflag8;
                          gendecay(r,ximp,yimp,ri5.T,halflife8,0.,0.,ri8,bflag8,nflag8);
                          //!register beta
                          if (bflag8>0){
                              datatype betadecay;
                              copydata(&ri8,betadecay);
                              betadecay.type=13;
                              betadecay.evt=nimplant;
                              betaMap.insert(make_pair(betadecay.T,betadecay));
                          }
                          if (bflag8>0) hdecayall->Fill(ri8.T*1e9-ion.T*1e9);

                      }else if (nflag5==1){//delayed 1 neutron
                          //!register beta
                          if (bflag5>0){
                              datatype betadecay;
                              copydata(&ri5,betadecay);
                              betadecay.type=14;
                              betadecay.evt=nimplant;
                              betaMap.insert(make_pair(betadecay.T,betadecay));
                          }
                          if (bflag5>0) hdecayall->Fill(ri5.T*1e9-ion.T*1e9);

                          //!register neutron
                          datatype bd1n5;Int_t bd1nflag5;
                          genneutron(r,ri5.T,bd1n5,bd1nflag5);
                          if (bd1nflag5>0){
                              datatype neu;
                              copydata(&bd1n5,neu);
                              neu.type=14;
                              neu.type2=-1;
                              neu.evt=nimplant;
                              neuMap.insert(make_pair(neu.T,neu));
                          }
                          if (bflag5>0&&bd1nflag5>0) hdecayw1neu->Fill(ri5.T*1e9-ion.T*1e9);
                          //!
                          //!nothing go up to here
                          //!
                      }else{ cout<<"sth wrong"<<endl;}// I will leave it here and don't say anything

                  }else{ cout<<"sth wrong"<<endl;}// I will leave it here and don't say anything

              }else if (nflag1==1){//delayed 1 neutron
                  Double_t betat;
                  //!register beta
                  if (bflag1>0){
                      datatype betadecay;
                      copydata(&ri1,betadecay);
                      betadecay.type=4;
                      betadecay.evt=nimplant;
                      betaMap.insert(make_pair(betadecay.T,betadecay));
                      betat=betadecay.T;
                  }
                  if (bflag1>0) hdecayall->Fill(ri1.T*1e9-ion.T*1e9);

                  //!register neutron
                  datatype bd1n1;Int_t bd1nflag1;
                  genneutron(r,ri1.T,bd1n1,bd1nflag1);
                  if (bd1nflag1>0){
                      datatype neu;
                      copydata(&bd1n1,neu);
                      neu.type=4;
                      neu.type2=-1;
                      neu.evt=nimplant;
                      neuMap.insert(make_pair(neu.T,neu));
                      hdbetaimpwneu4->Fill(ri1.T-ion.T);
                      if (bflag1>0) hdbetaneu4->Fill(neu.T-betat);
                  }
                  if (bflag1>0&&bd1nflag1>0) {
                      hdecayw1neu->Fill(ri1.T*1e9-ion.T*1e9);
                      hdecayw1neu4->Fill(ri1.T*1e9-ion.T*1e9);
                  }

                  //! daugter decay in 1 neutron emission branch
                  datatype ri3;Int_t bflag3;Int_t nflag3;
                  gendecay(r,ximp,yimp,ri1.T,halflife3,p1n3,0.,ri3,bflag3,nflag3);
                  if (nflag3==0){//isobaric decay
                      //!register beta
                      if (bflag3>0){
                          datatype betadecay;
                          copydata(&ri3,betadecay);
                          betadecay.type=5;
                          betadecay.evt=nimplant;
                          betaMap.insert(make_pair(betadecay.T,betadecay));
                      }
                      if (bflag3>0) hdecayall->Fill(ri3.T*1e9-ion.T*1e9);

                      //! decay at 5 from 3
                      datatype ri5;Int_t bflag5;Int_t nflag5;
                      gendecay(r,ximp,yimp,ri3.T,halflife5,p1n5,0.,ri5,bflag5,nflag5);
                      if (nflag5==0){//isobaric decay
                          //!register beta
                          if (bflag5>0){
                              datatype betadecay;
                              copydata(&ri5,betadecay);
                              betadecay.type=9;
                              betadecay.evt=nimplant;
                              betaMap.insert(make_pair(betadecay.T,betadecay));
                          }
                          if (bflag5>0) hdecayall->Fill(ri5.T*1e9-ion.T*1e9);

                          //! decay at 8 from 5 (assume no n-emission)
                          datatype ri8;Int_t bflag8;Int_t nflag8;
                          gendecay(r,ximp,yimp,ri5.T,halflife8,0.,0.,ri8,bflag8,nflag8);
                          //!register beta
                          if (bflag8>0){
                              datatype betadecay;
                              copydata(&ri8,betadecay);
                              betadecay.type=13;
                              betadecay.evt=nimplant;
                              betaMap.insert(make_pair(betadecay.T,betadecay));
                          }
                          if (bflag8>0) hdecayall->Fill(ri8.T*1e9-ion.T*1e9);
                      }else{//delayed 1 neutron
                          //!register beta
                          if (bflag5>0){
                              datatype betadecay;
                              copydata(&ri5,betadecay);
                              betadecay.type=14;
                              betadecay.evt=nimplant;
                              betaMap.insert(make_pair(betadecay.T,betadecay));
                          }
                          if (bflag5>0) hdecayall->Fill(ri5.T*1e9-ion.T*1e9);

                          //!register neutron
                          datatype bd1n5;Int_t bd1nflag5;
                          genneutron(r,ri5.T,bd1n5,bd1nflag5);
                          if (bd1nflag5>0){
                              datatype neu;
                              copydata(&bd1n5,neu);
                              neu.type=14;
                              neu.type2=-1;
                              neu.evt=nimplant;
                              neuMap.insert(make_pair(neu.T,neu));
                          }
                          if (bflag5>0&&bd1nflag5>0) hdecayw1neu->Fill(ri5.T*1e9-ion.T*1e9);
                          //!
                          //!nothing up to here
                          //!
                      }
                  }else if (nflag3==1){// delayed 1 neutron
                      //!register beta
                      if (bflag3>0){
                          datatype betadecay;
                          copydata(&ri3,betadecay);
                          betadecay.type=10;
                          betadecay.evt=nimplant;
                          betaMap.insert(make_pair(betadecay.T,betadecay));
                      }
                      if (bflag3>0) hdecayall->Fill(ri3.T*1e9-ion.T*1e9);

                      //!register neutron
                      datatype bd1n3;Int_t bd1nflag3;
                      genneutron(r,ri3.T,bd1n3,bd1nflag3);
                      if (bd1nflag3>0){
                          datatype neu;
                          copydata(&bd1n3,neu);
                          neu.type=10;
                          neu.type2=-1;
                          neu.evt=nimplant;
                          neuMap.insert(make_pair(neu.T,neu));
                      }
                      if (bflag3>0&&bd1nflag3>0) hdecayw1neu->Fill(ri3.T*1e9-ion.T*1e9);
                      //! go into 2n emission branch
                      datatype ri9;Int_t bflag9;Int_t nflag9;
                      gendecay(r,ximp,yimp,ri3.T,halflife9,p1n9,0.,ri9,bflag9,nflag9);
                      if (nflag9==0){//isobaric decay
                          //!register beta
                          if (bflag9>0){
                              datatype betadecay;
                              copydata(&ri9,betadecay);
                              betadecay.type=15;
                              betadecay.evt=nimplant;
                              betaMap.insert(make_pair(betadecay.T,betadecay));
                          }
                          if (bflag9>0) hdecayall->Fill(ri9.T*1e9-ion.T*1e9);

                      }else{
                          //!register beta
                          if (bflag9>0){
                              datatype betadecay;
                              copydata(&ri9,betadecay);
                              betadecay.type=17;
                              betadecay.evt=nimplant;
                              betaMap.insert(make_pair(betadecay.T,betadecay));
                          }
                          if (bflag9>0) hdecayall->Fill(ri9.T*1e9-ion.T*1e9);

                          //!register neutron
                          datatype bd1n9;Int_t bd1nflag9;
                          genneutron(r,ri9.T,bd1n9,bd1nflag9);
                          if (bd1nflag9>0){
                              datatype neu;
                              copydata(&bd1n9,neu);
                              neu.type=17;
                              neu.type2=-1;
                              neu.evt=nimplant;
                              neuMap.insert(make_pair(neu.T,neu));
                          }
                          if (bflag9>0&&bd1nflag9>0) hdecayw1neu->Fill(ri9.T*1e9-ion.T*1e9);
                      }

                  }else{ cout<<"sth wrong"<<endl;}// I will leave it here and don't say anything

              }else if (nflag1==2){//! delayed 2 neutron
                  //!register beta
                  if (bflag1>0){
                      datatype betadecay;
                      copydata(&ri1,betadecay);
                      betadecay.type=6;
                      betadecay.evt=nimplant;
                      betaMap.insert(make_pair(betadecay.T,betadecay));
                  }
                  if (bflag1>0) hdecayall->Fill(ri1.T*1e9-ion.T*1e9);

                  //!register 1st neutron
                  datatype bd2n1_1;Int_t bd2nflag1_1;
                  genneutron(r,ri1.T,bd2n1_1,bd2nflag1_1);
                  if (bd2nflag1_1>0){
                      datatype neu;
                      copydata(&bd2n1_1,neu);
                      neu.type=6;
                      neu.type2=-1;
                      neu.evt=nimplant;
                      neuMap.insert(make_pair(neu.T,neu));
                  }
                  //!register 2nd neutron
                  datatype bd2n1_2;Int_t bd2nflag1_2;
                  genneutron(r,ri1.T,bd2n1_2,bd2nflag1_2);
                  if (bd2nflag1_2>0){
                      datatype neu;
                      copydata(&bd2n1_2,neu);
                      neu.type=6;
                      neu.type2=-1;
                      neu.evt=nimplant;
                      neuMap.insert(make_pair(neu.T,neu));
                  }

                  if (bd2nflag1_1>0&&bd2nflag1_2>0&&bflag1>0) hdecayw2neu->Fill(ri1.T*1e9-ion.T*1e9);
                  if ((bd2nflag1_1>0&&bd2nflag1_2<=0&&bflag1>0)||(bd2nflag1_1<=0&&bd2nflag1_2>0&&bflag1>0)) hdecayw1neu2neu->Fill(ri1.T*1e9-ion.T*1e9);

                  //! continue decay
                  datatype ri6;Int_t bflag6;Int_t nflag6;
                  gendecay(r,ximp,yimp,ri1.T,halflife6,p1n6,0.,ri6,bflag6,nflag6);
                  if (nflag6==0){//isobaric decay
                      //!register beta
                      if (bflag6>0){
                          datatype betadecay;
                          copydata(&ri6,betadecay);
                          betadecay.type=11;
                          betadecay.evt=nimplant;
                          betaMap.insert(make_pair(betadecay.T,betadecay));
                      }
                      if (bflag6>0) hdecayall->Fill(ri6.T*1e9-ion.T*1e9);

                      //! continue continue decay
                      datatype ri9;Int_t bflag9;Int_t nflag9;
                      gendecay(r,ximp,yimp,ri6.T,halflife9,p1n9,0.,ri9,bflag9,nflag9);
                      if (nflag9==0){//isobaric decay
                          //!register beta
                          if (bflag9>0){
                              datatype betadecay;
                              copydata(&ri9,betadecay);
                              betadecay.type=15;
                              betadecay.evt=nimplant;
                              betaMap.insert(make_pair(betadecay.T,betadecay));
                          }
                          if (bflag9>0) hdecayall->Fill(ri9.T*1e9-ion.T*1e9);

                      }else{
                          //!register beta
                          if (bflag9>0){
                              datatype betadecay;
                              copydata(&ri9,betadecay);
                              betadecay.type=17;
                              betadecay.evt=nimplant;
                              betaMap.insert(make_pair(betadecay.T,betadecay));
                          }
                          if (bflag9>0) hdecayall->Fill(ri9.T*1e9-ion.T*1e9);

                          //!register neutron
                          datatype bd1n9;Int_t bd1nflag9;
                          genneutron(r,ri9.T,bd1n9,bd1nflag9);
                          if (bd1nflag9>0){
                              datatype neu;
                              copydata(&bd1n9,neu);
                              neu.type=17;
                              neu.type2=-1;
                              neu.evt=nimplant;
                              neuMap.insert(make_pair(neu.T,neu));
                          }
                          if (bflag9>0&&bd1nflag9>0) hdecayw1neu->Fill(ri9.T*1e9-ion.T*1e9);
                      }
                  }else if (nflag6==1){
                      //!register beta
                      if (bflag6>0){
                          datatype betadecay;
                          copydata(&ri6,betadecay);
                          betadecay.type=16;
                          betadecay.evt=nimplant;
                          betaMap.insert(make_pair(betadecay.T,betadecay));
                      }
                      if (bflag6>0) hdecayall->Fill(ri6.T*1e9-ion.T*1e9);

                      //! register neutron
                      datatype bd1n6;Int_t bd1nflag6;
                      genneutron(r,ri6.T,bd1n6,bd1nflag6);
                      if (bd1nflag6>0){
                          datatype neu;
                          copydata(&bd1n6,neu);
                          neu.type=16;
                          neu.type2=-1;
                          neu.evt=nimplant;
                          neuMap.insert(make_pair(neu.T,neu));
                      }
                      if (bflag6>0&&bd1nflag6>0) hdecayw1neu->Fill(ri6.T*1e9-ion.T*1e9);
                  }else{ cout<<"sth wrong"<<endl;}// I will leave it here and don't say anything
              }

              nimplant++;
          }//isotope ratio
      }//implant within dssd
  }

  cout<<"Finished step 1"<<endl;
  //! final time
  Double_t tend=0;
  if (tsimp>tsbeta) tend = tsimp; else tend=tsbeta;

  //! generate beta backgroud with gausian spatial distribution
  Double_t tsbetabkgg=tsoffset;
  Int_t iloopbkgg=0;
  Int_t nassociatedneu=0;
  Int_t nassociated2neu=0;
  while (tsbetabkgg<tend){
      Double_t xbetabkgg=r->Gaus(xbetabkgmean,xbetabkgsigma);
      Double_t ybetabkgg=r->Gaus(ybetabkgmean,ybetabkgsigma);

      if (isSpatialDistFromHist!=0) {
          xbetabkgg=hhx->GetRandom();
          ybetabkgg=hhy->GetRandom();
      }

      if (xbetabkgg<xmax&&ybetabkgg<ymax&&xbetabkgg>=xmin&&ybetabkgg>=ymin){
          Double_t dtbetabkgg=r->Exp(1/betabkgrateg);
          tsbetabkgg=dtbetabkgg+tsbetabkgg;
          iloopbkgg++;

          hxybetabkgg->Fill(xbetabkgg,ybetabkgg);

          datatype aida;
          aida.x=xbetabkgg;
          aida.y=ybetabkgg;
          aida.z=0;
          aida.T=tsbetabkgg;
          aida.Tcorr=-1;
          aida.type=20;
          aida.evt=-1;
          betaMap.insert(make_pair(tsbetabkgg,aida));

          //! single neutron associated with background beta
          Int_t bbkgd1nflag;
          datatype bd1nbkg;
          genneutronwithbkg(percentage_bkgnb,r,tsbetabkgg,bd1nbkg,bbkgd1nflag);
          if (bbkgd1nflag>0){
              datatype neu;
              copydata(&bd1nbkg,neu);
              neu.Tcorr=-1;
              neu.type=21;
              neu.type2=-1;
              neu.evt=-2;
              neuMap.insert(make_pair(neu.T,neu));
              nassociatedneu++;
          }

          //! 2 neutron associated with background beta
          Int_t bbkgd2nflag;
          datatype bd2nbkg1;
          datatype bd2nbkg2;
          gen2neutronwithbkg(percentage_bkg2nb,r,tsbetabkgg,bd2nbkg1,bd2nbkg2,bbkgd2nflag);
          if (bbkgd2nflag>0){
              //! register 1st neutron
              datatype neu1;
              copydata(&bd2nbkg1,neu1);
              neu1.Tcorr=-1;
              neu1.type=22;
              neu1.type2=-1;
              neu1.evt=-2;
              neuMap.insert(make_pair(neu1.T,neu1));
              //! register 2nd neutron
              datatype neu2;
              copydata(&bd2nbkg2,neu2);
              neu2.Tcorr=-1;
              neu2.type=22;
              neu2.type2=-1;
              neu2.evt=-2;
              neuMap.insert(make_pair(neu2.T,neu2));
              nassociated2neu++;
          }

      }
  }

  cout<<"bkg rate gaussian = "<<(Double_t)iloopbkgg/(tend-tsoffset)<<endl;
  cout<<"neutron rate asssociated with bkg rate gaussian = "<<(Double_t)nassociatedneu/(tend-tsoffset)<<endl;
  cout<<"2neutron rate asssociated with bkg rate gaussian = "<<(Double_t)nassociated2neu/(tend-tsoffset)<<endl;

  //! generate beta backgroud with uniform distribution
  Double_t tsbetabkgu=tsoffset;
  Int_t iloopbkgu=0;
  while (tsbetabkgu<tend){
      Double_t xbetabkgu=xmin+r->Rndm()*(xmax-xmin);
      Double_t ybetabkgu=ymin+r->Rndm()*(ymax-ymin);
      Double_t dtbetabkgu=r->Exp(1/betabkgrateu);
      tsbetabkgu=dtbetabkgu+tsbetabkgu;
      iloopbkgu++;
      hxybetabkgu->Fill(xbetabkgu,ybetabkgu);
      datatype aida;
      aida.x=xbetabkgu;
      aida.y=ybetabkgu;
      aida.z=0;
      aida.T=tsbetabkgu;
      aida.Tcorr=-1;
      aida.type=21;
      aida.evt=-1;
      betaMap.insert(make_pair(tsbetabkgu,aida));
  }
  cout<<"bkg rate uniform = "<<(Double_t)iloopbkgu/(tend-tsoffset)<<endl;

  //! generate random neutron backgroud
  Double_t tsneubkg=tsoffset;
  Int_t iloopneubkg=0;

  while (tsneubkg<tend){
      Double_t dtneubkg=r->Exp(1/neurndbkgrate);
      tsneubkg=dtneubkg+tsneubkg;
      iloopneubkg++;
      datatype aida;
      aida.x=0;
      aida.y=0;
      aida.z=0;
      aida.T=tsneubkg;
      aida.Tcorr=-1;
      aida.type=20;
      aida.type2=-1;
      aida.evt=-2;//associated with nothing
      neuMap.insert(make_pair(tsneubkg,aida));
  }
  cout<<"bkg rate random neutron = "<<(Double_t)iloopneubkg/(tend-tsoffset)<<endl;

  //! generate random 2 neutron backgroud with unknown source
  tsneubkg=tsoffset;
  iloopneubkg=0;

  while (tsneubkg<tend){
      //! timming of the unknown source
      Double_t dtneubkg=r->Exp(1/r2neurndbkgrate);
      tsneubkg=dtneubkg+tsneubkg;
      iloopneubkg++;
      //! correlated neutron

      //! 1st neutron
      Double_t dtneu=r->Exp(beamneutronmodtime/TMath::Log(2));//assume unknown source from upstream having different moderation time
      datatype simneu1;
      simneu1.T=tsneubkg+dtneu;
      simneu1.Tcorr=dtneu;
      simneu1.x=0;
      simneu1.y=0;
      simneu1.z=0;
      simneu1.type=23;
      simneu1.type2=-1;
      simneu1.evt=-2;//associated with nothing
      neuMap.insert(make_pair(simneu1.T,simneu1));

      //! 2nd neutron
      dtneu=r->Exp(beamneutronmodtime/TMath::Log(2));
      datatype simneu2;
      simneu2.T=tsneubkg+dtneu;
      simneu2.Tcorr=dtneu;
      simneu2.x=0;
      simneu2.y=0;
      simneu2.z=0;
      simneu2.type=23;
      simneu2.type2=-1;
      simneu2.evt=-2;//associated with nothing
      neuMap.insert(make_pair(simneu2.T,simneu2));
  }
  cout<<"bkg rate random 2 neutron = "<<(Double_t)iloopneubkg/(tend-tsoffset)<<endl;


  //! go through all entries in BETA map and fill tree
  for (betaMap_it=betaMap.begin();betaMap_it!=betaMap.end();betaMap_it++){
      datatype mbeta = betaMap_it->second;
      beta.x=mbeta.x;
      beta.y=mbeta.y;
      beta.z=mbeta.z;
      beta.T=mbeta.T;
      beta.Tcorr=mbeta.Tcorr;
      beta.type=mbeta.type;
      beta.type2=-1;
      beta.evt=mbeta.evt;
      if (beta.type==20||beta.type==21)
          treebeta->Fill();
      else{
          if (enableFlag[beta.type-1]) treebeta->Fill();
      }
  }
  cout<<"betatree filled "<<treebeta->GetEntries()<<" entries"<<endl;

  //! go through all entries in NEUTRON map and fill tree
  for (neuMap_it=neuMap.begin();neuMap_it!=neuMap.end();neuMap_it++){
      datatype mneu = neuMap_it->second;
      neu.x=mneu.x;
      neu.y=mneu.y;
      neu.z=mneu.z;
      neu.T=mneu.T;
      neu.Tcorr=mneu.Tcorr;
      neu.type=mneu.type;
      neu.type2=mneu.type2;
      neu.evt=mneu.evt;
      if (neutronEnableFlag[neu.type]) treeneu->Fill();
  }

  cout<<"\ntotal number of implant = "<<nimplant<<endl;
  cout<<"total number of beta = "<<treebeta->GetEntries()<<endl;
  cout<<"total number of neutron = "<<treeneu->GetEntries()<<endl;


  /*
  TCanvas* c1=new TCanvas("c1","c1",900,900);
  c1->Divide(3,3);
  c1->cd(1);
  hxyimp->Draw("colz");
  c1->cd(2);
  hdtimp->Draw();
  c1->cd(3);
  hpimp->Draw();
  c1->cd(4);
  hbetadxy->Draw("colz");
  c1->cd(5);
  hxybeta->Draw("colz");
  c1->cd(6);
  hpbeta->Draw();
  c1->cd(7);
  hdtbetaimp->Draw();
  c1->cd(8);
  hxybetabkgg->Draw("colz");
  c1->cd(9);
  hxybetabkgu->Draw("colz");


  TCanvas* c2=new TCanvas("c2","c2",900,700);
  c2->Divide(3,2);
  c2->cd(1);
  hdbemneu->Draw();
  c2->cd(2);
  hdbetaneu3->Draw();
  c2->cd(3);
  hdbetaneu4->Draw();
  c2->cd(4);
  hdbetaneu61->Draw();
  c2->cd(5);
  hdbetaneu62->Draw();


  c1->Write();
  c2->Write();
  */

  hxyimp->Write();
  hbetadxy->Write();
  hxybetabkgg->Write();
  hxybetabkgu->Write();


  hdecayall->Write();
  hdecayw1neu->Write();
  hdecayw1neu4->Write();
  hdecayw2neu->Write();
  hdecayw1neu2neu->Write();


  hdbetaimpwneu4->Write();
  hdtbetaimp->Write();
  treeion->Write();
  treebeta->Write();
  treeneu->Write();
  ofile->Close();
}


