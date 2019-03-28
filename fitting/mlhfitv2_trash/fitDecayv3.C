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



const Int_t knri=9;
const Int_t kmaxpar=5;
const Int_t kmaxndecay=10;
const Int_t kmaxpaths=100;
//Double_t neueff=0.66*(100-0.8)/100;
Double_t neueff=0.68;


//! variable
Int_t npaths=100;
Int_t ndecay[kmaxpaths];
Int_t decaymap[kmaxpaths][kmaxndecay];
Int_t nneu[kmaxpaths][kmaxndecay];

void getparms(Double_t* parms,Double_t* parmserr, Double_t* parmsmax, Double_t* parmsmin, Bool_t* isparmsfix,char* infile,Double_t nsig=5.)
{
    std::ifstream ifs(infile);
    Int_t rino;
    Double_t temp;
    Bool_t flagfix[knri][3];
    Double_t decayparms[knri][3];
    Double_t decayparms_err[knri][3];

    for (int i=0;i<knri;i++){
        ifs>>rino;
        for (int j=0;j<3;j++){
            ifs>>temp;
            if (temp>=0){
                flagfix[i][j]=true;
                decayparms[i][j]=temp;
            }else{
                flagfix[i][j]=false;
                decayparms[i][j]=(-temp);
            }
            ifs>>temp;decayparms_err[i][j]=temp;
        }
    }

    for (int i=0;i<knri;i++){
        for (int j=0;j<3;j++){
            if (j==0){//for half-life
                decayparms_err[i][j]=log(2)/decayparms[i][j]/decayparms[i][j]*decayparms_err[i][j];
                decayparms[i][j]=log(2)/decayparms[i][j];
            }
        }
    }


    //! calculate output
    for (int i=0;i<knri;i++){
        for (int j=0;j<2;j++){
            isparmsfix[j*knri+i]=flagfix[i][j];
            parms[j*knri+i]=decayparms[i][j];
            parmserr[j*knri+i]=decayparms_err[i][j];
            parmsmax[j*knri+i]=decayparms[i][j]+decayparms_err[i][j]*nsig;
            if ((decayparms[i][j]-decayparms_err[i][j]*nsig)>0)
                parmsmin[j*knri+i]=decayparms[i][j]-decayparms_err[i][j]*nsig;
            else
                parmsmin[j*knri+i]=log(2)/100000000000;
            //for pn
            if (j!=0) {parmsmin[j*knri+i]=0;parmsmax[j*knri+i]=1;}
        }
    }

    //! read p2n parrent
    parms[knri*2]=decayparms[0][2];
    parmserr[knri*2]=decayparms_err[0][2];
    parms[knri*2]=decayparms[0][2];
    parmsmin[knri*2]=0;parmsmax[knri*2]=1;
    isparmsfix[knri*2]=flagfix[0][2];

    //! read initial activity, background and random coincidence factor
    for (int i=0;i<5;i++){
        ifs>>parms[knri*2+i+1]>>parmsmin[knri*2+i+1]>>parmsmax[knri*2+i+1];
        parmserr[knri*2+i+1]=parms[knri*2+i+1]-parmsmin[knri*2+i+1];
        cout<<knri*2+i+1<<"\t"<<parms[knri*2+i+1]<<"\t"<<parmsmin[knri*2+i+1]<<"\t"<<parmsmax[knri*2+i+1]<<endl;
        isparmsfix[knri*2+i+1]=false;
        if (i==4) {
            parmsmin[knri*2+i+1]=0;parmsmax[knri*2+i+1]=1;
            if (parms[knri*2+i+1]>=0){
                isparmsfix[knri*2+i+1]=true;
            }else{
                parms[knri*2+i+1]=-parms[knri*2+i+1];
                isparmsfix[knri*2+i+1]=false;
            }
        }
    }
}


//! Global Bateaman function
Double_t corefcn(Int_t ndecay,Int_t* decaymap,Int_t* nneu, Double_t* b1n,Double_t* b2n,Double_t* lamda,Double_t N0,Double_t t){
    Double_t fcnret=0;

    Double_t factor1=1.;
    //only parrent decay p2n
    for (int i=0;i<ndecay-1;i++){
        if (nneu[i]==0){
            factor1=factor1 * (1-b1n[decaymap[i]-1]-b2n[decaymap[i]-1])*lamda[decaymap[i]-1];
        }else if (nneu[i]==1){
            factor1=factor1 * b1n[decaymap[i]-1]*lamda[decaymap[i]-1];
        }else{
            factor1=factor1 * b2n[decaymap[i]-1]*lamda[decaymap[i]-1];
        }
    }
    Double_t factor2=0;
    for (int i=0;i<ndecay;i++){
        Double_t factor2i=exp(-lamda[decaymap[i]-1]*t);
        Double_t factor2ij=1;
        for (int j=0;j<ndecay;j++)
            if (j!=i) factor2ij=factor2ij*(lamda[decaymap[j]-1]-lamda[decaymap[i]-1]);
        factor2=factor2+factor2i/factor2ij;
    }

    fcnret=factor1*N0*factor2;
    return fcnret;
}


//! Global function
Double_t fcn_gen(Double_t *x, Double_t *par) {
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];

    //! Parent nuclei
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0]);

    for (Int_t i=0;i<npaths;i++){
        returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0]);
    }
    return returnval;
}

//! Global function
Double_t fcn_gen_wneutron(Double_t *x, Double_t *par) {
    Double_t randcoinf=par[knri*2+3];
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];

    //! Parent nuclei
    //returnval+=lamda[0]*N0*exp(-lamda[0]*x[0]);

    for (Int_t i=0;i<npaths;i++){
        returnval+=neueff*pn[decaymap[i][ndecay[i]-1]-1]*lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(1-randcoinf);
    }

    //! decay with neutron part
    returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf);
    returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf);
    return returnval;
}

//! Global function
Double_t fcn_gen_w2neutron(Double_t *x, Double_t *par) {
    Double_t randcoinf=par[knri*2+3];
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];

    for (Int_t i=0;i<npaths;i++){
        returnval+=neueff*pn[decaymap[i][ndecay[i]-1]-1]*lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf;
    }
    //! decay with neutron part
    returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf;
    returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf;
    //! decay with 2 neutron part
    returnval+=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0]);
    return returnval;
}

// definition of shared parameter
// background function
int iparB[21] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};

// signal + background function
int iparSB[22] ={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,23};

// signal 2n + background function
int iparSB2[22] ={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,22, 23};

// Create the GlobalCHi2 structure


struct GlobalChi2 {
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2,
                ROOT::Math::IMultiGenFunction & f3) :
      fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}

   // parameter vector is first background (in common 1 and 2)
   // and then is signal (only in 2)
   double operator() (const double *par) const {
      double p1[21];
      for (int i = 0; i < 21; ++i) p1[i] = par[iparB[i] ];

      double p2[22];
      for (int i = 0; i < 22; ++i) p2[i] = par[iparSB[i] ];

      double p3[22];
      for (int i = 0; i < 22; ++i) p3[i] = par[iparSB2[i] ];

      return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3);
   }
   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
   const  ROOT::Math::IMultiGenFunction * fChi2_3;
};
Double_t parms[knri*2+6];
Double_t parmserr[knri*2+6];
Double_t parmsmax[knri*2+6];
Double_t parmsmin[knri*2+6];
Bool_t isparmsfix[knri*2+6];
Double_t mcparms[knri*2+6];

void mc(TRandom3* rseed)
{
    for (int i=0;i<knri*2;i++){
        if (isparmsfix[i]&&(parmserr[i]>0)){
            mcparms[i]=rseed->Gaus(parms[i],parmserr[i]);
        }
    }

    //! mc for neutron efficiency?
}

void fitDecay(char* fitname,char* infile,char* parmsfile, Int_t ninterations, char* outfile="out.root", Int_t binning=900,Int_t islastdssd=1){

    Double_t lowerlimit=0.05;
    Double_t upperlimit=9.05;
    Double_t nsigma=2.;

    TRandom3* rseed=new TRandom3;
    //! book file and tree
    TFile* fout=new TFile(outfile,"recreate");
    Double_t outparms[knri*2+6];
    Double_t outparmserr[knri*2+6];
    Int_t isvary[knri*2+6];

    for (int i=0;i<knri*2+6;i++){
        outparms[i]=0;
        outparmserr[i]=0;
        isvary[i]=0;
        if (isparmsfix[i]&&(parmserr[i]>0)){
            if(i<knri*2) isvary[i]=1;
        }
    }

    char tempchar1[1000];
    sprintf(tempchar1,"tree%s",fitname);
    TTree* treeout=new TTree(tempchar1,tempchar1);
    treeout->Branch("outparms",outparms,Form("outparms[%d]/D",knri*2+6));
    treeout->Branch("outparmserr",outparmserr,Form("outparmserr[%d]/D",knri*2+6));
    treeout->Branch("isvary",outparmserr,Form("isvary[%d]/I",knri*2+6));

    //! define decay map
    npaths=12;
    //! path 1(go for ri2)
    ndecay[0]=2;
    decaymap[0][0]=1;decaymap[0][1]=2;
    nneu[0][0]=0;
    //! path 2(go for ri4)
    ndecay[1]=3;
    decaymap[1][0]=1;decaymap[1][1]=2;decaymap[1][2]=4;
    nneu[1][0]=0;nneu[1][1]=0;
    //! path 3(go for ri7)
    ndecay[2]=4;
    decaymap[2][0]=1;decaymap[2][1]=2;decaymap[2][2]=4;decaymap[2][3]=7;
    nneu[2][0]=0;nneu[2][1]=0;nneu[2][2]=0;
    //! path4(go for ri3)
    ndecay[3]=2;
    decaymap[3][0]=1;decaymap[3][1]=3;
    nneu[3][0]=1;
    //!path5(go for ri6)
    ndecay[4]=2;
    decaymap[4][0]=1;decaymap[4][1]=6;
    nneu[4][0]=2;

    //! path6(go for ri5-route 1)
    ndecay[5]=3;
    decaymap[5][0]=1;decaymap[5][1]=2;decaymap[5][2]=5;
    nneu[5][0]=0;nneu[5][1]=1;
    //! path7(go for ri5-route 2)
    ndecay[6]=3;
    decaymap[6][0]=1;decaymap[6][1]=3;decaymap[6][2]=5;
    nneu[6][0]=1;nneu[6][1]=0;

    //! path8 (go for ri9-route 1)
    ndecay[7]=3;
    decaymap[7][0]=1;decaymap[7][1]=3;decaymap[7][2]=9;
    nneu[7][0]=1;nneu[7][1]=1;
    //! path9 (go for ri9-route 2)
    ndecay[8]=3;
    decaymap[8][0]=1;decaymap[8][1]=6;decaymap[8][2]=9;
    nneu[8][0]=2;nneu[8][1]=0;

    //! path10 (go for ri8-route 1)
    ndecay[9]=4;
    decaymap[9][0]=1;decaymap[9][1]=2;decaymap[9][2]=4;decaymap[9][3]=8;
    nneu[9][0]=0;nneu[9][1]=0;nneu[9][2]=1;

    //! path11 (go for ri8-route 2)
    ndecay[10]=4;
    decaymap[10][0]=1;decaymap[10][1]=2;decaymap[10][2]=5;decaymap[10][3]=8;
    nneu[10][0]=0;nneu[10][1]=1;nneu[10][2]=0;

    //! path12 (go for ri8-route 3)
    ndecay[11]=4;
    decaymap[11][0]=1;decaymap[11][1]=3;decaymap[11][2]=5;decaymap[11][3]=8;
    nneu[11][0]=1;nneu[11][1]=0;nneu[11][2]=0;


    getparms(parms,parmserr,parmsmax,parmsmin,isparmsfix,parmsfile,nsigma);

    //!******************************************Define All BETA decay function
    TF1* fB=new TF1("fB",fcn_gen,lowerlimit,upperlimit,21);
    fB->SetNpx(500);
    fB->SetLineWidth(2);
    fB->SetLineColor(8);
    //!******************************************Define Delayed Neutron decay function
    TF1* fSB=new TF1("fSB",fcn_gen_wneutron,lowerlimit,upperlimit,22);
    fSB->SetNpx(500);
    fSB->SetLineWidth(2);
    fSB->SetLineColor(8);
    //!******************************************Define Delayed 2 Neutron decay function
    TF1* fSB2=new TF1("fSB2",fcn_gen_w2neutron,lowerlimit,upperlimit,22);
    fSB2->SetNpx(500);
    fSB2->SetLineWidth(2);
    fSB2->SetLineColor(8);


    //!****************************************** All BETA decay function
    //! read input and get parameters
    for (int i=0;i<knri*2+2;i++){
        fB->SetParameter(i,parms[i]);
        fB->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (parmserr[i]==0||isparmsfix[i]){
           cout<<"fixed value p"<<i<<"\tval="<<parms[i]<<endl;
           fB->FixParameter(i,parms[i]);
        }else{
           cout<<"variable value p"<<i<<"\tval "<<parms[i]<<"\tmin="<<parmsmin[i]<<"\tmax="<<parmsmax[i]<<endl;
        }
    }
    fB->SetParameter(knri*2+2,parms[knri*2+2]);
    fB->SetParLimits(knri*2+2,parmsmin[knri*2+2],parmsmax[knri*2+2]);
    cout<<fB->Eval(1.)<<endl;

    //!******************************************Delayed Neutron decay function
    //! read input and get parameters
    for (int i=0;i<knri*2+4;i++){
        fSB->FixParameter(i,parms[i]);
        fSB->SetParameter(i,parms[i]);
        fSB->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (parmserr[i]==0||isparmsfix[i]){
           cout<<"sb fixed value p"<<i<<"\tval="<<parms[i]<<endl;
           fSB->FixParameter(i,parms[i]);
        }else{
           cout<<"sb variable value p"<<i<<"\tval "<<parms[i]<<"\tmin="<<parmsmin[i]<<"\tmax="<<parmsmax[i]<<endl;
        }
    }
    fSB->SetParameter(knri*2+2,parms[knri*2+3]);
    fSB->SetParLimits(knri*2+2,parmsmin[knri*2+3],parmsmax[knri*2+3]);
    fSB->SetParameter(knri*2+3,0.5);//random coincidence factor
    fSB->SetParLimits(knri*2+3,0.,1.);//random coincidence factor
    cout<<fSB->Eval(1.)<<endl;


    //!******************************************Delayed 2 Neutron decay function
    //! read input and get parameters
    for (int i=0;i<knri*2+2;i++){
        fSB2->SetParameter(i,parms[i]);
        fSB2->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (parmserr[i]==0||isparmsfix[i]){
           cout<<"fixed value p"<<i<<"\tval="<<parms[i]<<endl;
           fSB2->FixParameter(i,parms[i]);
        }else{
           cout<<"variable value p"<<i<<"\tval "<<parms[i]<<"\tmin="<<parmsmin[i]<<"\tmax="<<parmsmax[i]<<endl;
        }
    }
    fSB2->SetParameter(knri*2+2,parms[knri*2+4]);
    fSB2->SetParLimits(knri*2+2,parmsmin[knri*2+4],parmsmax[knri*2+4]);

    fSB2->SetParameter(knri*2+3,0.5);//random coincidence factor
    fSB2->SetParLimits(knri*2+3,0.,1.);//random coincidence factor

    cout<<fSB2->Eval(1.)<<endl;




    //!****************************GET HISTOGRAM FROM FILE********************
    //!
    TFile *f = TFile::Open(infile);
    sprintf(tempchar1,"tree");
    TTree* treeb=(TTree*) f->Get(tempchar1);
    if (islastdssd){
        treeb->Draw(Form("x>>hb(%d,%f,%f)",binning,lowerlimit,upperlimit),"","goff");
        treeb->Draw(Form("x>>hb1n(%d,%f,%f)",binning,lowerlimit,upperlimit),"y==1","goff");
        treeb->Draw(Form("x>>hb2n(%d,%f,%f)",binning,lowerlimit,upperlimit),"y==2","goff");
    }else{
        treeb->Draw(Form("x>>hb(%d,%f,%f)",binning,lowerlimit,upperlimit),"z<5","goff");
        treeb->Draw(Form("x>>hb1n(%d,%f,%f)",binning,lowerlimit,upperlimit),"y==1&&z<5","goff");
        treeb->Draw(Form("x>>hb2n(%d,%f,%f)",binning,lowerlimit,upperlimit),"y==2&&z<5","goff");
    }
    sprintf(tempchar1,"hb");
    TH1F* hdecay=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hb1n");
    TH1F* hdecay1n=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hb2n");
    TH1F* hdecay2n=(TH1F*) gDirectory->Get(tempchar1);

   TH1F * hB = (TH1F*) hdecay->Clone();
   TH1F * hSB = (TH1F*) hdecay1n->Clone();
   TH1F * hSB2 = (TH1F*) hdecay2n->Clone();

   /*
   if (rebin!=0) {
       hB->Rebin(rebin);
       hSB->Rebin(rebin);
       hSB2->Rebin(rebin);
   }
   */
   //!****************************GET HISTOGRAM FROM FILE********************


   ROOT::Math::WrappedMultiTF1 wfB(*fB,1);
   ROOT::Math::WrappedMultiTF1 wfSB(*fSB,1);
   ROOT::Math::WrappedMultiTF1 wfSB2(*fSB2,1);


   ROOT::Fit::DataOptions opt;
   ROOT::Fit::DataRange rangeB;
   // set the data range
   rangeB.SetRange(lowerlimit,upperlimit);
   ROOT::Fit::BinData dataB(opt,rangeB);
   ROOT::Fit::FillData(dataB, hB);

   ROOT::Fit::DataRange rangeSB;
   rangeSB.SetRange(lowerlimit,upperlimit);
   ROOT::Fit::BinData dataSB(opt,rangeSB);
   ROOT::Fit::FillData(dataSB, hSB);


   ROOT::Fit::DataRange rangeSB2;
   rangeSB2.SetRange(lowerlimit,upperlimit);
   ROOT::Fit::BinData dataSB2(opt,rangeSB2);
   ROOT::Fit::FillData(dataSB2, hSB2);

   ROOT::Fit::Chi2Function chi2_B(dataB, wfB);
   ROOT::Fit::Chi2Function chi2_SB(dataSB, wfSB);
   ROOT::Fit::Chi2Function chi2_SB2(dataSB2, wfSB2);
   GlobalChi2 globalChi2(chi2_B, chi2_SB, chi2_SB2);

   ROOT::Fit::Fitter fitter;

   // create before the parameter settings in order to fix or set range on them

   cout<<"\n\n*****SETTING PARAMETERS......\n"<<endl;
   fitter.Config().SetParamsSettings(knri*2+6,parms);
   for (int i=0;i<knri*2+6;i++){
       if (parmserr[i]==0||isparmsfix[i]){
          cout<<"fixed value p"<<i<<"\tval="<<parms[i]<<endl;
          fitter.Config().ParSettings(i).Fix();
       }else{
          fitter.Config().ParSettings(i).SetLimits(parmsmin[i],parmsmax[i]);
          cout<<"variable value p"<<i<<"\tval "<<parms[i]<<"\tmin="<<parmsmin[i]<<"\tmax="<<parmsmax[i]<<endl;
       }
   }

   fitter.Config().MinimizerOptions().SetPrintLevel(0);
   fitter.Config().SetMinimizer("Minuit2","Migrad");


   cout<<"\n\n************FITTING.........."<<endl;
   //! fit FCN function directly
   // (specify optionally data size and flag to indicate that is a chi2 fit)
   for (int i=0;i<knri*2+6;i++) mcparms[i]=parms[i];
   for(int i=0;i<ninterations;i++) {
       cout<<"mc "<<i+1<<endl;
       mc(rseed);
       fitter.Config().SetParamsSettings(knri*2+6,mcparms);
       fitter.FitFCN(24,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),true);
       const Double_t* resultpar=fitter.Result().GetParams();
       const Double_t* resulterr=fitter.Result().GetErrors();
       for (unsigned int i=0;i<knri*2+6;i++){
           outparms[i]=resultpar[i];
           outparmserr[i]=resulterr[i];
       }
       treeout->Fill();
   }
   fitter.Config().SetParamsSettings(knri*2+6,parms);
   fitter.FitFCN(24,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),true);




   ROOT::Fit::FitResult result = fitter.Result();
   cout<<"\n*******PRINTING RESULT**********\n"<<endl;
   //result.Print(std::cout);
   const Double_t* resultpar=result.GetParams();
   const Double_t* resulterr=result.GetErrors();
   for (unsigned int i=0;i<result.NPar();i++){
       cout<<"p"<<i<<" = "<<resultpar[i]<<" +/- "<<resulterr[i]<<endl;
   }


   cout<<"\n*****************\n"<<endl;
   Double_t nimplants;
   sprintf(tempchar1,"treeimp");
   TTree* treeimp=(TTree*) f->Get(tempchar1);
   if (islastdssd){
       nimplants=treeimp->Draw("","","goff");
   }else{
       nimplants=treeimp->Draw("","z<5","goff");
   }


   cout<<"Number of implantations = "<<nimplants<<endl;
   cout<<"Efficiency = "<<resultpar[19]/hB->GetBinWidth(1)/resultpar[0]/nimplants*100.<<" %"<<endl;
   cout<<"t1/2\terr\tp1n\terr\tp2n\terr"<<endl;
   cout<<log(2)/resultpar[0]<<"\t"<<log(2)/resultpar[0]/resultpar[0]*resulterr[0]<<"\t"<<resultpar[9]<<"\t"<<resulterr[9]<<"\t"<<resultpar[18]<<"\t"<<resulterr[18]<<endl;
   cout<<"\n*****************\n"<<endl;
   cout<<log(2)/resultpar[0]<<"\t"<<log(2)/resultpar[0]/resultpar[0]*resulterr[0]<<endl;
   cout<<resultpar[9]<<"\t"<<resulterr[9]<<endl;
   cout<<resultpar[18]<<"\t"<<resulterr[18]<<endl;
   cout<<"\n*****************\n"<<endl;
   cout<<log(2)/resultpar[0]<<"\t"<<log(2)/resultpar[0]/resultpar[0]*resulterr[0]<<"\t";
   cout<<resultpar[9]<<"\t"<<resulterr[9]<<"\t";
   cout<<resultpar[18]<<"\t"<<resulterr[18]<<"\t"<<resultpar[19]/hB->GetBinWidth(1)/resultpar[0]/nimplants*100<<endl;
   cout<<log(2)/resultpar[0]<<"\t"<<log(2)/resultpar[0]/resultpar[0]*resulterr[0]<<"\t";
   cout<<resultpar[9]<<"\t"<<resulterr[9]<<"\t";
   cout<<resultpar[18]<<"\t"<<resulterr[18]<<"\t"<<resultpar[19]/hB->GetBinWidth(1)/resultpar[0]/nimplants*100<<endl;
   cout<<"\n*****************\n"<<endl;

   //! write to text file
   std::ofstream ofs("x2fitresult.txt", std::ofstream::out | std::ofstream::app);
   ofs<<fitname<<",";
   ofs<<log(2)/resultpar[0]<<","<<log(2)/resultpar[0]/resultpar[0]*resulterr[0]<<",";
   ofs<<resultpar[9]<<","<<resulterr[9]<<",";
   //ofs<<resultpar[18]<<","<<resulterr[18]<<","<<resultpar[19]/hB->GetBinWidth(1)/resultpar[0]/nimplants*100<<","<<binning<<","<<lowerlimit<<","<<upperlimit<<endl;
   ofs<<resultpar[18]<<","<<resulterr[18]<<","<<resultpar[19]/hB->GetBinWidth(1)/resultpar[0]/nimplants*100<<","<<binning<<","<<lowerlimit<<","<<upperlimit<<","<<resultpar[20]<<","<<resultpar[19]/resultpar[20]<<endl;


   sprintf(tempchar1,"Simfit%s",fitname);
   TCanvas * c1 = new TCanvas(tempchar1,"Simultaneous fit of two histograms",
                              10,10,1500,900);
   c1->Divide(2,3);
   c1->cd(1);
   gStyle->SetOptFit(1111);

   fB->SetFitResult( result, iparB);
   fB->SetRange(rangeB().first, rangeB().second);
   fB->SetLineColor(kBlue);
   hB->GetListOfFunctions()->Add(fB);
   hB->GetXaxis()->SetRangeUser(lowerlimit,upperlimit);
   hB->SetMarkerStyle(20);
   hB->SetMarkerSize(1);
   hB->Draw("P0");


   c1->cd(3);
   fSB->SetFitResult( result, iparSB);
   fSB->SetRange(rangeSB().first, rangeSB().second);
   fSB->SetLineColor(kRed);
   hSB->GetListOfFunctions()->Add(fSB);
   hSB->GetXaxis()->SetRangeUser(lowerlimit,upperlimit);
   hSB->SetMarkerStyle(20);
   hSB->SetMarkerSize(1);
   hSB->Draw("P0");


   c1->cd(5);
   fSB2->SetFitResult( result, iparSB2);
   fSB2->SetRange(rangeSB2().first, rangeSB2().second);
   fSB2->SetLineColor(kGreen);
   hSB2->GetListOfFunctions()->Add(fSB2);
   hSB2->Draw();
   hSB2->GetXaxis()->SetRangeUser(lowerlimit,upperlimit);
   hSB2->SetMarkerStyle(20);
   hSB2->SetMarkerSize(1);
   hSB2->Draw("P0");

   TH1F* histcomphB=new TH1F("residual_decay","residual decay",hB->GetNbinsX(),hB->GetXaxis()->GetXmin(),hB->GetXaxis()->GetXmax());
   for (Int_t i=0;i<hB->GetNbinsX();i++){
       double res;
       if (hB->GetBinError(i+1)<0.001)
       res=  0;
       else
       res=  (hB->GetBinContent(i+1)- fB->Eval( hB->GetBinCenter(i+1) ) )/hB->GetBinError(i+1);

       histcomphB->SetBinContent(i+1,res);
       histcomphB->SetBinError(i+1,1);
   }
   c1->cd(2);
   //histcomphB->Rebin(5);
   histcomphB->GetXaxis()->SetRangeUser(lowerlimit,upperlimit);
   histcomphB->Draw("e1c");


   TH1F* histcomphSB=new TH1F("residual_decay1neu","residual decay 1 neutron",hSB->GetNbinsX(),hSB->GetXaxis()->GetXmin(),hSB->GetXaxis()->GetXmax());
   for (Int_t i=0;i<hSB->GetNbinsX();i++){
       double res;
       if (hSB->GetBinError(i+1)<0.001)
        res =  0;
       else
        res =  (hSB->GetBinContent(i+1)- fSB->Eval( hSB->GetBinCenter(i+1) ) )/hSB->GetBinError(i+1);
       histcomphSB->SetBinContent(i+1,res);
       histcomphSB->SetBinError(i+1,1);
   }
   c1->cd(4);
   //histcomphSB->Rebin(5);
   histcomphSB->GetXaxis()->SetRangeUser(lowerlimit,upperlimit);
   histcomphSB->Draw("e1c");

   TH1F* histcomphSB2=new TH1F("residual_decay2neu","residual decay 2 neutrons",hSB2->GetNbinsX(),hSB2->GetXaxis()->GetXmin(),hSB2->GetXaxis()->GetXmax());
   for (Int_t i=0;i<hSB2->GetNbinsX();i++){
       double res;
       if (hSB2->GetBinError(i+1)<0.001)
        res =  0;
       else
        res =  (hSB2->GetBinContent(i+1)- fSB2->Eval( hSB2->GetBinCenter(i+1) ) )/hSB2->GetBinError(i+1);
       histcomphSB2->SetBinContent(i+1,res);
       histcomphSB2->SetBinError(i+1,1);
   }
   c1->cd(6);
   //histcomphSB2->Rebin(5);
   histcomphSB2->GetXaxis()->SetRangeUser(lowerlimit,upperlimit);
   histcomphSB2->Draw("e1c");

   fout->cd();
   c1->Write();
   treeout->Write();
   fout->Close();

   c1->Close();
}

