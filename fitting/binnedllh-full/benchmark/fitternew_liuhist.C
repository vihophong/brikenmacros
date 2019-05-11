/*Fitter with Full decay chain
 *A Linear background is adoptted
 *Nov 6, 2018
 * Update fitting the negative background
 * Plot for residual and decay components
 * Dec. 4. Added fitting code for half-life only, developed for WASABI PTEP paper
 * Fix a bug on (isparmsfix[i]>0) , replaced with (!isparmsfix[i])
 * This fitter use flat background
*/

#include "TFrame.h"
#include "TBox.h"
#include "TRatioPlot.h"
#include "makepath.C"
//Double_t neueff=0.66*(100-0.8)/100;
Double_t neueff=0.605;//changed to 62 %
Double_t neueff_mean=0.605;
Double_t neueff_err=0.053;

Bool_t reject=false;
Double_t rejectrange=0.05;//first 10 ms

//! Global Bateaman function
Double_t corefcn(Int_t ndecay,Int_t* decaymap,Int_t* nneu, Double_t* b1n,Double_t* b2n,Double_t* lamda,Double_t N0,Double_t t){
    Double_t fcnret=0;

    Double_t factor1=1.;
    //only parrent decay p2n
    for (int i=0;i<ndecay-1;i++){
        if (nneu[i]==0){
            factor1=factor1 * (1-b1n[decaymap[i]]-b2n[decaymap[i]])*lamda[decaymap[i]];
        }else if (nneu[i]==1){
            factor1=factor1 * b1n[decaymap[i]]*lamda[decaymap[i]];
        }else{
            factor1=factor1 * b2n[decaymap[i]]*lamda[decaymap[i]];
        }
    }

    Double_t factor2=0;
    for (int i=0;i<ndecay;i++){
        Double_t factor2i=exp(-lamda[decaymap[i]]*t);
        Double_t factor2ij=1;
        for (int j=0;j<ndecay;j++)
            if (j!=i) factor2ij=factor2ij*(lamda[decaymap[j]]-lamda[decaymap[i]]);
        factor2=factor2+factor2i/factor2ij;
    }

    fcnret=factor1*N0*factor2;
    return fcnret;
}

//! Global function
Double_t fcn_decay(Double_t *x, Double_t *par) {
    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t* p2n=&par[knri*2];

    //bkg
    if (x[0]<0) return par[knri*3+1]-par[knri*3+2]*x[0];
    Double_t returnval=par[knri*3+1]+par[knri*3+2]*x[0];
    //init
    Double_t N0=par[knri*3]/par[0];

    //! Parent nuclei
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0]);

    for (Int_t i=0;i<npaths;i++){
        returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0]);
    }

    /*
    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    */
    return returnval;
}
//! parent component of decay without neutron gaet
Double_t fcn_decay_parent(Double_t *x, Double_t *par) {
    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t* p2n=&par[knri*2];

    //bkg
    if (x[0]<0) return par[knri*3+1]-par[knri*3+2]*x[0];
    Double_t returnval=par[knri*3+1]+par[knri*3+2]*x[0];
    //init
    Double_t N0=par[knri*3]/par[0];

    //! Parent nuclei
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0]);

    /*
    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    */
    return returnval;
}
//! descendants components of decay without neutron gate
Double_t fcn_decay_des(Double_t *x, Double_t *par) {
    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t* p2n=&par[knri*2];

    //bkg
    if (x[0]<0) return par[knri*3+1]-par[knri*3+2]*x[0];
    Double_t returnval=par[knri*3+1]+par[knri*3+2]*x[0];
    //init
    Double_t N0=par[knri*3]/par[0];

    for (Int_t i=0;i<npaths;i++){
        returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0]);
    }

    /*
    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    */
    return returnval;
}

//! Function for decay with 1 neutron emission
Double_t fcn_1ndecay(Double_t *x, Double_t *par) {

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t* p2n=&par[knri*2];

    Double_t randcoinfgt0n=par[knri*3+3];
    Double_t randcoinf1n=par[knri*3+2];
    //bkg
    if (x[0]<0) return par[knri*3+1]-par[knri*3+4]*x[0];
    Double_t returnval=par[knri*3+1]+par[knri*3+4]*x[0];
    //init
    Double_t N0=par[knri*3]/par[0];

    //! Parent nuclei

    //! random coinc of beta decay of parent
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;
    //! decay with 1 neutron of parent
    returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf1n-randcoinfgt0n);

    //! decay with 1 neutron of parent from p2n
    returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf1n-randcoinfgt0n);
    //! decay with 2 neutron of parent (not random 1 neutron)
    returnval-=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;

    for (Int_t i=0;i<npaths;i++){
        //! random coinc of beta decay of daugter
        returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf1n;
        //! decay with 1 neutron of daugter nuclei
        returnval+=neueff*pn[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(1-randcoinf1n-randcoinfgt0n);

        //! decay with 1 neutron of daugter from p2n
        returnval+=2*(neueff*(1-neueff))*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(1-randcoinf1n-randcoinfgt0n);
        //! decay with 2 neutron of parent (not random 1 neutron)
        returnval-=neueff*neueff*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf1n;
    }

    /*
    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    */
    return returnval;

}

//! Function for parent component of the decay with 1 neutron emission
Double_t fcn_1ndecay_parent(Double_t *x, Double_t *par) {

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t* p2n=&par[knri*2];

    Double_t randcoinfgt0n=par[knri*3+3];
    Double_t randcoinf1n=par[knri*3+2];
    //bkg
    if (x[0]<0) return par[knri*3+1]-par[knri*3+4]*x[0];
    Double_t returnval=par[knri*3+1]+par[knri*3+4]*x[0];
    //init
    Double_t N0=par[knri*3]/par[0];

    //! Parent nuclei

    //! random coinc of beta decay of parent
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;
    //! decay with 1 neutron of parent
    returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf1n-randcoinfgt0n);

    //! decay with 1 neutron of parent from p2n
    returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf1n-randcoinfgt0n);
    //! decay with 2 neutron of parent (not random 1 neutron)
    returnval-=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;

}

//! Function for descendant components of decay with 1 neutron emission
Double_t fcn_1ndecay_des(Double_t *x, Double_t *par) {

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t* p2n=&par[knri*2];

    Double_t randcoinfgt0n=par[knri*3+3];
    Double_t randcoinf1n=par[knri*3+2];
    //bkg
    if (x[0]<0) return par[knri*3+1]-par[knri*3+4]*x[0];
    Double_t returnval=par[knri*3+1]+par[knri*3+4]*x[0];
    //init
    Double_t N0=par[knri*3]/par[0];

    for (Int_t i=0;i<npaths;i++){
        //! random coinc of beta decay of daugter
        returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf1n;
        //! decay with 1 neutron of daugter nuclei
        returnval+=neueff*pn[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(1-randcoinf1n-randcoinfgt0n);

        //! decay with 1 neutron of daugter from p2n
        returnval+=2*(neueff*(1-neueff))*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(1-randcoinf1n-randcoinfgt0n);
        //! decay with 2 neutron of parent (not random 1 neutron)
        returnval-=neueff*neueff*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf1n;
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;

}


//! Function for decay with 2 neutron emission
Double_t fcn_2ndecay(Double_t *x, Double_t *par) {
    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t* p2n=&par[knri*2];

    Double_t randcoinf2n=par[knri*3+4];
    Double_t randcoinfgt0n=par[knri*3+3];
    Double_t randcoinf1n=par[knri*3+2];
    //bkg
    if (x[0]<0) return par[knri*3+1]-par[knri*3+5]*x[0];
    Double_t returnval=par[knri*3+1]+par[knri*3+5]*x[0];
    //init
    Double_t N0=par[knri*3]/par[0];

    //! parent
    //! decay with 2 neutron from P2n of parent
    returnval+=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf2n-randcoinfgt0n);

    //! random coinc of beta decay of parent
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf2n;
    //! random 1n decay of parent
    returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

    //! decay with 1 neutron from P2n of parent - randomly correlated
    returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

    for (Int_t i=0;i<npaths;i++){
        //! decay with 2 neutron from P2n of daugter
        returnval+=neueff*neueff*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(1-randcoinf2n-randcoinfgt0n);
        //! random coinc of beta decay of daugter
        returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf2n;
        //! random 1n decay of daugter
        returnval+=neueff*pn[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

        //! decay with 1 neutron from P2n of daugter - randomly correlated
        returnval+=2*(neueff*(1-neueff))*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);
    }

    /*
    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    */
    return returnval;
}

//! Function for parent component of decay with 2 neutron emission
Double_t fcn_2ndecay_parent(Double_t *x, Double_t *par) {
    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t* p2n=&par[knri*2];

    Double_t randcoinf2n=par[knri*3+4];
    Double_t randcoinfgt0n=par[knri*3+3];
    Double_t randcoinf1n=par[knri*3+2];
    //bkg
    if (x[0]<0) return par[knri*3+1]-par[knri*3+5]*x[0];
    Double_t returnval=par[knri*3+1]+par[knri*3+5]*x[0];
    //init
    Double_t N0=par[knri*3]/par[0];

    //! parent
    //! decay with 2 neutron from P2n of parent
    returnval+=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf2n-randcoinfgt0n);

    //! random coinc of beta decay of parent
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf2n;
    //! random 1n decay of parent
    returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

    //! decay with 1 neutron from P2n of parent - randomly correlated
    returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);


    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

//! Function for descendant components of decay with 2 neutron emission
Double_t fcn_2ndecay_des(Double_t *x, Double_t *par) {
    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t* p2n=&par[knri*2];

    Double_t randcoinf2n=par[knri*3+4];
    Double_t randcoinfgt0n=par[knri*3+3];
    Double_t randcoinf1n=par[knri*3+2];
    //bkg
    if (x[0]<0) return par[knri*3+1]-par[knri*3+5]*x[0];
    Double_t returnval=par[knri*3+1]+par[knri*3+5]*x[0];
    //init
    Double_t N0=par[knri*3]/par[0];

    for (Int_t i=0;i<npaths;i++){
        //! decay with 2 neutron from P2n of daugter
        returnval+=neueff*neueff*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(1-randcoinf2n-randcoinfgt0n);
        //! random coinc of beta decay of daugter
        returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf2n;
        //! random 1n decay of daugter
        returnval+=neueff*pn[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

        //! decay with 1 neutron from P2n of daugter - randomly correlated
        returnval+=2*(neueff*(1-neueff))*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);
    }
    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

// Definition of shared parameter

// Decay parameters
int iparB[kmaxparms];

// Decay 1 neutron parameters
int iparSB[kmaxparms];

// Decay 2 neutrons parameters
int iparSB2[kmaxparms];

// Create the GlobalCHi2 structure

struct GlobalChi2 {
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2,
                ROOT::Math::IMultiGenFunction & f3) :
      fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}

   // parameter vector is first background (in common 1 and 2)
   // and then is signal (only in 2)
   double operator() (const double *par) const {
      double p1[knri*3+3];
      for (int i = 0; i < knri*3+3; ++i) p1[i] = par[iparB[i] ];

      double p2[knri*3+5];
      for (int i = 0; i < knri*3+5; ++i) p2[i] = par[iparSB[i] ];

      double p3[knri*3+6];
      for (int i = 0; i < knri*3+6; ++i) p3[i] = par[iparSB2[i] ];

      return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3);
   }
   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
   const  ROOT::Math::IMultiGenFunction * fChi2_3;
};



//! Monte carlo variation
void mc(TRandom3* rseed)
{
    for (int i=0;i<knri*3;i++){
        if (isparmsfix[i]==1&&parmserr[i]>0.){
            mcparms[i]=rseed->Gaus(parms[i],parmserr[i]);
        }
    }
    //! mc for neutron efficiency?
    if (neueff_err>0){// only if there is error
      //neueff=neueff_mean-neueff_err+rseed->Rndm()*neueff_err*2;//for an uniform distribution
      neueff=rseed->Gaus(neueff_mean,neueff_err);//for a gausian distribution
    }
}


void fitter(char* infilename,char* parmsfilename,char* outfilename,Double_t neueff_i=0.605, Double_t neuefferr_i=0.053, Int_t ninterations=0,Int_t rebin=1)
{
    neueff_err=neuefferr_i;
    neueff_mean=neueff_i;
    neueff=neueff_i;

    Double_t lowerlimit=-5;
    Double_t upperlimit=5;

    //! input decay parameters and make decay path
    makepath(parmsfilename);

    TRandom3* rseed=new TRandom3;
    //! construct params

    Double_t bkgactmaxmin=0.50; //100% of max min bkg or initial activity
    Double_t plotrange[]={0.05,5.,-2};

    //! input default value for other parameters
    parms[knri*3]=1000;//initial activity
    parms[knri*3+1]=2000;//background of no neutron gate curve
    parms[knri*3+2]=500;//background of 1n gate curve
    parms[knri*3+3]=100;//background of 2n gate curve

    parms[knri*3+4]=0.02;//random 1 neutron factor
    parms[knri*3+5]=0.022;//random gt0 neutron factor
    parms[knri*3+6]=0.005;//random 2 neutrons factor

    parms[knri*3+7]=0.;//background slope of no neutron gate curve
    parms[knri*3+8]=0.;//background slope of 1n gate curve
    parms[knri*3+9]=0.;//background slope of 2n gate curve

    parmsmin[knri*3]=0;parmsmax[knri*3]=10000;
    parmsmin[knri*3+1]=0;parmsmax[knri*3+1]=10000;
    parmsmin[knri*3+2]=0;parmsmax[knri*3+2]=10000;
    parmsmin[knri*3+3]=0;parmsmax[knri*3+3]=10000;

    parmsmin[knri*3+4]=0;parmsmax[knri*3+4]=1;
    parmsmin[knri*3+5]=0;parmsmax[knri*3+5]=1;
    parmsmin[knri*3+6]=0;parmsmax[knri*3+6]=1;

    parmsmin[knri*3+7]=-1000;parmsmax[knri*3+7]=1000;
    parmsmin[knri*3+8]=-1000;parmsmax[knri*3+8]=1000;
    parmsmin[knri*3+9]=-1000;parmsmax[knri*3+9]=1000;

    isparmsfix[knri*3]=0;//inital activity - vary
    isparmsfix[knri*3+1]=1;//background of no neutron gate curve ?fix
    isparmsfix[knri*3+2]=1;//background of 1n gate curve ?fix
    isparmsfix[knri*3+3]=1;//background of 2n gate curve ?fix

    isparmsfix[knri*3+4]=1;
    isparmsfix[knri*3+5]=1;
    isparmsfix[knri*3+6]=1;

    isparmsfix[knri*3+7]=1;//background slope of no neutron gate curve
    isparmsfix[knri*3+8]=1;//background slope of 1n gate curve
    isparmsfix[knri*3+9]=1;//background slope of 2n gate curve

    //!****************************GET HISTOGRAM FROM FILE********************
    //!
    //!
    TFile *f = TFile::Open(infilename);
    char tempchar1[1000];
    sprintf(tempchar1,"decaytime");
    TH1F* hdecay=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"decaytimen1");
    TH1F* hdecay1n=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"decaytimen2");
    TH1F* hdecay2n=(TH1F*) gDirectory->Get(tempchar1);

    sprintf(tempchar1,"decaytimen1_min_n1");
    TH1F* hdecay1nbwd=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"decaytimen2_min_n2");
    TH1F* hdecay2nbwd=(TH1F*) gDirectory->Get(tempchar1);

    TH1F* hdecaygt0nbwd=(TH1F*) hdecay1nbwd->Clone();
    hdecaygt0nbwd->Add(hdecay2nbwd);
    hdecaygt0nbwd->SetName("Sn136ImplantBetagt0nBW");



    //! Calculate random coincidence paramters
    Double_t n1nbwd=(Double_t) hdecay1nbwd->GetEntries();
    Double_t gt0nbwd=(Double_t) hdecaygt0nbwd->GetEntries();
    Double_t n2nbwd=(Double_t) hdecay2nbwd->GetEntries();
    Double_t nball=(Double_t) hdecay->GetEntries();

    parms[knri*3+4]=n1nbwd/nball;
    parms[knri*3+5]=gt0nbwd/nball;
    parms[knri*3+6]=n2nbwd/nball;

    //! Binning
    Int_t binning=hdecay->GetNbinsX();
    //!Rebin
    hdecay->Rebin(rebin);
    hdecay1n->Rebin(rebin);
    hdecay2n->Rebin(rebin);
    TH1F * hB = (TH1F*) hdecay->Clone();
    TH1F * hSB = (TH1F*) hdecay1n->Clone();
    TH1F * hSB2 = (TH1F*) hdecay2n->Clone();

    //! background parameter estimation
    //! background as average of several first bin

    //background of decay no neutron gate

    hdecay->Fit("pol0","LQRE0","goff",lowerlimit,0.);

    parms[knri*3+1]=hdecay->GetFunction("pol0")->GetParameter(0);
    parmserr[knri*3+1]=hdecay->GetFunction("pol0")->GetParError(0);
    parmsmin[knri*3+1]=parms[knri*3+1]-parms[knri*3+1]*bkgactmaxmin;
    parmsmax[knri*3+1]=parms[knri*3+1]+parms[knri*3+1]*bkgactmaxmin;

    //parms[knri*3+7]=-hdecay->GetFunction("pol1")->GetParameter(1);
    //parmserr[knri*3+7]=-hdecay->GetFunction("pol1")->GetParError(1);

    //activity
    parms[knri*3]=hdecay->GetBinContent(hdecay->GetXaxis()->FindBin(rejectrange))-parms[knri*3+1];
    parmserr[knri*3]=parms[knri*3]*bkgactmaxmin;
    parmsmin[knri*3]=parms[knri*3]-parms[knri*3]*bkgactmaxmin;
    parmsmax[knri*3]=parms[knri*3]*2+parms[knri*3]*bkgactmaxmin;

    //background of decay 1 neutron gate

    hdecay1n->Fit("pol0","LQRE0","goff",lowerlimit,0.);

    parms[knri*3+2]=hdecay1n->GetFunction("pol0")->GetParameter(0);
    parmserr[knri*3+2]=hdecay1n->GetFunction("pol0")->GetParError(0);
    parmsmin[knri*3+2]=parms[knri*3+2]-parms[knri*3+2]*bkgactmaxmin;
    parmsmax[knri*3+2]=parms[knri*3+2]+parms[knri*3+2]*bkgactmaxmin;

    //parms[knri*3+8]=-hdecay1n->GetFunction("pol1")->GetParameter(1);
    //parmserr[knri*3+8]=-hdecay1n->GetFunction("pol1")->GetParError(1);

    //background of decay 2 neutrons gate

    hdecay2n->Fit("pol0","LQRE0","goff",lowerlimit,0.);

    parms[knri*3+3]=hdecay2n->GetFunction("pol0")->GetParameter(0);
    parmserr[knri*3+3]=hdecay2n->GetFunction("pol0")->GetParError(0);
    parmsmin[knri*3+3]=parms[knri*3+3]-parms[knri*3+3]*bkgactmaxmin;
    parmsmax[knri*3+3]=parms[knri*3+3]+parms[knri*3+3]*bkgactmaxmin;

    //parms[knri*3+9]=-hdecay2n->GetFunction("pol1")->GetParameter(1);
    //parmserr[knri*3+9]=-hdecay2n->GetFunction("pol1")->GetParError(1);

    //! ********** Define FITTING FUNCTION
    //! Define function without neutron gate
    TF1* fB=new TF1("fB",fcn_decay,lowerlimit,upperlimit,knri*3+3);
    fB->SetNpx(2000);
    fB->SetLineWidth(2);
    fB->SetLineColor(8);

    //! initializing parameters
    for (Int_t i=0;i<knri*3;i++){
        fB->SetParameter(i,parms[i]);
        fB->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (isparmsfix[i]>0) fB->FixParameter(i,parms[i]);
    }
    fB->SetParameter(knri*3,parms[knri*3]);//inital activity
    fB->SetParLimits(knri*3,parmsmin[knri*3],parmsmax[knri*3]);

    fB->SetParameter(knri*3+1,parms[knri*3+1]);//background
    fB->SetParLimits(knri*3+1,parmsmin[knri*3+1],parmsmax[knri*3+1]);

    fB->SetParameter(knri*3+2,parms[knri*3+7]);//background slope
    fB->SetParLimits(knri*3+2,parmsmin[knri*3+7],parmsmax[knri*3+7]);

    //! Define function with 1 neutron gate
    TF1* fSB=new TF1("fSB",fcn_1ndecay,lowerlimit,upperlimit,knri*3+5);
    fSB->SetNpx(2000);
    fSB->SetLineWidth(2);
    fSB->SetLineColor(8);

    // initializing parameters
    for (Int_t i=0;i<knri*3;i++){
        fSB->SetParameter(i,parms[i]);
        fSB->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (isparmsfix[i]>0) fSB->FixParameter(i,parms[i]);
    }
    fSB->SetParameter(knri*3,parms[knri*3]);//inital activity
    fSB->SetParLimits(knri*3,parmsmin[knri*3],parmsmax[knri*3]);

    fSB->SetParameter(knri*3+1,parms[knri*3+2]);//background
    fSB->SetParLimits(knri*3+1,parmsmin[knri*3+2],parmsmax[knri*3+2]);

    fSB->SetParameter(knri*3+2,parms[knri*3+4]);//random 1 neutron factor
    fSB->SetParLimits(knri*3+2,parmsmin[knri*3+4],parmsmax[knri*3+4]);

    fSB->SetParameter(knri*3+3,parms[knri*3+5]);//random gt0 neutron factor
    fSB->SetParLimits(knri*3+3,parmsmin[knri*3+5],parmsmax[knri*3+5]);

    fSB->SetParameter(knri*3+4,parms[knri*3+8]);//background slope
    fSB->SetParLimits(knri*3+4,parmsmin[knri*3+8],parmsmax[knri*3+8]);

    //! Define function with 2 neutrons gate
    TF1* fSB2=new TF1("fSB2",fcn_2ndecay,lowerlimit,upperlimit,knri*3+6);
    fSB2->SetNpx(2000);
    fSB2->SetLineWidth(2);
    fSB2->SetLineColor(8);

    // initializing parameters
    for (Int_t i=0;i<knri*3;i++){
        fSB2->SetParameter(i,parms[i]);
        fSB2->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (isparmsfix[i]>0) fSB2->FixParameter(i,parms[i]);
    }
    fSB2->SetParameter(knri*3,parms[knri*3]);//inital activity
    fSB2->SetParLimits(knri*3,parmsmin[knri*3],parmsmax[knri*3]);

    fSB2->SetParameter(knri*3+1,parms[knri*3+3]);//background
    fSB2->SetParLimits(knri*3+1,parmsmin[knri*3+3],parmsmax[knri*3+3]);

    fSB2->SetParameter(knri*3+2,parms[knri*3+4]);//random 1 neutron factor
    fSB2->SetParLimits(knri*3+2,parmsmin[knri*3+4],parmsmax[knri*3+4]);

    fSB2->SetParameter(knri*3+3,parms[knri*3+5]);//random gt0 neutron factor
    fSB2->SetParLimits(knri*3+3,parmsmin[knri*3+5],parmsmax[knri*3+5]);

    fSB2->SetParameter(knri*3+4,parms[knri*3+6]);//random 2 neutrons factor
    fSB2->SetParLimits(knri*3+4,parmsmin[knri*3+6],parmsmax[knri*3+6]);

    fSB2->SetParameter(knri*3+5,parms[knri*3+9]);//background slope
    fSB2->SetParLimits(knri*3+5,parmsmin[knri*3+9],parmsmax[knri*3+9]);

    cout<<fB->Eval(1.)<<endl;
    cout<<fSB->Eval(1.)<<endl;
    cout<<fSB2->Eval(1.)<<endl;


    //! Book output file and tree output for systematic half-life estimation
    TFile* outfile=new TFile(outfilename,"recreate");

    Double_t outparms[knri*3+10];
    Double_t outparmserr[knri*3+10];
    Int_t iparms[knri*3+10];
    Int_t isvary[knri*3+10];
    for (int i=0;i<knri*3+10;i++){
        outparms[i]=0;
        outparmserr[i]=0;
        iparms[i]=i;
        isvary[i]=isparmsfix[i];
    }
    sprintf(tempchar1,"treemc");
    TTree* treeout=new TTree(tempchar1,tempchar1);
    treeout->Branch("outparms",outparms,Form("outparms[%d]/D",knri*3+10));
    treeout->Branch("outparmserr",outparmserr,Form("outparmserr[%d]/D",knri*3+10));
    treeout->Branch("iparms",iparms,Form("iparms[%d]/I",knri*3+10));
    treeout->Branch("neueff",&neueff,"neueff/D");
    treeout->Branch("isvary",outparmserr,Form("isvary[%d]/I",knri*3+10));


    //! Define Simultaneous fitting function
    //! initialzing shared prameter matrix for simultaneous fitting
    for (int i = 0; i < knri*3; i++) {
        iparB[i]=i;
        iparSB[i]=i;
        iparSB2[i]=i;
    }
    iparB[knri*3]=knri*3;//initial activity
    iparB[knri*3+1]=knri*3+1;//background of no neutron gate curve
    iparB[knri*3+2]=knri*3+7;//background slope of no neutron gate curve

    iparSB[knri*3]=knri*3;//initial activity
    iparSB[knri*3+1]=knri*3+2;//background of 1n gate curve
    iparSB[knri*3+4]=knri*3+8;//background slope of 1n gate curve

    iparSB2[knri*3]=knri*3;//initial activity
    iparSB2[knri*3+1]=knri*3+3;//background of 1n gate curve
    iparSB2[knri*3+5]=knri*3+9;//background slope of 1n gate curve

    iparSB[knri*3+2]=knri*3+4;//random 1 neutron factor
    iparSB[knri*3+3]=knri*3+5;//random gt0 neutron factor

    iparSB2[knri*3+2]=knri*3+4;//random 1 neutron factor
    iparSB2[knri*3+3]=knri*3+5;//random gt0 neutron factor
    iparSB2[knri*3+4]=knri*3+6;//random 2 neutrons factor

    ROOT::Math::WrappedMultiTF1 wfB(*fB,1);
    ROOT::Math::WrappedMultiTF1 wfSB(*fSB,1);
    ROOT::Math::WrappedMultiTF1 wfSB2(*fSB2,1);
    ROOT::Fit::DataOptions opt;

    // limit within the fitting range
    opt.fUseRange  =true;
    // set the data range
    ROOT::Fit::DataRange rangeB;
    rangeB.SetRange(rejectrange,upperlimit);
    ROOT::Fit::BinData dataB(opt,rangeB);
    ROOT::Fit::FillData(dataB, hB);

    ROOT::Fit::DataRange rangeSB;
    rangeSB.SetRange(rejectrange,upperlimit);
    ROOT::Fit::BinData dataSB(opt,rangeSB);
    ROOT::Fit::FillData(dataSB, hSB);

    ROOT::Fit::DataRange rangeSB2;
    rangeSB2.SetRange(rejectrange,upperlimit);
    ROOT::Fit::BinData dataSB2(opt,rangeSB2);
    ROOT::Fit::FillData(dataSB2, hSB2);

    ROOT::Fit::PoissonLLFunction chi2_B(dataB, wfB);
    ROOT::Fit::PoissonLLFunction chi2_SB(dataSB, wfSB);
    ROOT::Fit::PoissonLLFunction chi2_SB2(dataSB2, wfSB2);
    GlobalChi2 globalChi2(chi2_B, chi2_SB, chi2_SB2);

    ROOT::Fit::Fitter fitter;

    //! ***********SETTING PARAMETERS************

    cout<<"\n\n*****SETTING PARAMETERS......\n"<<endl;
    fitter.Config().SetParamsSettings(knri*3+10,parms);
    for (int i=0;i<knri*3+10;i++){
        fitter.Config().ParSettings(i).SetLimits(parmsmin[i],parmsmax[i]);
        if (isparmsfix[i]>0){
           cout<<"fixed value p"<<i<<"\tval="<<parms[i]<<endl;
           fitter.Config().ParSettings(i).Fix();
        }else{
           cout<<"variable value p"<<i<<"\tval "<<parms[i]<<"\tmin="<<parmsmin[i]<<"\tmax="<<parmsmax[i]<<endl;
        }
    }

    //! Setting MINOS package
    fitter.Config().SetMinimizer("Minuit2","Migrad");
    //fitter.Config().SetMinosErrors();
    if (fitter.Config().MinosErrors()) cout<<"minos enabled"<<endl;

    //!Perform the fit
    fitter.FitFCN(knri*3+10,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),false);


    //! Printting Results
    ROOT::Fit::FitResult result = fitter.Result();
    const Double_t* resultpar=result.GetParams();
    const Double_t* resulterr=result.GetErrors();


    //! fix monte carlo parameter center
    for (int i=0;i<knri*3+10;i++) mcparms[i]=parms[i];
    for(int i=0;i<ninterations;i++) {
        cout<<"mc "<<i+1<<endl;
        mc(rseed);
        fitter.Config().SetParamsSettings(knri*3+10,mcparms);
        for (unsigned int j=0;j<knri*3+10;j++){
            fitter.Config().ParSettings(j).SetLimits(parmsmin[j],parmsmax[j]);
            if (isparmsfix[j]>0){
               fitter.Config().ParSettings(j).Fix();
            }
        }

        fitter.FitFCN(knri*3+10,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),false);

        const Double_t* resultparmc=fitter.Result().GetParams();
        const Double_t* resulterrmc=fitter.Result().GetErrors();
        for (unsigned int i=0;i<knri*3+10;i++){
            outparms[i]=resultparmc[i];
            outparmserr[i]=resulterrmc[i];
        }
        treeout->Fill();
    }

    //!Perform the fit again
    //fitter.Config().SetParamsSettings(knri*3+10,parms);
    //fitter.FitFCN(knri*3+10,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),false);

    //! Filling the output tree
    for (unsigned int i=0;i<knri*3+10;i++){
      outparms[i]=resultpar[i];
      outparmserr[i]=resulterr[i];
    }
    treeout->Fill();



    cout<<"\n*******PRINTING RESULT**********\n"<<endl;
    result.Print(std::cout);

    cout<<"\n*******PRINTING MAIN RESULTS**********\n"<<endl;
    cout<<"t1/2\terr\tp1n\terr\tp2n\terr"<<endl;
    cout<<log(2)/resultpar[0]<<"\t"<<log(2)/resultpar[0]/resultpar[0]*resulterr[0]<<"\t"<<resultpar[knri]<<"\t"<<resulterr[knri]<<"\t"<<resultpar[knri*2]<<"\t"<<resulterr[knri*2]<<endl;


    cout<<"\n*****************\n"<<endl;
    cout<<log(2)/resultpar[0]<<"\t"<<log(2)/resultpar[0]/resultpar[0]*resulterr[0]<<"\t";
    cout<<resultpar[knri]<<"\t"<<resulterr[knri]<<"\t";
    cout<<resultpar[knri*2]<<"\t"<<resulterr[knri*2]<<endl;
    cout<<log(2)/resultpar[0]<<"\t"<<log(2)/resultpar[0]/resultpar[0]*resulterr[0]<<"\t";
    cout<<resultpar[knri]<<"\t"<<resulterr[knri]<<"\t";
    cout<<resultpar[knri*2]<<"\t"<<resulterr[knri*2]<<endl;


    //! write to text file
    std::ofstream ofs("fitresults.txt", std::ofstream::out | std::ofstream::app);
    ofs<<parmsfilename<<"\t"<<log(2)/resultpar[0]<<"\t"<<log(2)/resultpar[0]/resultpar[0]*resulterr[0]<<"\t"<<resultpar[knri]<<"\t"<<resulterr[knri]<<"\t"<<resultpar[knri*2]<<"\t"<<resulterr[knri*2]<<endl;


    TH1F* histcomphB=new TH1F("residual_decay","Residual of decay curve",hB->GetNbinsX(),hB->GetXaxis()->GetXmin(),hB->GetXaxis()->GetXmax());
    for (Int_t i=0;i<hB->GetNbinsX();i++){
        double res;
        if (hB->GetBinError(i+1)<0.001)
          res=  0;
        else
        res=  (hB->GetBinContent(i+1)- fB->Eval( hB->GetBinCenter(i+1) ) )/hB->GetBinError(i+1);

        histcomphB->SetBinContent(i+1,res);
        histcomphB->SetBinError(i+1,0);
    }

    TH1F* histcomphSB=new TH1F("residual_decay1neu","Residual of decay curve with one neutron gate",hSB->GetNbinsX(),hSB->GetXaxis()->GetXmin(),hSB->GetXaxis()->GetXmax());
    for (Int_t i=0;i<hSB->GetNbinsX();i++){
        double res;
        if (hSB->GetBinError(i+1)<0.001)
         res =  0;
        else
         res =  (hSB->GetBinContent(i+1)- fSB->Eval( hSB->GetBinCenter(i+1) ) )/hSB->GetBinError(i+1);
        histcomphSB->SetBinContent(i+1,res);
        histcomphSB->SetBinError(i+1,0);
    }

    TH1F* histcomphSB2=new TH1F("residual_decay2neu","Residual of decay curve with two neutron gate",hSB2->GetNbinsX(),hSB2->GetXaxis()->GetXmin(),hSB2->GetXaxis()->GetXmax());
    for (Int_t i=0;i<hSB2->GetNbinsX();i++){
        double res;
        if (hSB2->GetBinError(i+1)<0.001)
         res =  0;
        else
         res =  (hSB2->GetBinContent(i+1)- fSB2->Eval( hSB2->GetBinCenter(i+1) ) )/hSB2->GetBinError(i+1);
        histcomphSB2->SetBinContent(i+1,res);
        histcomphSB2->SetBinError(i+1,0);
    }

    //! Input descendant decay
    TF1* fbBparent=new TF1("fbBparent",fcn_decay_parent,plotrange[0],plotrange[1],knri*3+10);
    TF1* fbBdes=new TF1("fbBdes",fcn_decay_des,plotrange[0],plotrange[1],knri*3+10);

    TF1* fbSBparent=new TF1("fbSBparent",fcn_1ndecay_parent,plotrange[0],plotrange[1],knri*3+10);
    TF1* fbSBdes=new TF1("fbSBdes",fcn_1ndecay_des,plotrange[0],plotrange[1],knri*3+10);

    TF1* fbSB2parent=new TF1("fbSB2parent",fcn_2ndecay_parent,plotrange[0],plotrange[1],knri*3+10);
    TF1* fbSB2des=new TF1("fbSB2des",fcn_2ndecay_des,plotrange[0],plotrange[1],knri*3+10);


    TCanvas* c1=new TCanvas("c1","c1",900,1200);

    gStyle->SetOptStat(11111);
    c1->Divide(1,3);
    c1->cd(1);
    fB->SetFitResult( result, iparB);
    fB->SetRange(plotrange[0], plotrange[1]);
    fB->SetLineColor(kRed);
    fB->SetNpx(binning);

    hB->SetTitle("Decay curve");
    hB->GetListOfFunctions()->Add(fB);
    hB->GetXaxis()->SetRangeUser(plotrange[2],plotrange[1]);
    hB->GetXaxis()->SetTitle("Time (s)");
    hB->GetYaxis()->SetTitle("Counts");
    hB->GetXaxis()->SetTitleSize(0.05);
    hB->GetYaxis()->SetTitleSize(0.05);
    hB->SetMarkerStyle(20);
    hB->SetMarkerSize(0.8);
    hB->GetXaxis()->SetLabelSize(0.05);
    hB->GetYaxis()->SetLabelSize(0.05);
    hB->Draw("P0 E");

    fbBparent->SetFitResult( result, iparB);
    fbBparent->SetLineColor(kBlue);
    fbBparent->SetNpx(binning);
    fbBparent->Draw("same");

    fbBdes->SetFitResult( result, iparB);
    fbBdes->SetLineColor(kGreen);
    fbBdes->SetNpx(binning);
    fbBdes->Draw("same");

    hdecay->GetFunction("pol0")->SetLineColor(kYellow);
    hdecay->GetFunction("pol0")->Draw("same");

    c1->cd(2);


    fSB->SetFitResult( result, iparSB);
    fSB->SetRange(plotrange[0], plotrange[1]);
    fSB->SetLineColor(kRed);
    fSB->SetNpx(binning);

    hSB->SetTitle("Decay curve with one neutron gate");
    hSB->GetListOfFunctions()->Add(fSB);
    hSB->GetXaxis()->SetRangeUser(plotrange[2],plotrange[1]);
    hSB->GetXaxis()->SetTitle("Time (s)");
    hSB->GetYaxis()->SetTitle("Counts");
    hSB->GetXaxis()->SetTitleSize(0.05);
    hSB->GetYaxis()->SetTitleSize(0.05);
    hSB->SetMarkerStyle(20);
    hSB->SetMarkerSize(0.8);
    hSB->GetXaxis()->SetLabelSize(0.05);
    hSB->GetYaxis()->SetLabelSize(0.05);
    hSB->Draw("P0 E");

    fbSBparent->SetFitResult( result, iparSB);
    fbSBparent->SetLineColor(kBlue);
    fbSBparent->SetNpx(binning);
    fbSBparent->Draw("same");

    fbSBdes->SetFitResult( result, iparSB);
    fbSBdes->SetLineColor(kGreen);
    fbSBdes->SetNpx(binning);
    fbSBdes->Draw("same");

    hdecay1n->GetFunction("pol0")->SetLineColor(kYellow);
    hdecay1n->GetFunction("pol0")->Draw("same");

    c1->cd(3);
    fSB2->SetFitResult( result, iparSB2);
    fSB2->SetRange(plotrange[0], plotrange[1]);
    fSB2->SetLineColor(kRed);
    fSB2->SetNpx(binning);
    hSB2->GetXaxis()->SetTitle("Time (s)");
    hSB2->GetYaxis()->SetTitle("Counts");
    hSB2->GetXaxis()->SetTitleSize(0.05);
    hSB2->GetYaxis()->SetTitleSize(0.05);
    hSB2->GetXaxis()->SetLabelSize(0.05);
    hSB2->GetYaxis()->SetLabelSize(0.05);
    hSB2->SetTitle("Decay curve with two neutron gate");
    hSB2->GetListOfFunctions()->Add(fSB2);
    hSB2->GetXaxis()->SetRangeUser(plotrange[2],plotrange[1]);

    hSB2->SetMarkerStyle(20);
    hSB2->SetMarkerSize(0.8);
    hSB2->GetXaxis()->SetLabelSize(0.05);
    hSB2->GetYaxis()->SetLabelSize(0.05);
    hSB2->Draw("P0 E");

    hdecay2n->GetFunction("pol0")->SetLineColor(kYellow);
    hdecay2n->GetFunction("pol0")->Draw("same");

    hB->Write();
    fB->Write();
    fbBparent->Write();
    fbBdes->Write();
    histcomphB->Write();

    hSB->Write();
    fSB->Write();
    fbSBparent->Write();
    fbSBdes->Write();
    histcomphSB->Write();

    hSB2->Write();
    fSB2->Write();
    fbSB2parent->Write();
    fbSB2des->Write();
    histcomphSB2->Write();

    c1->Write();


    outfile->Close();
}


void fitterHalfLifeOnly(char* infilename,char* parmsfilename,char* outfilename,Int_t ninterations=0,Int_t rebin=1)
{
    Double_t lowerlimit=-10;
    Double_t upperlimit=10;

    //! input decay parameters and make decay path
    makepath(parmsfilename);

    TRandom3* rseed=new TRandom3;
    //! construct params

    Double_t bkgactmaxmin=0.50; //100% of max min bkg or initial activity
    Double_t plotrange[]={rejectrange,10.,-10};

    //! input default value for other parameters
    parms[knri*3]=1000;//initial activity
    parms[knri*3+1]=2000;//background of no neutron gate curve
    parms[knri*3+2]=500;//background of 1n gate curve
    parms[knri*3+3]=100;//background of 2n gate curve

    parms[knri*3+4]=0.02;//random 1 neutron factor
    parms[knri*3+5]=0.022;//random gt0 neutron factor
    parms[knri*3+6]=0.005;//random 2 neutrons factor

    parms[knri*3+7]=0.;//background slope of no neutron gate curve
    parms[knri*3+8]=0.;//background slope of 1n gate curve
    parms[knri*3+9]=0.;//background slope of 2n gate curve

    parmsmin[knri*3]=0;parmsmax[knri*3]=10000;
    parmsmin[knri*3+1]=0;parmsmax[knri*3+1]=10000;
    parmsmin[knri*3+2]=0;parmsmax[knri*3+2]=10000;
    parmsmin[knri*3+3]=0;parmsmax[knri*3+3]=10000;

    parmsmin[knri*3+4]=0;parmsmax[knri*3+4]=1;
    parmsmin[knri*3+5]=0;parmsmax[knri*3+5]=1;
    parmsmin[knri*3+6]=0;parmsmax[knri*3+6]=1;

    parmsmin[knri*3+7]=-1000;parmsmax[knri*3+7]=1000;
    parmsmin[knri*3+8]=-1000;parmsmax[knri*3+8]=1000;
    parmsmin[knri*3+9]=-1000;parmsmax[knri*3+9]=1000;

    isparmsfix[knri*3]=0;//inital activity - vary
    isparmsfix[knri*3+1]=0;//background of no neutron gate curve ?fix
    isparmsfix[knri*3+2]=1;//background of 1n gate curve ?fix
    isparmsfix[knri*3+3]=1;//background of 2n gate curve ?fix

    isparmsfix[knri*3+4]=1;
    isparmsfix[knri*3+5]=1;
    isparmsfix[knri*3+6]=1;

    isparmsfix[knri*3+7]=0;//background slope of no neutron gate curve
    isparmsfix[knri*3+8]=1;//background slope of 1n gate curve
    isparmsfix[knri*3+9]=1;//background slope of 2n gate curve

    //!****************************GET HISTOGRAM FROM FILE********************
    //!
    //!
    TFile *f = TFile::Open(infilename);
    char tempchar1[1000];
    sprintf(tempchar1,"hdecay");
    TH1F* hdecay=(TH1F*) gDirectory->Get(tempchar1);

    //! Binning
    Int_t binning=hdecay->GetNbinsX();
    //!Rebin
    hdecay->Rebin(rebin);
    TH1F * hB = (TH1F*) hdecay->Clone();

    //! background parameter estimation
    //! background as average of several first bin

    //background of decay no neutron gate

    hdecay->Fit("pol1","LQRE0","goff",lowerlimit,0.);

    parms[knri*3+1]=hdecay->GetFunction("pol1")->GetParameter(0);
    parmserr[knri*3+1]=hdecay->GetFunction("pol1")->GetParError(0);
    parmsmin[knri*3+1]=parms[knri*3+1]-parms[knri*3+1]*bkgactmaxmin;
    parmsmax[knri*3+1]=parms[knri*3+1]+parms[knri*3+1]*bkgactmaxmin;

    parms[knri*3+7]=-hdecay->GetFunction("pol1")->GetParameter(1);
    parmserr[knri*3+7]=-hdecay->GetFunction("pol1")->GetParError(1);

    //activity
    parms[knri*3]=hdecay->GetBinContent(hdecay->GetXaxis()->FindBin(rejectrange))-parms[knri*3+1];
    parmserr[knri*3]=parms[knri*3]*bkgactmaxmin;
    parmsmin[knri*3]=parms[knri*3]-parms[knri*3]*bkgactmaxmin;
    parmsmax[knri*3]=parms[knri*3]*2+parms[knri*3]*bkgactmaxmin;

    //! ********** Define FITTING FUNCTION
    //! Define function without neutron gate
    TF1* fB=new TF1("fB",fcn_decay,lowerlimit,upperlimit,knri*3+3);
    fB->SetNpx(2000);
    fB->SetLineWidth(2);
    fB->SetLineColor(8);

    //! initializing parameters
    for (Int_t i=0;i<knri*3;i++){
        fB->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (isparmsfix[i]>0) fB->FixParameter(i,parms[i]);
        else fB->SetParameter(i,parms[i]);
    }
    fB->SetParameter(knri*3,parms[knri*3]);//inital activity
    fB->SetParLimits(knri*3,parmsmin[knri*3],parmsmax[knri*3]);

    fB->SetParLimits(knri*3+1,parmsmin[knri*3+1],parmsmax[knri*3+1]);
    fB->FixParameter(knri*3+1,parms[knri*3+1]);//background
    //fB->SetParameter(knri*3+1,parms[knri*3+1]);//background

    fB->SetParLimits(knri*3+2,parmsmin[knri*3+7],parmsmax[knri*3+7]);
    fB->FixParameter(knri*3+2,parms[knri*3+7]);//background slope
    //fB->SetParameter(knri*3+2,parms[knri*3+7]);//background slope

    cout<<fB->Eval(1.)<<endl;


    //! Book output file and tree output for systematic half-life estimation
    TFile* outfile=new TFile(outfilename,"recreate");

    Double_t outparms[knri*3+10];
    Double_t outparmserr[knri*3+10];
    Int_t iparms[knri*3+10];
    Int_t isvary[knri*3+10];
    for (int i=0;i<knri*3+10;i++){
        outparms[i]=0;
        outparmserr[i]=0;
        iparms[i]=i;
        isvary[i]=isparmsfix[i];
    }
    sprintf(tempchar1,"treemc");
    TTree* treeout=new TTree(tempchar1,tempchar1);
    treeout->Branch("outparms",outparms,Form("outparms[%d]/D",knri*3+10));
    treeout->Branch("outparmserr",outparmserr,Form("outparmserr[%d]/D",knri*3+10));
    treeout->Branch("iparms",iparms,Form("iparms[%d]/I",knri*3+10));
    treeout->Branch("neueff",&neueff,"neueff/D");
    treeout->Branch("isvary",outparmserr,Form("isvary[%d]/I",knri*3+10));


    //!Perform the fit
    TFitResult* resul=hB->Fit(fB,"LVRE0","goff",rejectrange,upperlimit).Get();

    /*
    //! fix monte carlo parameter center
    for (int i=0;i<knri*3+10;i++) mcparms[i]=parms[i];
    for(int i=0;i<ninterations;i++) {
        cout<<"mc "<<i+1<<endl;
        mc(rseed);
        fitter.Config().SetParamsSettings(knri*3+10,mcparms);
        for (unsigned int j=0;j<knri*3+10;j++){
            fitter.Config().ParSettings(j).SetLimits(parmsmin[j],parmsmax[j]);
            if (isparmsfix[j]){
               fitter.Config().ParSettings(j).Fix();
            }
        }

        fitter.FitFCN(knri*3+10,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),false);

        const Double_t* resultparmc=fitter.Result().GetParams();
        const Double_t* resulterrmc=fitter.Result().GetErrors();
        for (unsigned int i=0;i<knri*3+10;i++){
            outparms[i]=resultparmc[i];
            outparmserr[i]=resulterrmc[i];
        }
        treeout->Fill();
    }

    //!Perform the fit again
    //fitter.Config().SetParamsSettings(knri*3+10,parms);
    //fitter.FitFCN(knri*3+10,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),false);

    //! Filling the output tree
    for (unsigned int i=0;i<knri*3+10;i++){
      outparms[i]=resultpar[i];
      outparmserr[i]=resulterr[i];
    }
    treeout->Fill();
    */


    cout<<"\n*******PRINTING RESULT**********\n"<<endl;

    cout<<"\n*******PRINTING MAIN RESULTS**********\n"<<endl;
    cout<<"t1/2\terr\tp1n\terr\tp2n\terr"<<endl;



    cout<<"\n*****************\n"<<endl;
    cout<<log(2)/fB->GetParameter(0)<<"\t"<<log(2)/fB->GetParameter(0)/fB->GetParameter(0)*fB->GetParError(0)<<"\n";


    //! write to text file
    std::ofstream ofs("fitresults.txt", std::ofstream::out | std::ofstream::app);
    //ofs<<parmsfilename<<"\t"<<log(2)/resultpar[0]<<"\t"<<log(2)/resultpar[0]/resultpar[0]*resulterr[0]<<"\t"<<resultpar[knri]<<"\t"<<resulterr[knri]<<"\t"<<resultpar[knri*2]<<"\t"<<resulterr[knri*2]<<endl;


    TH1F* histcomphB=new TH1F("histcomphB","residual",hB->GetNbinsX(),hB->GetXaxis()->GetXmin(),hB->GetXaxis()->GetXmax());
    for (Int_t i=0;i<hB->GetNbinsX();i++){
        double res;
        if (hB->GetBinError(i+1)<0.001)
          res=  0;
        else
        res=  (hB->GetBinContent(i+1)- fB->Eval( hB->GetBinCenter(i+1) ) )/hB->GetBinError(i+1);
        histcomphB->SetBinContent(i+1,res);
        histcomphB->SetBinError(i+1,0);
    }
    //! Input descendant decay
    TF1* fbBparent=new TF1("fbBparent",fcn_decay_parent,plotrange[0],plotrange[1],knri*3+10);
    TF1* fbBdes=new TF1("fbBdes",fcn_decay_des,plotrange[0],plotrange[1],knri*3+10);
    for (Int_t i=0;i<knri*3+10;i++) {
        fbBparent->SetParameter(i,fB->GetParameter(i));
        fbBdes->SetParameter(i,fB->GetParameter(i));
    }
    auto c1 = new TCanvas("c1", "fit residual simple");
    gPad->SetFrameFillStyle(0);
    gStyle->SetOptStat(11111);
    fB->SetRange(plotrange[0], plotrange[1]);
    fB->SetLineColor(kRed);
    fB->SetNpx(binning);

    hB->SetTitle("Decay curve");
    hB->GetListOfFunctions()->Add(fB);
    hB->GetXaxis()->SetRangeUser(plotrange[2],plotrange[1]);
    hB->GetXaxis()->SetTitle("Time (s)");
    hB->GetYaxis()->SetTitle("Counts");
    hB->GetXaxis()->SetTitleSize(0.05);
    hB->GetYaxis()->SetTitleSize(0.05);
    hB->SetMarkerStyle(20);
    hB->SetMarkerSize(0.8);
    hB->GetXaxis()->SetLabelSize(0.05);
    hB->GetYaxis()->SetLabelSize(0.05);
    hB->Draw("P0 E");


    fbBparent->SetLineColor(kBlue);
    fbBparent->SetNpx(binning);
    fbBparent->Draw("same");

    fbBdes->SetLineColor(kGreen);
    fbBdes->SetNpx(binning);
    fbBdes->Draw("same");

    hdecay->GetFunction("pol1")->SetLineColor(kYellow);
    hdecay->GetFunction("pol1")->Draw("same");

    c1->Clear();

    /*
        auto rp1 = new TRatioPlot(hB, "",resul);
        rp1->SetGraphDrawOpt("L");
        rp1->SetSeparationMargin(0.0);
        rp1->Draw();
        rp1->GetLowerRefGraph()->SetMinimum(-2);
        rp1->GetLowerRefGraph()->SetMaximum(2);
        c1->Update();
*/

    hB->Write();
    fB->Write();
    fbBparent->Write();
    fbBdes->Write();
    histcomphB->Write();

    c1->Write();
    outfile->Close();
}


