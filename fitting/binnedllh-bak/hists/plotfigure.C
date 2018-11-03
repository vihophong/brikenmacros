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



const Int_t knri=9;
const Int_t kmaxpar=5;
const Int_t kmaxndecay=10;
const Int_t kmaxpaths=100;
//Double_t neueff=0.66*(100-0.8)/100;
Double_t neueff=0.62;//changed to 62 %

Bool_t reject=false;
Double_t rejectrange=0.05;//first 50 ms

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

    //! read initial activity, background and random coincidence factors
    for (int i=0;i<7;i++){
        ifs>>parms[knri*2+i+1]>>parmsmin[knri*2+i+1]>>parmsmax[knri*2+i+1];
        parmserr[knri*2+i+1]=parms[knri*2+i+1]-parmsmin[knri*2+i+1];
        cout<<knri*2+i+1<<"\t"<<parms[knri*2+i+1]<<"\t"<<parmsmin[knri*2+i+1]<<"\t"<<parmsmax[knri*2+i+1]<<endl;
        isparmsfix[knri*2+i+1]=false;
        if (i>3) {
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

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}


//! Global function
Double_t fcn_gen_bkg(Double_t *x, Double_t *par) {
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;
    return returnval;
}



Double_t fcn_gennuc1(Double_t *x, Double_t *par) {
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];

    //! Parent nuclei
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0]);

    //for (Int_t i=0;i<npaths;i++){
        //returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0]);
    //}

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

Double_t fcn_gennuc2(Double_t *x, Double_t *par) {
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];


    for (Int_t i=0;i<npaths;i++){
        if (i==0)
        returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

Double_t fcn_gennuc3(Double_t *x, Double_t *par) {
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];


    for (Int_t i=0;i<npaths;i++){
        if (i==3)
        returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

Double_t fcn_gennuc6(Double_t *x, Double_t *par) {
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];


    for (Int_t i=0;i<npaths;i++){
        if (i==4)
        returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}


Double_t fcn_gennuc4(Double_t *x, Double_t *par) {
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];


    for (Int_t i=0;i<npaths;i++){
        if (i==1)
        returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

Double_t fcn_gennuc5(Double_t *x, Double_t *par) {
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];


    for (Int_t i=0;i<npaths;i++){
        if (i==5||i==6)
        returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

Double_t fcn_gennuc9(Double_t *x, Double_t *par) {
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];

    for (Int_t i=0;i<npaths;i++){
        if (i==7||i==8)
        returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}




//! Global function
Double_t fcn_gen_wneutron(Double_t *x, Double_t *par) {
    Double_t randcoinfgt0n=par[knri*2+4];
    Double_t randcoinf1n=par[knri*2+3];
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];

    //! Parent nuclei
    //! decay with 1 neutron of parent
    returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf1n-randcoinfgt0n);
    //! decay with 1 neutron of parent
    returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf1n-randcoinfgt0n);

    //! decay with 2 neutron of parent (not random 1 neutron)
    returnval-=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;

    //! random coinc part
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;

    for (Int_t i=0;i<npaths;i++){
        //! random coinc of beta decay of daugter
        returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf1n;
        //! decay with 1 neutron part of daugter nuclei
        returnval+=neueff*pn[decaymap[i][ndecay[i]-1]-1]*lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(1-randcoinf1n-randcoinfgt0n);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

Double_t fcn_gen_wneutronc2(Double_t *x, Double_t *par) {
    Double_t randcoinfgt0n=par[knri*2+4];
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];

    //! Parent nuclei
    //! decay with 1 neutron of parent
    returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinfgt0n);

    for (Int_t i=0;i<npaths;i++){
        //! decay with 1 neutron part of daugter nuclei
        returnval+=neueff*pn[decaymap[i][ndecay[i]-1]-1]*lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(1-randcoinfgt0n);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

Double_t fcn_gen_wneutronc3(Double_t *x, Double_t *par) {
    Double_t randcoinfgt0n=par[knri*2+4];
    Double_t randcoinf1n=par[knri*2+3];
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];

    //! decay with 1 neutron of parent
    returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinfgt0n);

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

Double_t fcn_gen_wneutronc1(Double_t *x, Double_t *par) {
    Double_t randcoinfgt0n=par[knri*2+4];
    Double_t randcoinf1n=par[knri*2+3];
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];

    //! Parent nuclei
    //! decay with 1 neutron of parent
    returnval-=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;
    //! decay with 1 neutron of parent
    returnval-=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;

    //! decay with 2 neutron of parent (not random 1 neutron)
    returnval-=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;

    //! random coinc part
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;

    for (Int_t i=0;i<npaths;i++){
        //! random coinc of beta decay of daugter
        returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf1n;
        //! decay with 1 neutron part of daugter nuclei
        returnval-=neueff*pn[decaymap[i][ndecay[i]-1]-1]*lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf1n;
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}
Double_t fcn_gen_wneutron_bkg(Double_t *x, Double_t *par) {
    Double_t randcoinfgt0n=par[knri*2+4];
    Double_t randcoinf1n=par[knri*2+3];
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;
    return returnval;
}



//! Global function
Double_t fcn_gen_w2neutron(Double_t *x, Double_t *par) {
    Double_t randcoinf2n=par[knri*2+5];
    Double_t randcoinfgt0n=par[knri*2+4];
    Double_t randcoinf1n=par[knri*2+3];
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];

    //! parent
    //! decay with 2 neutron from P2n of parent
    returnval+=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf2n-randcoinfgt0n);

    //! random 1n decay of parent
    returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);
    //! decay with 1 neutron from P2n of parent
    returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

    //! random coinc part
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf2n;

    for (Int_t i=0;i<npaths;i++){
        //! random coinc of beta decay of daugter
        returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf2n;
        //! random 1n decay of daugter
        returnval+=neueff*pn[decaymap[i][ndecay[i]-1]-1]*lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

//! component 1 of delayed two neutron
Double_t fcn_gen_w2neutronc1(Double_t *x, Double_t *par) {
    Double_t randcoinf2n=par[knri*2+5];
    Double_t randcoinfgt0n=par[knri*2+4];
    Double_t randcoinf1n=par[knri*2+3];
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];

    //! parent
    //! decay with 2 neutron from P2n of parent
    returnval-=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf2n;

    //! random 1n decay of parent
    returnval-=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf2n;
    //! decay with 1 neutron from P2n of parent
    returnval-=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf2n;

    //! random coinc part
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf2n;

    for (Int_t i=0;i<npaths;i++){
        //! random coinc of beta decay of daugter
        returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf2n;
        //! random 1n decay of daugter
        returnval-=neueff*pn[decaymap[i][ndecay[i]-1]-1]*lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf2n;
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

//! component 2 of delayed two neutron
Double_t fcn_gen_w2neutronc2(Double_t *x, Double_t *par) {
    Double_t randcoinf2n=par[knri*2+5];
    Double_t randcoinfgt0n=par[knri*2+4];
    Double_t randcoinf1n=par[knri*2+3];
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];

    //! parent
    //! decay with 2 neutron from P2n of parent
    returnval+=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinfgt0n);

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}


//! component 3 of delayed 2 neutron
Double_t fcn_gen_w2neutronc3(Double_t *x, Double_t *par) {
    Double_t randcoinf2n=par[knri*2+5];
    Double_t randcoinfgt0n=par[knri*2+4];
    Double_t randcoinf1n=par[knri*2+3];
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];

    //! parent

    //! random 1n decay of parent
    returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(randcoinf1n*(1-randcoinfgt0n));

    for (Int_t i=0;i<npaths;i++){
        //! random 1n decay of daugter
        returnval+=neueff*pn[decaymap[i][ndecay[i]-1]-1]*lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n));
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

//! Global function
Double_t fcn_gen_w2neutronc4(Double_t *x, Double_t *par) {
    Double_t randcoinf2n=par[knri*2+5];
    Double_t randcoinfgt0n=par[knri*2+4];
    Double_t randcoinf1n=par[knri*2+3];
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knri*2+1]/par[0];


    //! decay with 1 neutron from P2n of parent
    returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(randcoinf1n*(1-randcoinfgt0n));

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

Double_t fcn_gen_w2neutron_bkg(Double_t *x, Double_t *par) {
    Double_t randcoinf2n=par[knri*2+5];
    Double_t randcoinfgt0n=par[knri*2+4];
    Double_t randcoinf1n=par[knri*2+3];
    Double_t bkg=par[knri*2+2];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;
    return returnval;
}


// definition of shared parameter
// background function
int iparB[21] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};

// signal + background function
int iparSB[23] ={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,23,24};

// signal 2n + background function
int iparSB2[24] ={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,22, 23,24,25};

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

      double p2[23];
      for (int i = 0; i < 23; ++i) p2[i] = par[iparSB[i] ];

      double p3[24];
      for (int i = 0; i < 24; ++i) p3[i] = par[iparSB2[i] ];

      return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3);
   }
   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
   const  ROOT::Math::IMultiGenFunction * fChi2_3;
};


Double_t parms[knri*2+8];
Double_t parmserr[knri*2+8];
Double_t parmsmax[knri*2+8];
Double_t parmsmin[knri*2+8];
Bool_t isparmsfix[knri*2+8];
Double_t mcparms[knri*2+8];

void mc(TRandom3* rseed)
{
    for (int i=0;i<knri*2;i++){
        if (isparmsfix[i]&&(parmserr[i]>0)){
            mcparms[i]=rseed->Gaus(parms[i],parmserr[i]);
        }
    }
    //! mc for neutron efficiency?
}

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

