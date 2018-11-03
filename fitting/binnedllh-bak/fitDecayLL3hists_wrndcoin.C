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

void fitDecayLL3hists_wrndcoin(char* fitname,char* infile,char* parmsfile, Int_t ninterations, char* outfile, Int_t entrybegin, Int_t nentries, Int_t binning=2000){

    Double_t lowerlimit=-10;
    Double_t upperlimit=10;
    Double_t nsigma=2.;
    Double_t bkgactmaxmin=0.20; //*100% of max min bkg or initial activity

    Double_t plotrange[]={0.,10.};

    TRandom3* rseed=new TRandom3;

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

    //!****************************GET HISTOGRAM FROM FILE********************
    //!
    //!
    TFile *f = TFile::Open(infile);
    char tempchar1[1000];
    sprintf(tempchar1,"treeb");
    TTree* treeb=(TTree*) f->Get(tempchar1);
    cout<<entrybegin<<"\t EEE\t"<<nentries<<endl;
    if (nentries<0) nentries=treeb->GetEntries();
    treeb->Draw(Form("x>>hdecay(%d,%f,%f)",binning,-10.,10.),"","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>hdecaynuc1(%d,%f,%f)",binning,-10.,10.),"breal!=0&&(btype==4||btype==6||btype==1)","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>hdecaynuc2(%d,%f,%f)",binning,-10.,10.),"breal!=0&&(btype==2||btype==3)","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>hdecaynuc3(%d,%f,%f)",binning,-10.,10.),"breal!=0&&(btype==5||btype==10)","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>hdecaynuc4(%d,%f,%f)",binning,-10.,10.),"breal!=0&&(btype==7||btype==8)","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>hdecaynuc5(%d,%f,%f)",binning,-10.,10.),"breal!=0&&(btype==9||btype==14)","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>hdecaynuc6(%d,%f,%f)",binning,-10.,10.),"breal!=0&&(btype==11||btype==16)","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>hdecaynuc9(%d,%f,%f)",binning,-10.,10.),"breal!=0&&(btype==15||btype==17)","goff",nentries,entrybegin);

    treeb->Draw(Form("x>>hdecaybkg(%d,%f,%f)",binning,-10.,10.),"breal==0","goff",nentries,entrybegin);

    treeb->Draw(Form("x>>hdecay1n(%d,%f,%f)",binning,-10.,10.),"nfwd==1","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>hdecay2n(%d,%f,%f)",binning,-10.,10.),"nfwd==2","goff",nentries,entrybegin);

    treeb->Draw(Form("x>>hdecay1nbwd(%d,%f,%f)",binning,-10.,10.),"nbwd==1","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>hdecaygt0nbwd(%d,%f,%f)",binning,-10.,10.),"nbwd>0","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>hdecay2nbwd(%d,%f,%f)",binning,-10.,10.),"nbwd==2","goff",nentries,entrybegin);

    treeb->Draw(Form("x>>h1ncomp1(%d,%f,%f)",binning,-10.,10.),"breal==1&&nrealflag==0&&nfwd==1","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>h1ncomp2(%d,%f,%f)",binning,-10.,10.),"breal==1&&nrealflag==1&&btype!=6&&nfwd==1","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>h1ncomp3(%d,%f,%f)",binning,-10.,10.),"breal==1&&nrealflag==1&&btype==6&&nfwd==1","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>h1nbkg(%d,%f,%f)",binning,-10.,10.),"breal==0&&nfwd==1","goff",nentries,entrybegin);


    treeb->Draw(Form("x>>h2ncomp1(%d,%f,%f)",binning,-10.,10.),"breal==1&&nrealflag==0&&nfwd==2","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>h2ncomp2(%d,%f,%f)",binning,-10.,10.),"breal==1&&nrealflag==2&&btype==6&&nfwd==2","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>h2ncomp3(%d,%f,%f)",binning,-10.,10.),"breal==1&&nrealflag==1&&btype!=6&&nfwd==2","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>h2ncomp4(%d,%f,%f)",binning,-10.,10.),"breal==1&&nrealflag==1&&btype==6&&nfwd==2","goff",nentries,entrybegin);
    treeb->Draw(Form("x>>h2nbkg(%d,%f,%f)",binning,-10.,10.),"breal==0&&nfwd==2","goff",nentries,entrybegin);


    sprintf(tempchar1,"hdecay");
    TH1F* hdecay=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecaynuc1");
    TH1F* hdecaynuc1=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecaynuc2");
    TH1F* hdecaynuc2=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecaynuc3");
    TH1F* hdecaynuc3=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecaynuc4");
    TH1F* hdecaynuc4=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecaynuc5");
    TH1F* hdecaynuc5=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecaynuc6");
    TH1F* hdecaynuc6=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecaynuc9");
    TH1F* hdecaynuc9=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecaybkg");
    TH1F* hdecaybkg=(TH1F*) gDirectory->Get(tempchar1);


    sprintf(tempchar1,"hdecay1n");
    TH1F* hdecay1n=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecay2n");
    TH1F* hdecay2n=(TH1F*) gDirectory->Get(tempchar1);

    sprintf(tempchar1,"h1ncomp1");
    TH1F* h1ncomp1=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"h1ncomp2");
    TH1F* h1ncomp2=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"h1ncomp3");
    TH1F* h1ncomp3=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"h1nbkg");
    TH1F* h1nbkg=(TH1F*) gDirectory->Get(tempchar1);


    sprintf(tempchar1,"h2ncomp1");
    TH1F* h2ncomp1=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"h2ncomp2");
    TH1F* h2ncomp2=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"h2ncomp3");
    TH1F* h2ncomp3=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"h2ncomp4");
    TH1F* h2ncomp4=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"h2nbkg");
    TH1F* h2nbkg=(TH1F*) gDirectory->Get(tempchar1);

    sprintf(tempchar1,"hdecay1nbwd");
    TH1F* hdecay1nbwd=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecaygt0nbwd");
    TH1F* hdecaygt0nbwd=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecay2nbwd");
    TH1F* hdecay2nbwd=(TH1F*) gDirectory->Get(tempchar1);

    Double_t n1nbwd=(Double_t) hdecay1nbwd->GetEntries();
    Double_t gt0nbwd=(Double_t) hdecaygt0nbwd->GetEntries();
    Double_t n2nbwd=(Double_t) hdecay2nbwd->GetEntries();
    Double_t nball=(Double_t) hdecay->GetEntries();

    parms[knri*2+5]=n1nbwd/nball;
    parms[knri*2+6]=gt0nbwd/nball;
    parms[knri*2+7]=n2nbwd/nball;

    
    //! background parameter estimation
    //! background as average of several first bin
    Int_t bincnt=0;
    Double_t bkgsum=0;
    for (Int_t i=hdecay->GetXaxis()->FindBin(-5);i<hdecay->GetXaxis()->FindBin(-0.5);i++){
        bkgsum+=hdecay->GetBinContent(i);
        bincnt++;
    }
    parms[knri*2+2]=bkgsum/((Double_t)bincnt);
    parmserr[knri*2+2]=parms[knri*2+2]*bkgactmaxmin;
    parmsmin[knri*2+2]=parms[knri*2+2]-parms[knri*2+2]*bkgactmaxmin;
    parmsmax[knri*2+2]=parms[knri*2+2]+parms[knri*2+2]*bkgactmaxmin;

    //activity
    parms[knri*2+1]=hdecay->GetBinContent(hdecay->GetXaxis()->FindBin(rejectrange))-parms[knri*2+2];
    parmsmin[knri*2+1]=parms[knri*2+1]-parms[knri*2+1]*bkgactmaxmin;
    parmsmax[knri*2+1]=parms[knri*2+1]*2+parms[knri*2+1]*bkgactmaxmin;

    bincnt=0;
    bkgsum=0;
    for (Int_t i=hdecay1n->GetXaxis()->FindBin(-5);i<hdecay1n->GetXaxis()->FindBin(-0.5);i++){
        bkgsum+=hdecay1n->GetBinContent(i);
        bincnt++;
    }
    parms[knri*2+3]=bkgsum/((Double_t)bincnt);
    parmsmin[knri*2+3]=parms[knri*2+3]-parms[knri*2+3]*bkgactmaxmin;
    parmsmax[knri*2+3]=parms[knri*2+3]+parms[knri*2+3]*bkgactmaxmin;

    bincnt=0;
    bkgsum=0;
    for (Int_t i=hdecay2n->GetXaxis()->FindBin(-5);i<hdecay2n->GetXaxis()->FindBin(-0.5);i++){
        bkgsum+=hdecay2n->GetBinContent(i);
        bincnt++;
    }
    parms[knri*2+4]=bkgsum/((Double_t)bincnt);
    parmsmin[knri*2+4]=parms[knri*2+4]-parms[knri*2+4]*bkgactmaxmin;
    parmsmax[knri*2+4]=parms[knri*2+4]+parms[knri*2+4]*bkgactmaxmin;
    

    //! error of this parameters:
    //parmserr[knri*2+5]=n1nbwd/nball*sqrt(1/n1nbwd/n1nbwd+1/nball/nball);
    //parmserr[knri*2+6]=gt0nbwd/nball*sqrt(1/gt0nbwd/gt0nbwd+1/nball/nball);
    //parmserr[knri*2+7]=n2nbwd/nball*sqrt(1/n2nbwd/n2nbwd+1/nball/nball);

    cout<<"count = "<<n1nbwd<<"\t"<<gt0nbwd<<"\t"<<n2nbwd<<"\t"<<nball<<endl;
    cout<<"parameter error = "<<parmserr[knri*2+5]<<"\t"<<parmserr[knri*2+6]<<"\t"<<parmserr[knri*2+7]<<endl;


    //!******************************************Define All BETA decay function
    TF1* fB=new TF1("fB",fcn_gen,lowerlimit,upperlimit,21);
    fB->SetNpx(2000);
    fB->SetLineWidth(2);
    fB->SetLineColor(8);
    //!******************************************Define Delayed Neutron decay function
    TF1* fSB=new TF1("fSB",fcn_gen_wneutron,lowerlimit,upperlimit,23);
    fSB->SetNpx(2000);
    fSB->SetLineWidth(2);
    fSB->SetLineColor(8);
    //!******************************************Define Delayed 2 Neutron decay function
    TF1* fSB2=new TF1("fSB2",fcn_gen_w2neutron,lowerlimit,upperlimit,24);
    fSB2->SetNpx(2000);
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
    fSB->SetParameter(knri*2+4,0.5);//random coincidence 2 neutron factor
    fSB->SetParLimits(knri*2+4,0.,1.);//random coincidence 2 neutron factor

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

    fSB2->SetParameter(knri*2+4,0.5);//random coincidence 2 neutron factor
    fSB2->SetParLimits(knri*2+4,0.,1.);//random coincidence 2 neutron factor

    fSB2->SetParameter(knri*2+5,0.5);//random coincidence 3 neutron factor
    fSB2->SetParLimits(knri*2+5,0.,1.);//random coincidence 3 neutron factor

    cout<<fSB2->Eval(1.)<<endl;


   TH1F * hB = (TH1F*) hdecay->Clone();
   TH1F * hSB = (TH1F*) hdecay1n->Clone();
   TH1F * hSB2 = (TH1F*) hdecay2n->Clone();


   //!****************************GET HISTOGRAM FROM FILE********************
   //! book file and tree
   TFile* fout=new TFile(outfile,"recreate");
   Double_t outparms[knri*2+8];
   Double_t outparmserr[knri*2+8];
   Int_t isvary[knri*2+6];

   for (int i=0;i<knri*2+6;i++){
       outparms[i]=0;
       outparmserr[i]=0;
       isvary[i]=0;
       if (isparmsfix[i]&&(parmserr[i]>0)){
           if(i<knri*2) isvary[i]=1;
       }
   }

   sprintf(tempchar1,"tree%s",fitname);
   TTree* treeout=new TTree(tempchar1,tempchar1);
   treeout->Branch("outparms",outparms,Form("outparms[%d]/D",knri*2+8));
   treeout->Branch("outparmserr",outparmserr,Form("outparmserr[%d]/D",knri*2+8));
   treeout->Branch("isvary",outparmserr,Form("isvary[%d]/I",knri*2+8));


   ROOT::Math::WrappedMultiTF1 wfB(*fB,1);
   ROOT::Math::WrappedMultiTF1 wfSB(*fSB,1);
   ROOT::Math::WrappedMultiTF1 wfSB2(*fSB2,1);

   ROOT::Fit::DataOptions opt;

   //! limit within the fitting range
   opt.fUseRange  =true;
   // set the data range

   ROOT::Fit::DataRange rangeB;
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

   ROOT::Fit::PoissonLLFunction chi2_B(dataB, wfB);
   ROOT::Fit::PoissonLLFunction chi2_SB(dataSB, wfSB);
   ROOT::Fit::PoissonLLFunction chi2_SB2(dataSB2, wfSB2);
   GlobalChi2 globalChi2(chi2_B, chi2_SB, chi2_SB2);


   ROOT::Fit::Fitter fitter;

   // create before the parameter settings in order to fix or set range on them

   cout<<"\n\n*****SETTING PARAMETERS......\n"<<endl;
   fitter.Config().SetParamsSettings(knri*2+8,parms);
   for (int i=0;i<knri*2+8;i++){
       if (parmserr[i]==0||isparmsfix[i]){
          cout<<"fixed value p"<<i<<"\tval="<<parms[i]<<endl;
          fitter.Config().ParSettings(i).Fix();
       }else{
          fitter.Config().ParSettings(i).SetLimits(parmsmin[i],parmsmax[i]);
          cout<<"variable value p"<<i<<"\tval "<<parms[i]<<"\tmin="<<parmsmin[i]<<"\tmax="<<parmsmax[i]<<endl;
       }
   }

   fitter.Config().MinimizerOptions().SetPrintLevel(1);
   //fitter.Config().SetMinimizer("Minuit2","Migrad");
   //fitter.Config().SetMinosErrors();


   if (fitter.Config().MinosErrors()) cout<<"minos enabled"<<endl;

   cout<<"\n\n************FITTING.........."<<endl;
   //! fit FCN function directly
   // (specify optionally data size and flag to indicate that is a chi2 fit)
   for (int i=0;i<knri*2+8;i++) mcparms[i]=parms[i];

   for(int i=0;i<ninterations;i++) {
       cout<<"mc "<<i+1<<endl;
       mc(rseed);       
       fitter.Config().SetParamsSettings(knri*2+8,mcparms);
       for (unsigned int j=0;j<knri*2+8;i++){
           if (parmserr[j]==0||isparmsfix[j]){
              fitter.Config().ParSettings(i).Fix();
           }else{
              fitter.Config().ParSettings(i).SetLimits(parmsmin[j],parmsmax[j]);
           }
       }
       fitter.FitFCN(knri*2+8,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),false);
       const Double_t* resultpar=fitter.Result().GetParams();
       const Double_t* resulterr=fitter.Result().GetErrors();
       for (unsigned int i=0;i<knri*2+8;i++){
           outparms[i]=resultpar[i];
           outparmserr[i]=resulterr[i];
       }
       treeout->Fill();
   }


   fitter.Config().SetParamsSettings(knri*2+8,parms);
   fitter.FitFCN(knri*2+8,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),false);

   ROOT::Fit::FitResult result = fitter.Result();
   cout<<"\n*******PRINTING RESULT**********\n"<<endl;
   //result.Print(std::cout);
   const Double_t* resultpar=result.GetParams();
   const Double_t* resulterr=result.GetErrors();
   for (unsigned int i=0;i<result.NPar();i++){
       cout<<"p"<<i<<" = "<<resultpar[i]<<" +/- "<<resulterr[i]<<" le= "<<result.LowerError(i)<<" ue= "<<result.UpperError(i)<<endl;
   }

   cout<<"\n*****************\n"<<endl;
   Double_t nimplants=1;

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
   ofs<<resultpar[18]<<","<<resulterr[18]<<","<<resultpar[19]/hB->GetBinWidth(1)/resultpar[0]/nimplants*100<<","<<binning<<","<<lowerlimit<<","<<upperlimit<<endl;


   TH1F* histcomphBcomb=new TH1F("hist_of_residual_decay","hist of residual decay",50,-10,10);
   TH1F* histcomphB=new TH1F("residual_decay","residual decay",hB->GetNbinsX(),hB->GetXaxis()->GetXmin(),hB->GetXaxis()->GetXmax());
   for (Int_t i=0;i<hB->GetNbinsX();i++){
       double res;
       //if (hB->GetBinError(i+1)<0.001)
       //res=  0;
       //else
       res=  (hB->GetBinContent(i+1)- fB->Eval( hB->GetBinCenter(i+1) ) )/hB->GetBinError(i+1);

       histcomphB->SetBinContent(i+1,res);
       histcomphB->SetBinError(i+1,0);
       if (hB->GetBinCenter(i+1)<plotrange[1]&&hB->GetBinCenter(i+1)>plotrange[0])
       histcomphBcomb->Fill(res);
   }
   TH1F* histcomphSBcomb=new TH1F("hist_of_residual_decay1neu","hist of residual decay 1neu",50,-10,10);
   TH1F* histcomphSB=new TH1F("residual_decay1neu","residual decay 1 neutron",hSB->GetNbinsX(),hSB->GetXaxis()->GetXmin(),hSB->GetXaxis()->GetXmax());
   for (Int_t i=0;i<hSB->GetNbinsX();i++){
       double res;
       //if (hSB->GetBinError(i+1)<0.001)
       // res =  0;
       //else
        res =  (hSB->GetBinContent(i+1)- fSB->Eval( hSB->GetBinCenter(i+1) ) )/hSB->GetBinError(i+1);
       histcomphSB->SetBinContent(i+1,res);
       histcomphSB->SetBinError(i+1,0);
       if (hSB->GetBinCenter(i+1)<plotrange[1]&&hSB->GetBinCenter(i+1)>plotrange[0])
       histcomphSBcomb->Fill(res);
   }

   TH1F* histcomphSB2comb=new TH1F("hist_of_residual_decay2neu","hist of residual decay 2neu",50,-10,10);
   TH1F* histcomphSB2=new TH1F("residual_decay2neu","residual decay 2 neutrons",hSB2->GetNbinsX(),hSB2->GetXaxis()->GetXmin(),hSB2->GetXaxis()->GetXmax());
   for (Int_t i=0;i<hSB2->GetNbinsX();i++){
       double res;
       //if (hSB2->GetBinError(i+1)<0.001)
        //res =  0;
       //else
        res =  (hSB2->GetBinContent(i+1)- fSB2->Eval( hSB2->GetBinCenter(i+1) ) )/hSB2->GetBinError(i+1);
       histcomphSB2->SetBinContent(i+1,res);
       histcomphSB2->SetBinError(i+1,0);
       if (hSB2->GetBinCenter(i+1)<plotrange[1]&&hSB2->GetBinCenter(i+1)>plotrange[0])
       histcomphSB2comb->Fill(res);
   }


   sprintf(tempchar1,"Simfit%s",fitname);
   TCanvas * c1 = new TCanvas(tempchar1,"Simultaneous fit of 3 histograms",
                              10,10,800,600);

   c1->Divide(1,2,0,0);
   c1->cd(1);
   gPad->SetBottomMargin(0.001);
   gPad->SetTopMargin(0.01);
   gPad->SetRightMargin(0.01);

   gStyle->SetOptFit(1111);
   fB->SetFitResult( result, iparB);
   fB->SetRange(rangeB().first, rangeB().second);
   fB->SetLineColor(kBlue);
   fB->SetNpx(binning);
   hB->GetListOfFunctions()->Add(fB);
   hB->GetXaxis()->SetRangeUser(plotrange[0],plotrange[1]);
   hB->SetMarkerStyle(20);
   hB->SetMarkerSize(1);
   hB->Draw("P0");

   c1->cd(2);
   gPad->SetBottomMargin(0.1);
   gPad->SetTopMargin(0.001);
   gPad->SetRightMargin(0.01);
   histcomphB->GetXaxis()->SetRangeUser(plotrange[0],plotrange[1]);
   histcomphB->Draw("hist");

   sprintf(tempchar1,"Simfit1%s",fitname);
   TCanvas * c2 = new TCanvas(tempchar1,"Simultaneous fit of 3 histograms",
                              10,10,800,600);

   c2->Divide(1,2,0,0);
   c2->cd(1);
   gPad->SetBottomMargin(0.001);
   gPad->SetTopMargin(0.01);
   gPad->SetRightMargin(0.01);

   fSB->SetFitResult( result, iparSB);
   fSB->SetRange(rangeSB().first, rangeSB().second);
   fSB->SetLineColor(kRed);
   fSB->SetNpx(binning);
   hSB->GetListOfFunctions()->Add(fSB);
   hSB->GetXaxis()->SetRangeUser(plotrange[0],plotrange[1]);
   hSB->SetMarkerStyle(20);
   hSB->SetMarkerSize(1);
   hSB->Draw("P0");

   c2->cd(2);
   gPad->SetBottomMargin(0.1);
   gPad->SetTopMargin(0.001);
   gPad->SetRightMargin(0.01);
   histcomphSB->GetXaxis()->SetRangeUser(plotrange[0],plotrange[1]);
   histcomphSB->Draw("hist");


   sprintf(tempchar1,"Simfit2%s",fitname);
   TCanvas * c3 = new TCanvas(tempchar1,"Simultaneous fit of 3 histograms",
                              10,10,800,600);

   c3->Divide(1,2,0,0);
   c3->cd(1);
   gPad->SetBottomMargin(0.001);
   gPad->SetTopMargin(0.01);
   gPad->SetRightMargin(0.01);

   fSB2->SetFitResult( result, iparSB2);
   fSB2->SetRange(rangeSB2().first, rangeSB2().second);
   fSB2->SetLineColor(kGreen);
   fSB2->SetNpx(binning);
   hSB2->GetListOfFunctions()->Add(fSB2);
   hSB2->Draw();
   hSB2->GetXaxis()->SetRangeUser(plotrange[0],plotrange[1]);
   hSB2->SetMarkerStyle(20);
   hSB2->SetMarkerSize(1);
   hSB2->Draw("P0");

   c3->cd(2);
   gPad->SetBottomMargin(0.1);
   gPad->SetTopMargin(0.001);
   gPad->SetRightMargin(0.01);
   //histcomphSB2->Rebin(5);
   histcomphSB2->GetXaxis()->SetRangeUser(plotrange[0],plotrange[1]);
   histcomphSB2->Draw("hist");


   sprintf(tempchar1,"Simfit3%s",fitname);
   TCanvas * c4 = new TCanvas(tempchar1,"Combined residuals histograms",
                              10,10,800,1200);


   c4->Divide(1,3);
   c4->cd(1);
   histcomphBcomb->Draw();
   c4->cd(2);
   histcomphSBcomb->Draw();
   c4->cd(3);
   histcomphSB2comb->Draw();

   fout->cd();
   c1->Write();
   c2->Write();
   c3->Write();
   c4->Write();
   histcomphBcomb->Write();
   histcomphB->Write();
   histcomphSBcomb->Write();
   histcomphSB->Write();
   histcomphSB2comb->Write();
   histcomphSB2->Write();

   hdecaynuc1->Write();
   hdecaynuc2->Write();
   hdecaynuc3->Write();
   hdecaynuc4->Write();
   hdecaynuc5->Write();
   hdecaynuc6->Write();
   hdecaynuc9->Write();
   hdecaybkg->Write();

   h1ncomp1->Write();
   h1ncomp2->Write();
   h1ncomp3->Write();
   h1nbkg->Write();

   h2ncomp1->Write();
   h2ncomp2->Write();
   h2ncomp3->Write();
   h2ncomp4->Write();
   h2nbkg->Write();

   TF1* fbBnuc1=new TF1("fbBnuc1",fcn_gennuc1,lowerlimit,upperlimit,knri*2+3);
   TF1* fbBnuc2=new TF1("fbBnuc2",fcn_gennuc2,lowerlimit,upperlimit,knri*2+3);
   TF1* fbBnuc3=new TF1("fbBnuc3",fcn_gennuc3,lowerlimit,upperlimit,knri*2+3);
   TF1* fbBnuc4=new TF1("fbBnuc4",fcn_gennuc4,lowerlimit,upperlimit,knri*2+3);
   TF1* fbBnuc5=new TF1("fbBnuc5",fcn_gennuc5,lowerlimit,upperlimit,knri*2+3);
   TF1* fbBnuc6=new TF1("fbBnuc6",fcn_gennuc6,lowerlimit,upperlimit,knri*2+3);
   TF1* fbBnuc9=new TF1("fbBnuc9",fcn_gennuc9,lowerlimit,upperlimit,knri*2+3);
   TF1* fBcbkg=new TF1("fBcbkg",fcn_gen_bkg,lowerlimit,upperlimit,knri*2+3);

   TF1* fSBc1=new TF1("fSBc1",fcn_gen_wneutronc1,lowerlimit,upperlimit,knri*2+5);
   TF1* fSBc2=new TF1("fSBc2",fcn_gen_wneutronc2,lowerlimit,upperlimit,knri*2+5);
   TF1* fSBc3=new TF1("fSBc3",fcn_gen_wneutronc3,lowerlimit,upperlimit,knri*2+5);
   TF1* fSBcbkg=new TF1("fSBcbkg",fcn_gen_wneutron_bkg,lowerlimit,upperlimit,knri*2+5);

   TF1* fSB2c1=new TF1("fSB2c1",fcn_gen_w2neutronc1,lowerlimit,upperlimit,knri*2+6);
   TF1* fSB2c2=new TF1("fSB2c2",fcn_gen_w2neutronc2,lowerlimit,upperlimit,knri*2+6);
   TF1* fSB2c3=new TF1("fSB2c3",fcn_gen_w2neutronc3,lowerlimit,upperlimit,knri*2+6);
   TF1* fSB2c4=new TF1("fSB2c4",fcn_gen_w2neutronc4,lowerlimit,upperlimit,knri*2+6);
   TF1* fSB2cbkg=new TF1("fSB2cbkg",fcn_gen_w2neutron_bkg,lowerlimit,upperlimit,knri*2+5);




   fbBnuc1->SetNpx(binning);
   fbBnuc2->SetNpx(binning);
   fbBnuc3->SetNpx(binning);
   fbBnuc4->SetNpx(binning);
   fbBnuc5->SetNpx(binning);
   fbBnuc6->SetNpx(binning);
   fbBnuc9->SetNpx(binning);

   fSBc1->SetNpx(binning);
   fSBc2->SetNpx(binning);
   fSBc3->SetNpx(binning);
   fSBcbkg->SetNpx(binning);

   fSB2c1->SetNpx(binning);
   fSB2c2->SetNpx(binning);
   fSB2c3->SetNpx(binning);
   fSB2c4->SetNpx(binning);

   for (Int_t i=0;i<knri*2+5;i++){
       if (i<knri*2+3){
           fbBnuc1->FixParameter(i,fB->GetParameter(i));
           fbBnuc2->FixParameter(i,fB->GetParameter(i));
           fbBnuc3->FixParameter(i,fB->GetParameter(i));
           fbBnuc4->FixParameter(i,fB->GetParameter(i));
           fbBnuc5->FixParameter(i,fB->GetParameter(i));
           fbBnuc6->FixParameter(i,fB->GetParameter(i));
           fbBnuc9->FixParameter(i,fB->GetParameter(i));
       }

       fSBc1->FixParameter(i,fSB->GetParameter(i));
       fSBc2->FixParameter(i,fSB->GetParameter(i));
       fSBc3->FixParameter(i,fSB->GetParameter(i));
       fSBcbkg->FixParameter(i,fSB->GetParameter(i));

       fSB2c1->FixParameter(i,fSB2->GetParameter(i));
       fSB2c2->FixParameter(i,fSB2->GetParameter(i));
       fSB2c3->FixParameter(i,fSB2->GetParameter(i));
       fSB2c4->FixParameter(i,fSB2->GetParameter(i));
       fSB2cbkg->FixParameter(i,fSB2->GetParameter(i));
   }
   fSB2c1->FixParameter(knri*2+5,fSB2->GetParameter(knri*2+5));
   fSB2c2->FixParameter(knri*2+5,fSB2->GetParameter(knri*2+5));
   fSB2c3->FixParameter(knri*2+5,fSB2->GetParameter(knri*2+5));
   fSB2c4->FixParameter(knri*2+5,fSB2->GetParameter(knri*2+5));
   fSB2cbkg->FixParameter(knri*2+5,fSB2->GetParameter(knri*2+5));

   fB->Write();
   fbBnuc1->Write();
   fbBnuc2->Write();
   fbBnuc3->Write();
   fbBnuc4->Write();
   fbBnuc5->Write();
   fbBnuc6->Write();
   fbBnuc9->Write();
   fBcbkg->Write();


   fSB->Write();
   fSBc1->Write();
   fSBc2->Write();
   fSBc3->Write();
   fSBcbkg->Write();

   fSB2c1->Write();
   fSB2c2->Write();
   fSB2c3->Write();
   fSB2c4->Write();
   fSB2cbkg->Write();


   fSB2->Write();

   hB->Write();
   hSB->Write();
   hSB2->Write();
   treeout->Write();
   fout->Close();
}

