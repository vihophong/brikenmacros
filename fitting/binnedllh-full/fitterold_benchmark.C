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
Double_t neueff=0.605;//changed to 62 %

Bool_t reject=false;
Double_t rejectrange=0.05;//first 50 ms

//! variable
Int_t npaths=100;
Int_t ndecay[kmaxpaths];
Int_t decaymap[kmaxpaths][kmaxndecay];
Int_t nneu[kmaxpaths][kmaxndecay];

Double_t parms[knri*2+8];
Double_t parmserr[knri*2+8];
Double_t parmsmax[knri*2+8];
Double_t parmsmin[knri*2+8];
Bool_t isparmsfix[knri*2+8];
Double_t mcparms[knri*2+8];

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

//! Function for decay with 2 neutron emission
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

//! Function for decay with 2 neutron emission
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

void fitter()
{
    Double_t lowerlimit=-10;
    Double_t upperlimit=10;
    Double_t nsigma=2.;

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


    getparms(parms,parmserr,parmsmax,parmsmin,isparmsfix,"testinputold.txt",nsigma);

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
    for (int i=0;i<knri*2+1;i++){
        fB->SetParameter(i,parms[i]);
        fB->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (parmserr[i]==0||isparmsfix[i]){
           cout<<"fixed value p"<<i<<"\tval="<<parms[i]<<endl;
           fB->FixParameter(i,parms[i]);
        }else{
           cout<<"variable value p"<<i<<"\tval "<<parms[i]<<"\tmin="<<parmsmin[i]<<"\tmax="<<parmsmax[i]<<endl;
        }
    }
    fB->SetParameter(knri*2+1,1000);
    fB->SetParLimits(knri*2+1,0,10000);
    fB->SetParameter(knri*2+2,2000);
    fB->SetParLimits(knri*2+2,0,10000);

    //!******************************************Delayed Neutron decay function
    //! read input and get parameters
    for (int i=0;i<knri*2+1;i++){
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
    fSB->SetParameter(knri*2+1,1000); //inital activity
    fSB->SetParLimits(knri*2+1,0,10000);
    fSB->SetParameter(knri*2+2,500); // background
    fSB->SetParLimits(knri*2+2,0,10000);
    fSB->SetParameter(knri*2+3,0.02);//random coincidence 1 neutron factor
    fSB->SetParLimits(knri*2+3,0.,1.);//random coincidence 1 neutron  factor
    fSB->SetParameter(knri*2+4,0.022);//random coincidence gt1 neutron factor
    fSB->SetParLimits(knri*2+4,0.,1.);//random coincidence gt1 neutron factor

    //!******************************************Delayed 2 Neutron decay function
    //! read input and get parameters
    for (int i=0;i<knri*2+1;i++){
        fSB2->SetParameter(i,parms[i]);
        fSB2->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (parmserr[i]==0||isparmsfix[i]){
           cout<<"fixed value p"<<i<<"\tval="<<parms[i]<<endl;
           fSB2->FixParameter(i,parms[i]);
        }else{
           cout<<"variable value p"<<i<<"\tval "<<parms[i]<<"\tmin="<<parmsmin[i]<<"\tmax="<<parmsmax[i]<<endl;
        }
    }

    fSB2->SetParameter(knri*2+1,1000); //inital activity
    fSB2->SetParLimits(knri*2+1,0,10000);
    fSB2->SetParameter(knri*2+2,100); // background
    fSB2->SetParLimits(knri*2+2,0,10000);
    fSB2->SetParameter(knri*2+3,0.02);//random coincidence 1 neutron factor
    fSB2->SetParLimits(knri*2+3,0.,1.);//random coincidence 1 neutron factor
    fSB2->SetParameter(knri*2+4,0.022);//random coincidence gt1 neutron factor
    fSB2->SetParLimits(knri*2+4,0.,1.);//random coincidence gt1 neutron factor
    fSB2->SetParameter(knri*2+5,0.005);//random coincidence 2 neutron factor
    fSB2->SetParLimits(knri*2+5,0.,1.);//random coincidence 2 neutron factor

    cout<<fB->Eval(1.)<<endl;
    cout<<fSB->Eval(1.)<<endl;
    cout<<fSB2->Eval(1.)<<endl;



    TCanvas* c1=new TCanvas("c1","c1",900,1200);
    c1->Divide(1,3);
    c1->cd(1);
    fB->Draw();
    c1->cd(2);
    fSB->Draw();
    c1->cd(3);
    fSB2->Draw();
}
