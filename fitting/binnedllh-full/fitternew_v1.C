#include "makepath.C"


//Double_t neueff=0.66*(100-0.8)/100;
Double_t neueff=0.62;//changed to 62 %
Bool_t reject=false;
Double_t rejectrange=0.05;//first 50 ms

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


    Double_t bkg=par[knri*3+1];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;
    Double_t N0=par[knri*3]/par[0];

    //! Parent nuclei
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0]);

    for (Int_t i=0;i<npaths;i++){
        returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

//! Function for decay with 1 neutron emission
Double_t fcn_1ndecay(Double_t *x, Double_t *par) {

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t* p2n=&par[knri*2];

    Double_t randcoinfgt0n=par[knri*3+3];
    Double_t randcoinf1n=par[knri*3+2];
    Double_t bkg=par[knri*3+1];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;
    Double_t N0=par[knri*3]/par[0];

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
        returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf1n;
        //! decay with 1 neutron part of daugter nuclei
        returnval+=neueff*pn[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(1-randcoinf1n-randcoinfgt0n);
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
    Double_t bkg=par[knri*3+1];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;
    Double_t N0=par[knri*3]/par[0];

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
        returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*randcoinf2n;
        //! random 1n decay of daugter
        returnval+=neueff*pn[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);
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
    //! construct params
    Double_t lowerlimit=-10;
    Double_t upperlimit=10;

    makepath("testinput.txt");
    //! Define function without neutron gate
    TF1* fB=new TF1("fB",fcn_decay,lowerlimit,upperlimit,knri*3+2);
    fB->SetNpx(2000);
    fB->SetLineWidth(2);
    fB->SetLineColor(8);

    // initializing parameters
    for (Int_t i=0;i<knri*3;i++){
        fB->SetParameter(i,parms[i]);
        fB->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (!isparmsfix[i]) fB->FixParameter(i,parms[i]);
    }
    fB->SetParameter(knri*3,1000);//inital activity
    fB->SetParLimits(knri*3,0,10000);

    fB->SetParameter(knri*3+1,2000);//background
    fB->SetParLimits(knri*3+1,0,10000);

    //! Define function with 1 neutron gate
    TF1* fSB=new TF1("fSB",fcn_1ndecay,lowerlimit,upperlimit,knri*3+4);
    fSB->SetNpx(2000);
    fSB->SetLineWidth(2);
    fSB->SetLineColor(8);

    // initializing parameters
    for (Int_t i=0;i<knri*3;i++){
        fSB->SetParameter(i,parms[i]);
        fSB->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (!isparmsfix[i]) fSB->FixParameter(i,parms[i]);
    }
    fSB->SetParameter(knri*3,1000);//inital activity
    fSB->SetParLimits(knri*3,0,10000);

    fSB->SetParameter(knri*3+1,500);//background
    fSB->SetParLimits(knri*3+1,0,10000);

    fSB->SetParameter(knri*3+2,0.02);//random 1 neutron factor
    fSB->SetParLimits(knri*3+2,0,1);

    fSB->SetParameter(knri*3+3,0.022);//random gt1 neutron factor
    fSB->SetParLimits(knri*3+3,0,1);


    //! Define function with 2 neutrons gate
    TF1* fSB2=new TF1("fSB2",fcn_2ndecay,lowerlimit,upperlimit,knri*3+5);
    fSB2->SetNpx(2000);
    fSB2->SetLineWidth(2);
    fSB2->SetLineColor(8);

    // initializing parameters
    for (Int_t i=0;i<knri*3;i++){
        fSB2->SetParameter(i,parms[i]);
        fSB2->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (!isparmsfix[i]) fSB2->FixParameter(i,parms[i]);
    }
    fSB2->SetParameter(knri*3,1000);//inital activity
    fSB2->SetParLimits(knri*3,0,10000);

    fSB2->SetParameter(knri*3+1,100);//background
    fSB2->SetParLimits(knri*3+1,0,10000);

    fSB2->SetParameter(knri*3+2,0.02);//random 1 neutron factor
    fSB2->SetParLimits(knri*3+2,0,1);

    fSB2->SetParameter(knri*3+3,0.022);//random gt1 neutron factor
    fSB2->SetParLimits(knri*3+3,0,1);

    fSB2->SetParameter(knri*3+4,0.005);//random 2 neutrons factor
    fSB2->SetParLimits(knri*3+4,0,1);

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

