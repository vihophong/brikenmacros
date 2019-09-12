/*Fitter with Full decay chain
 *A Linear background is adoptted
 *Nov 6, 2018
 * Update fitting the negative background
 * Plot for residual and decay components
*/

#include "TFrame.h"
#include "TBox.h"
#include "makepath.C"
//Double_t neueff=0.66*(100-0.8)/100;
Double_t neueff=0.62;//changed to 62 %

Bool_t flag_sum_isomer_ratio=true;

//! Global Bateaman function
Double_t corefcn(Int_t ndecay,Int_t*  idecaymap,Int_t*  inneu,Double_t* production_yield, Double_t* b1n,Double_t* b2n,Double_t* lamda,Double_t N0,Double_t t){
    Double_t fcnret=0;

    Double_t factor1=1.;
    //only parrent decay p2n
    for (int i=0;i<ndecay-1;i++){
        Double_t corrproductionyield = production_yield[idecaymap[i+1]];

        if (flag_sum_isomer_ratio){
            for (Int_t j=0;j<nisomers;j++){
                if (idecaymap[i+1] == groundstate[j]){
                    corrproductionyield = 1-production_yield[isomerstate[j]];
                }
            }
        }

        if (inneu[i]==0){
            factor1=factor1 * corrproductionyield*(1-b1n[idecaymap[i]]-b2n[idecaymap[i]])*lamda[idecaymap[i]];//branching here!
        }else if (inneu[i]==1){
            factor1=factor1 * corrproductionyield*b1n[idecaymap[i]]*lamda[idecaymap[i]];
        }else{
            factor1=factor1 * corrproductionyield*b2n[idecaymap[i]]*lamda[idecaymap[i]];
        }

    }

    Double_t factor2=0;
    for (int i=0;i<ndecay;i++){
        Double_t factor2i=exp(-lamda[idecaymap[i]]*t);

        Double_t factor2ij=1;
        for (int j=0;j<ndecay;j++)
            if (j!=i) factor2ij=factor2ij*(lamda[idecaymap[j]]-lamda[idecaymap[i]]);
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
    Double_t* production_yield=&par[knri*3];

    //bkg
    Double_t returnval=par[knri*4+1]+par[knri*4+2]*x[0];
    //init
    Double_t N0=par[knri*4]/par[0];

    //! Parent nuclei
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0]);

    for (Int_t i=0;i<npaths;i++){        
        returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0]);
    }
    return returnval;
}

//! Function for decay with 1 neutron emission
Double_t fcn_1ndecay(Double_t *x, Double_t *par) {

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t* p2n=&par[knri*2];
    Double_t* production_yield=&par[knri*3];

    //! pairs isomeric states


    Double_t randcoinfgt0n=par[knri*4+3];
    Double_t randcoinf1n=par[knri*4+2];
    //bkg
    Double_t returnval=par[knri*4+1]+par[knri*4+4]*x[0];
    //init
    Double_t N0=par[knri*4]/par[0];

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
        returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*randcoinf1n;
        //! decay with 1 neutron of daugter nuclei
        returnval+=neueff*pn[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(1-randcoinf1n-randcoinfgt0n);

        //! decay with 1 neutron of daugter from p2n
        returnval+=2*(neueff*(1-neueff))*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(1-randcoinf1n-randcoinfgt0n);
        //! decay with 2 neutron of parent (not random 1 neutron)
        returnval-=neueff*neueff*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*randcoinf1n;
    }

    return returnval;

}


//! Function for decay with 2 neutron emission
Double_t fcn_2ndecay(Double_t *x, Double_t *par) {
    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t* p2n=&par[knri*2];
    Double_t* production_yield=&par[knri*3];

    Double_t randcoinf2n=par[knri*4+4];
    Double_t randcoinfgt0n=par[knri*4+3];
    Double_t randcoinf1n=par[knri*4+2];
    //bkg
    Double_t returnval=par[knri*4+1]+par[knri*4+5]*x[0];
    //init
    Double_t N0=par[knri*4]/par[0];

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
        returnval+=neueff*neueff*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(1-randcoinf2n-randcoinfgt0n);
        //! random coinc of beta decay of daugter
        returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*randcoinf2n;
        //! random 1n decay of daugter
        returnval+=neueff*pn[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

        //! decay with 1 neutron from P2n of daugter - randomly correlated
        returnval+=2*(neueff*(1-neueff))*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);
    }
    return returnval;
}
