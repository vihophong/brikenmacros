/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

// Your description goes here...

#include "Riostream.h"

#include "fitF.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"


ClassImp(fitF)

 fitF::fitF(const char *name, const char *title,
                        RooAbsReal& _x,
                        RooAbsCategory& _y,
                        RooAbsReal& _ineueff,
                        char* parmsfilename
                        ) :
   RooAbsPdf(name,title),
   x("x","x",this,_x),
   y("y","y",this,_y),
   ineueff("ineueff","ineueff",this,_ineueff)
 {    
    makepath(parmsfilename);

    cout<<"\n\n\n********************************\nInitializing Unbinned likelihood P.D.F\n********************************\n\n\n"<<endl;

    //! default random coincidence factor
    parms[knri*4+4]=0.;//random 1 neutron factor
    parms[knri*4+5]=0.;//random gt0 neutron factor
    parms[knri*4+6]=0.;//random 2 neutrons factor
    parmsmin[knri*4+4]=0;parmsmax[knri*4+4]=1;
    parmsmin[knri*4+5]=0;parmsmax[knri*4+5]=1;
    parmsmin[knri*4+6]=0;parmsmax[knri*4+6]=1;
    isparmsfix[knri*4+4]=1;
    isparmsfix[knri*4+5]=1;
    isparmsfix[knri*4+6]=1;

    //! Decay parameters from input files
    int k=0;
    for (int i=0;i<knri*4;i++){
        if (isparmsfix[i]!=2){
            cout<<"Parms for unbinned fit, index1 = "<<i<<"\tindex2 = "<<k<<" : ";
            cout<<riname[i]<<"\t"<<parms[i]<<"\t"<<parmserr[i]<<"\t"<<parmsmin[i]<<"\t"<<parmsmax[i]<<"\t"<<isparmsfix[i]<<endl;
            //! initialize decay parameters (default value)
            _p[k]=new RooRealVar(Form("pp%d",k),Form("pp%d",k),parms[i],parmsmin[i],parmsmax[i]);
            if (isparmsfix[i]==1) _p[k]->setConstant();

            if (flag_sum_isomer_ratio){
                if (i>knri*3-1){
                    for (Int_t j=0;j<nisomers;j++){
                        if ((i-knri*3) == groundstate[j]){
                            _p[k]->setConstant();
                        }
                    }
                }
            }


            //! initialize roorealproxy
            p[k]=new RooRealProxy(Form("p%d",k),Form("p%d",k),this,*_p[k]);
            k++;
        }
    }
    for (Int_t i=0;i<k;i++){
        cout<<"parameter "<<i<<" isconstant"<<_p[k]->getAttribute("Constant")<<endl;
    }
    //! default values for random coincidence parameters
    //randcoinf1n
    _p[k]=new RooRealVar(Form("pp%d",k),Form("pp%d",k),0.000001,0.,1.);
    //randcoinfgt0n
    _p[k+1]=new RooRealVar(Form("pp%d",k+1),Form("pp%d",k+1),0.000001,0.,1.);
    //randcoinf2n
    _p[k+2]=new RooRealVar(Form("pp%d",k+2),Form("pp%d",k+2),0.000001,0.,1.);



    //! initialize roorealproxy
    p[k]=new RooRealProxy(Form("p%d",k),Form("p%d",k),this,*_p[k]);
    p[k+1]=new RooRealProxy(Form("p%d",k+1),Form("p%d",k+1),this,*_p[k+1]);
    p[k+2]=new RooRealProxy(Form("p%d",k+2),Form("p%d",k+2),this,*_p[k+2]);


    for (Int_t i=0;i<50;i++){
        _p[k+3+i]=new RooRealVar(Form("pp%d",k+3+i),Form("pp%d",k+3+i),0.000001,0.,1.);
        p[k+3+i]=new RooRealProxy(Form("p%d",k+3+i),Form("p%d",k+3+i),this,*_p[k+3+i]);
    }

 }

 fitF::fitF(const fitF& other, const char* name) :
   RooAbsPdf(other,name),
   x("x",this,other.x),
   y("y",this,other.y),
   ineueff("ineueff",this,other.ineueff)
 {
     for (Int_t i=0;i<nparmsactive+3;i++){
         p[i]=new RooRealProxy(Form("p%d",i),this,*other.p[i]);
     }
 }


 Double_t fitF::evaluate() const
 {     
     neueff = ineueff;
     Double_t par0[knri*4+3];
     for (Int_t i=0;i<knri*4;i++){
         if (isparmsfix[i]!=2){
             par0[i]=*p[indexparmsactive[i]];
         }else{
             par0[i]=parms[i];
         }
     }

     //! initial activity set to 1, background(off+slope) =0
     par0[knri*4] = 1;
     par0[knri*4+1] = 0;
     par0[knri*4+2] = 0;//slope

     Double_t par1[knri*4+5];
     for (Int_t i=0;i<knri*4;i++){
         if (isparmsfix[i]!=2){
             par1[i]=*p[indexparmsactive[i]];
         }else{
             par1[i]=parms[i];
         }
     }
     //! initial activity set to 1, background(off+slope) =0
     par1[knri*4] = 1;
     par1[knri*4+1] = 0;
     par1[knri*4+4] = 0;//slope

     //! default random coincidence factor
     par1[knri*4+2] = *p[nparmsactive];//randcoinf1n
     par1[knri*4+3] = *p[nparmsactive+1];//randcoinfgt0n


     Double_t par2[knri*4+6];
     for (Int_t i=0;i<knri*4;i++){
         if (isparmsfix[i]!=2){
             par2[i]=*p[indexparmsactive[i]];
         }else{
             par2[i]=parms[i];
         }
     }
     //! initial activity set to 1, background(off+slope) =0
     par2[knri*4] = 1;
     par2[knri*4+1] = 0;
     par2[knri*4+5] = 0;//slope

     //! default random coincidence factor
     par2[knri*4+2] = *p[nparmsactive];//randcoinf1n
     par2[knri*4+3] = *p[nparmsactive+1];//randcoinfgt0n
     par2[knri*4+4] = *p[nparmsactive+2];//randcoinf2n

     Double_t t[1];
     t[0]=x;

     Double_t ret=0;
     if (y==0){
         ret = fcn_decay(t,par0)-fcn_1ndecay(t,par1)-fcn_2ndecay(t,par2);
     }else if (y==1){
         ret = fcn_1ndecay(t,par1);
     }else{
         ret = fcn_2ndecay(t,par2);
     }
     return ret;
 }

 Double_t fitF::fcn_decay(Double_t *x, Double_t *par) const{
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
 Double_t fitF::fcn_1ndecay(Double_t *x, Double_t *par) const{

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
     returnval+=ineueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf1n-randcoinfgt0n);

     //! decay with 1 neutron of parent from p2n
     returnval+=2*(ineueff*(1-ineueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf1n-randcoinfgt0n);
     //! decay with 2 neutron of parent (not random 1 neutron)
     returnval-=ineueff*ineueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;

     for (Int_t i=0;i<npaths;i++){
         //! random coinc of beta decay of daugter
         returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*randcoinf1n;
         //! decay with 1 neutron of daugter nuclei
         returnval+=ineueff*pn[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(1-randcoinf1n-randcoinfgt0n);

         //! decay with 1 neutron of daugter from p2n
         returnval+=2*(ineueff*(1-ineueff))*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(1-randcoinf1n-randcoinfgt0n);
         //! decay with 2 neutron of parent (not random 1 neutron)
         returnval-=ineueff*ineueff*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*randcoinf1n;
     }

     return returnval;

 }


 //! Function for decay with 2 neutron emission
 Double_t fitF ::fcn_2ndecay(Double_t *x, Double_t *par) const{
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
     returnval+=ineueff*ineueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf2n-randcoinfgt0n);

     //! random coinc of beta decay of parent
     returnval+=lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf2n;
     //! random 1n decay of parent
     returnval+=ineueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

     //! decay with 1 neutron from P2n of parent - randomly correlated
     returnval+=2*(ineueff*(1-ineueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

     for (Int_t i=0;i<npaths;i++){
         //! decay with 2 neutron from P2n of daugter
         returnval+=ineueff*ineueff*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(1-randcoinf2n-randcoinfgt0n);
         //! random coinc of beta decay of daugter
         returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*randcoinf2n;
         //! random 1n decay of daugter
         returnval+=ineueff*pn[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

         //! decay with 1 neutron from P2n of daugter - randomly correlated
         returnval+=2*(ineueff*(1-ineueff))*p2n[decaymap[i][ndecay[i]-1]]*lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);
     }
     return returnval;
 }


 //! Global Bateaman function
 Double_t fitF::corefcn(Int_t ndecay,Int_t*  idecaymap,Int_t*  inneu,Double_t* production_yield, Double_t* b1n,Double_t* b2n,Double_t* lamda,Double_t N0,Double_t t) const{
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


