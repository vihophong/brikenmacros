/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

// Your description goes here...


#include "fitF.hh"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"

ClassImp(fitF)

 fitF::fitF(const char *name, const char *title,
                        RooAbsReal& _x,
                        RooAbsCategory& _y,
                        RooAbsReal& _neueff,
                        RooAbsReal *_pp[]) :
   RooAbsPdf(name,title),
   x("x","x",this,_x),
   y("y","y",this,_y),
   neueff("neueff","neueff",this,_neueff)
 {
    std::ifstream pathfile("path.txt");
    Int_t nri;
    pathfile>>nri;
    pathfile.close();
    //std::cout<<nri<<std::endl;
    for (Int_t i=0;i<nri*5+4;i++){
        p[i]=new RooRealProxy(Form("p%i",i),Form("p%i",i),this,*_pp[i]);
    }
 }


 fitF::fitF(const fitF& other, const char* name) :
   RooAbsPdf(other,name),
   x("x",this,other.x),
   y("y",this,other.y),
   neueff("neueff",this,other.neueff)
 {
     fpath=new path;
     std::ifstream pathfile("path.txt");
     pathfile>>fpath->nri;
     pathfile>>fpath->npaths;
     for (int i=0;i<fpath->npaths;i++){
         pathfile>>fpath->ndecay[i];
         for (int j=0;j<fpath->ndecay[i];j++){
             pathfile>>fpath->decaymap[i][j]>>fpath->nneu[i][j];
         }
     }
     pathfile.close();
//     std::cout<<fpath->nri<<std::endl;
//     std::cout<<fpath->npaths<<std::endl;
//     for (int i=0;i<fpath->npaths;i++){
//         std::cout<<fpath->ndecay[i]<<std::endl;
//         for (int j=0;j<fpath->ndecay[i];j++){
//             std::cout<<fpath->decaymap[i][j]<<"\t"<<fpath->nneu[i][j]<<std::endl;
//         }
//     }

     for (Int_t i=0;i<fpath->nri*5+4;i++)
         p[i]=new RooRealProxy(Form("p%i",i),this,*other.p[i]);
 }



 Double_t fitF::evaluate() const
 {

     Double_t par0[fpath->nri*5+4];
     Double_t par1[fpath->nri*5+4];
     Double_t par2[fpath->nri*5+4];
     for (Int_t i=0;i<fpath->nri*5+4;i++)
     {
         par0[i]=*p[i];
         par1[i]=*p[i];
         par2[i]=*p[i];
     }

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
     Double_t* pn=&par[fpath->nri];
     Double_t* lamda=par;
     Double_t* p2n=&par[fpath->nri*2];
     Double_t* production_yield=&par[fpath->nri*3];

     //init
     Double_t N0=par[fpath->nri*5]/par[0];

     //! Parent nuclei
     Double_t returnval=lamda[0]*N0*exp(-lamda[0]*x[0]);

     for (Int_t i=0;i<fpath->npaths;i++){
         returnval+=lamda[fpath->decaymap[i][fpath->ndecay[i]-1]]*corefcn(fpath->ndecay[i],fpath->decaymap[i],fpath->nneu[i],production_yield,pn,p2n,lamda,N0,x[0]);
     }
     return returnval;
 }

 //! Function for decay with 1 neutron emission
 Double_t fitF::fcn_1ndecay(Double_t *x, Double_t *par) const{
     Double_t* pn=&par[fpath->nri];
     Double_t* lamda=par;
     Double_t* p2n=&par[fpath->nri*2];
     Double_t* production_yield=&par[fpath->nri*3];
     Double_t* neuefffactor=&par[fpath->nri*4];
     Double_t randcoinfgt0n=par[fpath->nri*5+2];
     Double_t randcoinf1n=par[fpath->nri*5+1];
     //init
     Double_t N0=par[fpath->nri*5]/par[0];

     //! Parent nuclei

     //! random coinc of beta decay of parent
     Double_t returnval=lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;
     //! decay with 1 neutron of parent
     returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf1n-randcoinfgt0n);

     //! decay with 1 neutron of parent from p2n
     returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf1n-randcoinfgt0n);
     //! decay with 2 neutron of parent (not random 1 neutron)
     returnval-=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf1n;

     for (Int_t i=0;i<fpath->npaths;i++){
         Double_t neueffcorr=neueff*neuefffactor[fpath->decaymap[i][fpath->ndecay[i]-1]];
         //! random coinc of beta decay of daugter
         returnval+=lamda[fpath->decaymap[i][fpath->ndecay[i]-1]]*corefcn(fpath->ndecay[i],fpath->decaymap[i],fpath->nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*randcoinf1n;
         //! decay with 1 neutron of daugter nuclei
         returnval+=neueffcorr*pn[fpath->decaymap[i][fpath->ndecay[i]-1]]*lamda[fpath->decaymap[i][fpath->ndecay[i]-1]]*corefcn(fpath->ndecay[i],fpath->decaymap[i],fpath->nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(1-randcoinf1n-randcoinfgt0n);

         //! decay with 1 neutron of daugter from p2n
         returnval+=2*(neueffcorr*(1-neueffcorr))*p2n[fpath->decaymap[i][fpath->ndecay[i]-1]]*lamda[fpath->decaymap[i][fpath->ndecay[i]-1]]*corefcn(fpath->ndecay[i],fpath->decaymap[i],fpath->nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(1-randcoinf1n-randcoinfgt0n);
         //! decay with 2 neutron of parent (not random 1 neutron)
         returnval-=neueffcorr*neueffcorr*p2n[fpath->decaymap[i][fpath->ndecay[i]-1]]*lamda[fpath->decaymap[i][fpath->ndecay[i]-1]]*corefcn(fpath->ndecay[i],fpath->decaymap[i],fpath->nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*randcoinf1n;
     }

     return returnval;

 }


 //! Function for decay with 2 neutron emission
 Double_t fitF::fcn_2ndecay(Double_t *x, Double_t *par) const{
     Double_t* pn=&par[fpath->nri];
     Double_t* lamda=par;
     Double_t* p2n=&par[fpath->nri*2];
     Double_t* production_yield=&par[fpath->nri*3];
     Double_t* neuefffactor=&par[fpath->nri*4];
     Double_t randcoinf2n=par[fpath->nri*5+3];
     Double_t randcoinfgt0n=par[fpath->nri*5+2];
     Double_t randcoinf1n=par[fpath->nri*5+1];

     //init
     Double_t N0=par[fpath->nri*5]/par[0];

     //! parent
     //! decay with 2 neutron from P2n of parent
     Double_t returnval=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(1-randcoinf2n-randcoinfgt0n);

     //! random coinc of beta decay of parent
     returnval+=lamda[0]*N0*exp(-lamda[0]*x[0])*randcoinf2n;
     //! random 1n decay of parent
     returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

     //! decay with 1 neutron from P2n of parent - randomly correlated
     returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

     for (Int_t i=0;i<fpath->npaths;i++){
         Double_t neueffcorr=neueff*neuefffactor[fpath->decaymap[i][fpath->ndecay[i]-1]];

         //! decay with 2 neutron from P2n of daugter
         returnval+=neueffcorr*neueffcorr*p2n[fpath->decaymap[i][fpath->ndecay[i]-1]]*lamda[fpath->decaymap[i][fpath->ndecay[i]-1]]*corefcn(fpath->ndecay[i],fpath->decaymap[i],fpath->nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(1-randcoinf2n-randcoinfgt0n);
         //! random coinc of beta decay of daugter
         returnval+=lamda[fpath->decaymap[i][fpath->ndecay[i]-1]]*corefcn(fpath->ndecay[i],fpath->decaymap[i],fpath->nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*randcoinf2n;
         //! random 1n decay of daugter
         returnval+=neueffcorr*pn[fpath->decaymap[i][fpath->ndecay[i]-1]]*lamda[fpath->decaymap[i][fpath->ndecay[i]-1]]*corefcn(fpath->ndecay[i],fpath->decaymap[i],fpath->nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

         //! decay with 1 neutron from P2n of daugter - randomly correlated
         returnval+=2*(neueffcorr*(1-neueffcorr))*p2n[fpath->decaymap[i][fpath->ndecay[i]-1]]*lamda[fpath->decaymap[i][fpath->ndecay[i]-1]]*corefcn(fpath->ndecay[i],fpath->decaymap[i],fpath->nneu[i],production_yield,pn,p2n,lamda,N0,x[0])*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);
     }
     return returnval;
 }

 Double_t fitF::corefcn(Int_t ndecay,Int_t*  idecaymap,Int_t*  inneu,Double_t* production_yield, Double_t* b1n,Double_t* b2n,Double_t* lamda,Double_t N0,Double_t t) const{
     Double_t fcnret=0;

     Double_t factor1=1.;
     //only parrent decay p2n
     for (int i=0;i<ndecay-1;i++){
         Double_t corrproductionyield = production_yield[idecaymap[i+1]];

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

 ClassImp(fitFbkg)

  fitFbkg::fitFbkg(const char *name, const char *title,
                         RooAbsReal& _x,
                         RooAbsCategory& _y,
                         RooAbsReal& _bkg1,
                         RooAbsReal& _bkg2,
                         RooAbsReal& _slope1,
                         RooAbsReal& _slope2,
                         RooAbsReal& _slope3) :
    RooAbsPdf(name,title),
    x("x","x",this,_x),
    y("y","y",this,_y),
    bkg1("bkg1","bkg1",this,_bkg1),
    bkg2("bkg2","bkg2",this,_bkg2),
    slope1("slope1","slope1",this,_slope1),
    slope2("slope2","slope2",this,_slope2),
    slope3("slope3","slope3",this,_slope3)
  {
  }


  fitFbkg::fitFbkg(const fitFbkg& other, const char* name) :
    RooAbsPdf(other,name),
    x("x",this,other.x),
    y("y",this,other.y),
    bkg1("bkg1",this,other.bkg1),
    bkg2("bkg2",this,other.bkg2),
    slope1("slope1",this,other.slope1),
    slope2("slope2",this,other.slope2),
    slope3("slope3",this,other.slope3)
  {
  }
  Double_t fitFbkg::evaluate() const
  {
      Double_t ret=0;
      if (y==0){
          ret = slope3*x+1-(1+slope1*x)*bkg1-(1+slope2*x)*bkg2*bkg1;
      }else if (y==1){
          ret = bkg1*(1+slope1*x);
      }else{
          ret = bkg2*bkg1*(1+slope2*x);
      }
      return ret ;
  }


