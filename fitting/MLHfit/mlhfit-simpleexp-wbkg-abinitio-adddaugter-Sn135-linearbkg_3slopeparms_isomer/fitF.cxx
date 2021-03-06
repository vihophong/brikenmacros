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
                        RooAbsReal& _neueff,
                        RooAbsReal * _p) :
   RooAbsPdf(name,title),
   x("x","x",this,_x),
   y("y","y",this,_y),
   neueff("neueff","neueff",this,_neueff)
 {
    for (Int_t i=0;i<kmaxparms;i++){
        p[i] = new RooRealProxy(Form("p%d",i),Form("p%d",i),this,_p[i]);
    }
 }


 fitF::fitF(const fitF& other, const char* name) :
   RooAbsPdf(other,name),
   x("x",this,other.x),
   y("y",this,other.y),
   neueff("neueff",this,other.neueff)
 {
     for (Int_t i=0;i<kmaxparms;i++){
         p[i] = new RooRealProxy(Form("p%d",i),Form("p%d",i),this,other.p[i]);
     }
 }



 Double_t fitF::evaluate() const
 {     
     Double_t par0[19];
     for (Int_t i=0;i<19;i++){
        par0[0]=*p[0];
     }
     Double_t par1[21];
     for (Int_t i=0;i<21;i++){
        par1[0]=*p[0];
     }
     Double_t par2[22];
     for (Int_t i=0;i<22;i++){
        par2[0]=*p[0];
     }

     Double_t ret=0;
     if (y==0){
         ret = fcn_gen(x,par0)-fcn_gen_w1neutron(x,par1)-fcn_gen_w2neutron(x,par2);
     }else if (y==1){
         ret = fcn_gen_w1neutron(x,par1);
     }else{
         ret = fcn_gen_w2neutron(x,par2);
     }
     return ret ;
 }

 //! Global function
 Double_t fitF::fcn_gen(Double_t t, Double_t *par) const{
     Int_t npaths;
     Int_t ndecay[kmaxpaths];
     Int_t decaymap[kmaxpaths][kmaxndecay];
     Int_t nneu[kmaxpaths][kmaxndecay];
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

     Double_t returnval=0;

     Double_t* pn=&par[knri];
     Double_t* lamda=par;
     Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
     p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
     Double_t N0=1./par[0];

     //! Parent nuclei
     returnval+=lamda[0]*N0*exp(-lamda[0]*t);


     for (Int_t i=0;i<npaths;i++){
         returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,t);
     }

     return returnval;
 }

 Double_t fitF::fcn_gen_w1neutron(Double_t t, Double_t *par) const{
     Int_t npaths;
     Int_t ndecay[kmaxpaths];
     Int_t decaymap[kmaxpaths][kmaxndecay];
     Int_t nneu[kmaxpaths][kmaxndecay];
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

     Double_t randcoinfgt0n=par[knri*2+2];
     Double_t randcoinf1n=par[knri*2+1];

     Double_t returnval=0.;

     Double_t* pn=&par[knri];
     Double_t* lamda=par;
     Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
     p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
     Double_t N0=1./par[0];

     //! Parent nuclei
     //! decay with 1 neutron of parent
     returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*t)*(1-randcoinf1n-randcoinfgt0n);
     //! decay with 1 neutron of parent
     returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*t)*(1-randcoinf1n-randcoinfgt0n);

     //! decay with 2 neutron of parent (not random 1 neutron)
     returnval-=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*t)*randcoinf1n;

     //! random coinc part
     returnval+=lamda[0]*N0*exp(-lamda[0]*t)*randcoinf1n;


     for (Int_t i=0;i<npaths;i++){
         //! random coinc of beta decay of daugter
         returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,t)*randcoinf1n;
         //! decay with 1 neutron part of daugter nuclei
         returnval+=neueff*pn[decaymap[i][ndecay[i]-1]-1]*lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,t)*(1-randcoinf1n-randcoinfgt0n);
     }

     return returnval;
 }


 Double_t fitF::fcn_gen_w2neutron(Double_t t, Double_t *par) const{
     Int_t npaths;
     Int_t ndecay[kmaxpaths];
     Int_t decaymap[kmaxpaths][kmaxndecay];
     Int_t nneu[kmaxpaths][kmaxndecay];
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

     Double_t randcoinf2n=par[knri*2+3];
     Double_t randcoinfgt0n=par[knri*2+2];
     Double_t randcoinf1n=par[knri*2+1];
     Double_t returnval=0;

     Double_t* pn=&par[knri];
     Double_t* lamda=par;
     Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
     p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
     Double_t N0=1./par[0];

     //! parent
     //! decay with 2 neutron from P2n of parent
     returnval+=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*t)*(1-randcoinf2n-randcoinfgt0n);

     //! random 1n decay of parent
     returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*t)*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);
     //! decay with 1 neutron from P2n of parent
     returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*t)*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);

     //! random coinc part
     returnval+=lamda[0]*N0*exp(-lamda[0]*t)*randcoinf2n;


     for (Int_t i=0;i<npaths;i++){
         //! random coinc of beta decay of daugter
         returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,t)*randcoinf2n;
         //! random 1n decay of daugter
         returnval+=neueff*pn[decaymap[i][ndecay[i]-1]-1]*lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,t)*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);
     }

     return returnval;
 }



 Double_t fitF::corefcn(Int_t ndecay,Int_t* decaymap,Int_t* nneu, Double_t* b1n,Double_t* b2n,Double_t* lamda,Double_t N0,Double_t t) const{
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





