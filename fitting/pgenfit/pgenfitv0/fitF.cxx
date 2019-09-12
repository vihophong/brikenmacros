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
        cout<<"parameter "<<i<<" isconstant "<<_p[i]->getAttribute("Constant")<<endl;
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
     return ret ;

 }





