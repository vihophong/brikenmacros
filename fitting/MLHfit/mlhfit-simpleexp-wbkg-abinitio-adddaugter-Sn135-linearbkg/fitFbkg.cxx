/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

// Your description goes here...

#include "Riostream.h"

#include "fitFbkg.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"

ClassImp(fitFbkg)

 fitFbkg::fitFbkg(const char *name, const char *title,
                        RooAbsReal& _x,
                        RooAbsCategory& _y,
                        RooAbsReal& _bkg1,
                        RooAbsReal& _bkg2,
                        RooAbsReal& _slope) :
   RooAbsPdf(name,title),
   x("x","x",this,_x),
   y("y","y",this,_y),
   bkg1("bkg1","bkg1",this,_bkg1),
   bkg2("bkg2","bkg2",this,_bkg2),
   slope("slope","slope",this,_slope)
 {
 }


 fitFbkg::fitFbkg(const fitFbkg& other, const char* name) :
   RooAbsPdf(other,name),
   x("x",this,other.x),
   y("y",this,other.y),
   bkg1("bkg1",this,other.bkg1),
   bkg2("bkg2",this,other.bkg2),
   slope("slope",this,other.slope)
 {
 }
 Double_t fitFbkg::evaluate() const
 {
     Double_t ret=0;
     if (y==0){
         ret = (1-(1+slope*x)*(bkg1+bkg2*bkg1));
     }else if (y==1){
         ret = bkg1*(1+slope*x);
     }else{
         ret = bkg2*bkg1*(1+slope*x);
     }
     return ret ;
 }


