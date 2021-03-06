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


