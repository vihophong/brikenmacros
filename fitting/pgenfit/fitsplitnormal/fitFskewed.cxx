/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "fitFskewed.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 
#include <boost/math/distributions/skew_normal.hpp>

ClassImp(fitFskewed); 

 fitFskewed::fitFskewed(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _mean,
                        RooAbsReal& _sig1,
                        RooAbsReal& _sig2) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mean("mean","mean",this,_mean),
   sig1("sig1","sig1",this,_sig1),
   sig2("sig2","sig2",this,_sig2)
 { 
 } 


 fitFskewed::fitFskewed(const fitFskewed& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mean("mean",this,other.mean),
   sig1("sig1",this,other.sig1),
   sig2("sig2",this,other.sig2)
 { 
 } 



 Double_t fitFskewed::evaluate() const 
 { 
   //auto sn=boost::math::skew_normal_distribution<double>(mean,sig1,sig2);
   //return boost::math::quantile(sn,x);

   double a=TMath::Sqrt(2/TMath::Pi())*(sig1+sig2);
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   return a*exp(-0.5*pow((x-mean)/(sig1+(x<mean)*sig2*(x-mean)),2)) ;
 } 



