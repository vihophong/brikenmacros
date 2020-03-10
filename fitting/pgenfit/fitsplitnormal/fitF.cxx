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

ClassImp(fitF); 

 fitF::fitF(const char *name, const char *title, 
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


 fitF::fitF(const fitF& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mean("mean",this,other.mean),
   sig1("sig1",this,other.sig1),
   sig2("sig2",this,other.sig2)
 { 
 } 



 Double_t fitF::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   double a=TMath::Sqrt(2/TMath::Pi())*(sig1+sig2);
   double left=a*exp(-(x-mean)*(x-mean)/2/sig1/sig1);
   double right=a*exp(-(x-mean)*(x-mean)/2/sig2/sig2);
   if (x<mean) return left;
   else return right; 
 } 


