/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

#ifndef FITF
#define FITF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "fitFunctions.h"

const Int_t kmaxnparsinit = 0;
Int_t nparsinit = 0;
Int_t indexparsinit[kmaxnparsinit];

class fitF : public RooAbsPdf {
public:
  fitF() {}
  fitF(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsCategory& _y,
              RooAbsReal& _ineueff,
              list<RooAbsReal *> listparms
              );
  fitF(const fitF& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new fitF(*this,newname); }
  inline virtual ~fitF() { }

  int npar;
protected:

  RooRealProxy x ;
  RooCategoryProxy y ;
  RooRealProxy ineueff ;

  RooRealProxy *p[kmaxnparsinit];

  Double_t evaluate() const ; 

private:  

  ClassDef(fitF,1) // Your description goes here...
};

#endif
