/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

#ifndef FITFBKG
#define FITFBKG

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class fitFbkg : public RooAbsPdf {
public:
  fitFbkg() {} ;
  fitFbkg(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsCategory& _y,
              RooAbsReal& _bkg1,
              RooAbsReal& _bkg2);
  fitFbkg(const fitFbkg& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new fitFbkg(*this,newname); }
  inline virtual ~fitFbkg() { }

protected:

  RooRealProxy x ;
  RooCategoryProxy y ;
  RooRealProxy bkg1 ;
  RooRealProxy bkg2 ;

  Double_t evaluate() const ;

private:

  ClassDef(fitFbkg,1) // Your description goes here...
};

#endif
