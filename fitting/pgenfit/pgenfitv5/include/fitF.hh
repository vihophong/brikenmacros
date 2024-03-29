/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

#ifndef FITF
#define FITF


#include "Riostream.h"

#include <fstream>

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "common.hh"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"

#include "TStopwatch.h"

#define MAX_N_PARMS 200
class fitF : public RooAbsPdf {
public:
  fitF() {} ;
  fitF(const char *name, const char *title,
              RooAbsReal& _x, RooAbsCategory &_y,
              RooAbsReal* _pp[]);
  fitF(const fitF& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new fitF(*this,newname); }
  inline virtual ~fitF() {
      std::ifstream pathfile("path.txt");
      Int_t nri;
      pathfile>>nri;
      for (int i=0;i<nri*5+4;i++)
          delete p[i];
  }
  Double_t fcndecay(Double_t *x, Double_t *par) const;
  Double_t fcndecay1n(Double_t *x, Double_t *par) const;
  Double_t fcndecay2n(Double_t *x, Double_t *par) const;

  //! for plotting
  Double_t fcndecay_parent(Double_t *x, Double_t *par) const;
  Double_t fcndecay1n_parent(Double_t *x, Double_t *par) const;
  Double_t fcndecay2n_parent(Double_t *x, Double_t *par) const;
  Double_t fcndecay_daugter(Double_t *x, Double_t *par) const{
      return fcndecay(x, par)-fcndecay_parent(x, par)+par[fpath->nri5+8]+par[fpath->nri5+9]*x[0];
  }
  Double_t fcndecay1n_daugter(Double_t *x, Double_t *par) const{
      return fcndecay1n(x, par)-fcndecay1n_parent(x, par)+par[fpath->nri5+8]+par[fpath->nri5+9]*x[0];
  }
  Double_t fcndecay2n_daugter(Double_t *x, Double_t *par) const{
      return fcndecay2n(x, par)-fcndecay2n_parent(x, par)+par[fpath->nri5+8]+par[fpath->nri5+9]*x[0];
  }

  Double_t fcndecay1n_c1(Double_t *x, Double_t *par) const;
  Double_t fcndecay1n_c2(Double_t *x, Double_t *par) const;
  Double_t fcndecay1n_c3(Double_t *x, Double_t *par) const;
  Double_t fcndecay1n_c23(Double_t *x, Double_t *par) const;

  Double_t fcndecay2n_c1(Double_t *x, Double_t *par) const;
  Double_t fcndecay2n_c2(Double_t *x, Double_t *par) const;
  Double_t fcndecay2n_c3(Double_t *x, Double_t *par) const;
  Double_t fcndecay2n_c4(Double_t *x, Double_t *par) const;

  Double_t fcndecay2n_c134(Double_t *x, Double_t *par) const;


  void initPath() ;

protected:
  path* fpath;
  void calculateDecay(Double_t &fdecayall,Double_t &fparent, Double_t* fdaugters, Double_t *l,Double_t *e,Double_t *p1n,Double_t *p2n,Double_t *py,Double_t N0,Double_t be) const;
  Double_t calculateDecay1n(Double_t fparent, Double_t* fdaugters, Double_t *p1n,Double_t *p2n,Double_t *ne,Double_t randcoinf1n,Double_t randcoinfgt0n,Double_t be,Double_t b1ne,Double_t b2ne,Double_t n1n2ne) const;
  Double_t calculateDecay2n(Double_t fparent, Double_t* fdaugters, Double_t *p1n,Double_t *p2n,Double_t *ne,Double_t randcoinf1n,Double_t randcoinfgt0n,Double_t be,Double_t randcoinf2n,Double_t b1ne,Double_t b2ne,Double_t n1n2ne) const;

  RooRealProxy x ;
  RooCategoryProxy y ;
  RooRealProxy* p[MAX_N_PARMS];
  Double_t evaluate() const ;


private:
  ClassDef(fitF,1) // Your description goes here...
};

class fitFbkg : public RooAbsPdf {
public:
  fitFbkg() {} ;
  fitFbkg(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsCategory& _y,
              RooAbsReal& _bkg1,
              RooAbsReal& _bkg2,
              RooAbsReal& _slope1,
              RooAbsReal& _slope2,
              RooAbsReal& _slope3);
  fitFbkg(const fitFbkg& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new fitFbkg(*this,newname); }
  inline virtual ~fitFbkg() { }
protected:

  RooRealProxy x ;
  RooCategoryProxy y ;
  RooRealProxy bkg1 ;
  RooRealProxy bkg2 ;
  RooRealProxy slope1 ;
  RooRealProxy slope2 ;
  RooRealProxy slope3 ;
  Double_t evaluate() const ;

private:

  ClassDef(fitFbkg,1) // Your description goes here...
};

#endif
