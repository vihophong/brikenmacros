#include <fstream>
#include <sstream>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TRandom3.h>

Double_t fnc_ef(Double_t *x, Double_t *pars) {

    double f, g, r, f1, f2, f3, x1, x2, y1, y2, y3;
    double dg = 0.0, d1, df1 = 0.0, df2 = 0.0;

      /* this eval is for use with effit */
      /* calculate the fit using present values of the pars */

      x1 = log(x[0] / pars[7]);
      x2 = log(x[0] / pars[8]);
      g  = pars[6];
      f1 = pars[0] + pars[1]*x1 + pars[2]*x1*x1;
      f2 = pars[3] + pars[4]*x2 + pars[5]*x2*x2;
      if (f1 <= f2) {
        f = f1;
        r = f1 / f2;
      } else {
        f = f2;
        r = f2 / f1;
      }
      if (r <= 1e-6) {
        df1 = exp(f);
        return pars[9]*exp(f);

      } else {
        y1 = pow(r, g);
        y2 = y1 + 1.0;
        d1 = -1.0 / g;
        y3 = pow(y2, d1);
        f3 = exp(f * y3);
        return pars[9]*f3;
      }
}


Double_t errmc(TF1* fcn,Double_t x,Int_t niter=10000,Double_t sig=3)
{
    Int_t nbins=200;
    TF1* fcnclone=(TF1*)fcn->Clone();
    Double_t centroid=fcnclone->Eval(x);
    Double_t min=centroid*(1-sig);
    Double_t max=centroid*(1+sig);

    TH1F* hgaus=new TH1F("h1","",nbins,min,max);

    Int_t npar=fcn->GetNpar();
    Double_t mcparms[npar];
    Double_t parms[npar];
    Double_t parmserr[npar];

    for (Int_t i=0;i<npar;i++){
        parms[i]=fcn->GetParameter(i);
        mcparms[i]=fcn->GetParameter(i);
        parmserr[i]=fcn->GetParError(i);
    }

    TRandom3* rseed=new TRandom3();
    for (Int_t j=0;j<niter;j++){

        TF1* f1=new TF1("f1",fnc_ef,1,2000,10);

        for (Int_t i=0;i<npar;i++){
            mcparms[i]=rseed->Gaus(parms[i],parmserr[i]);
            f1->SetParameter(i,mcparms[i]);
        }
        hgaus->Fill(f1->Eval(x));
    }
    hgaus->Draw();
    hgaus->Fit("gaus","LRE0","",min,max);
    return hgaus->GetFunction("gaus")->GetParameter(2);
}
