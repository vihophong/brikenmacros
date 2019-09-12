#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooCategory.h"
#include "TRandom3.h"
#include "TFile.h"
#include "RooGaussian.h"
#include "RooConstVar.h"


#include "RooRandom.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooAddition.h"
#include "RooDataHist.h"
#include "RooPoisson.h"
#include "RooPlot.h"

#include "RooNLLVar.h"

#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "RooStats/ProofConfig.h"
#include "RooStats/ToyMCSampler.h"
#include "RooProfileLL.h"

#include "fitF.cxx"
#include "fitFbkg.cxx"


#include <fstream>

using namespace RooFit ;
using namespace RooStats;

const Int_t ncpu = 8;
Long64_t nentrieslimit = 5000;//set negative value for fitting to all entries

Double_t p_deadtime=0.;
Double_t p_timerange=10;

using namespace std;
int test()
{
    // Import tree to total decay data

    RooRealVar x("x","x",p_deadtime,p_timerange) ;
    RooRealVar xbkg("x","x",-p_timerange,-p_deadtime) ;

    //! define discrete variable y
    RooCategory y("y","y");
    y.defineType("0neu",0);
    y.defineType("1neu",1);
    y.defineType("2neu",2);

    RooRealVar ineueff("ineueff","ineueff",0.62,0.,1.) ;
    ineueff.setConstant(kTRUE);

    fitF totdecaymodel("totdecaymodel","totdecaymodel",x,y,ineueff,"In134short.txt");
    //! Toy MC dataset (not working!, just for testing)
    RooDataSet* data = totdecaymodel.generate(RooArgSet(x,y),5000) ;
    //cout<<totdecaymodel.test.id <<endl;
    return 0;
}
