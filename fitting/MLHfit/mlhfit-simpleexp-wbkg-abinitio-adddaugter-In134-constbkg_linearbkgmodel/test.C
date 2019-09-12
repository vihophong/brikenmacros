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
#include "fitF.cxx"
#include "fitFbkg.cxx"
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

#include <fstream>

using namespace RooFit ;
using namespace RooStats;


using namespace std;


void test()
{
     RooRealVar* p[22];
     for (int i=0;i<22;i++){
         p[i]=new RooRealVar(Form("p%d",i),Form("p%d",i),1,0,2);
     }
     RooRealVar x("x","x",0,10) ;
     //! define discrete variable y
     RooCategory y("y","y");
     y.defineType("0neu",0);
     y.defineType("1neu",1);
     y.defineType("2neu",2);
     RooRealVar neueff("neueff","neueff",0.62,0.,1.) ;
     neueff.setConstant(kTRUE);
     fitF totdecaymodel("totdecaymodel","totdecaymodel",x,y,neueff,*p);

     // Import tree to total decay data
     TFile *f=TFile::Open("outhist.root");
     TTree* tree;
     f->GetObject("tree",tree);
     //! prepare data set for fitting forward correlated data
     RooDataSet* data=new RooDataSet("data","data",RooArgSet(x,y),Import(*tree)) ;
     data->Print() ;
     // Fit
     RooFitResult* fitres;
     fitres=totdecaymodel.fitTo(*data,NumCPU(2),Save()) ;


}
