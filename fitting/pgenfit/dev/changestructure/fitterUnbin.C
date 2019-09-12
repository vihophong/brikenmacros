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
Long64_t nentrieslimit = -1;//set negative value for fitting to all entries

Double_t p_deadtime=0.;
Double_t p_timerange=10;

using namespace std;
void fitterUnbin()
{
    // Import tree to total decay data
    TFile *f=TFile::Open("../oldcompare-unbin/outhist.root");
    TTree* tree;
    f->GetObject("tree",tree);

    if (nentrieslimit>0) tree->SetEntries(nentrieslimit);

    RooRealVar x("x","x",p_deadtime,p_timerange) ;
    RooRealVar xbkg("x","x",-p_timerange,-p_deadtime) ;

    //! define discrete variable y
    RooCategory y("y","y");
    y.defineType("0neu",0);
    y.defineType("1neu",1);
    y.defineType("2neu",2);


    RooRealVar ineueff("ineueff","ineueff",0.62,0.,1.) ;
    ineueff.setConstant(kTRUE);



    makepath(parmsfilename);

    cout<<"\n\n\n********************************\nInitializing Unbinned likelihood P.D.F\n********************************\n\n\n"<<endl;

    RooRealVar* _p[kmaxnparsinit];
    //! Decay parameters from input files
    int k=0;
    for (int i=0;i<knri*4;i++){
        if (isparmsfix[i]!=2){
            cout<<"Parms for unbinned fit, index1 = "<<i<<"\tindex2 = "<<k<<" : ";
            if (i<knri) cout<<riname[i];
            cout<<"\t"<<parms[i]<<"\t"<<parmserr[i]<<"\t"<<parmsmin[i]<<"\t"<<parmsmax[i]<<"\t"<<isparmsfix[i]<<endl;
            //! initialize decay parameters (default value)
            _p[k]=new RooRealVar(Form("pp%d",k),Form("pp%d",k),parms[i],parmsmin[i],parmsmax[i]);
            if (isparmsfix[i]==1) _p[i]->setConstant();
            //! initialize roorealproxy
            p[k]=new RooRealProxy(Form("p%d",k),Form("p%d",k),this,*_p[k]);
            k++;
        }
    }
    //! default values for random coincidence parameters
    //randcoinf1n
    _p[k]=new RooRealVar(Form("pp%d",k),Form("pp%d",k),0.000001,0.,1.);
    //randcoinfgt0n
    _p[k+1]=new RooRealVar(Form("pp%d",k+1),Form("pp%d",k+1),0.000001,0.,1.);
    //randcoinf2n
    _p[k+2]=new RooRealVar(Form("pp%d",k+2),Form("pp%d",k+2),0.000001,0.,1.);

    list<RooAbsReal*> listofpars;
    listofpars.emplace(listofpars.end(),(RooAbsReal*)par1);
    listofpars.emplace(listofpars.end(),(RooAbsReal*)par2);

    fitF totdecaymodel("totdecaymodel","totdecaymodel",x,y,ineueff,listofpars);
    totdecaymodel.npar=5;
    nparsinit = listofpars.size();


    //! Toy MC dataset (not working!, just for testing)
    RooDataSet* data = totdecaymodel.generate(RooArgSet(x,y),100) ;

}
