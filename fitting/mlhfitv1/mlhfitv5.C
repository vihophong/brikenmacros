/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
///  'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #312
///
///  Performing fits in multiple (disjoint) ranges in one or more dimensions
///
/// \macro_output
/// \macro_code
/// \author 07/2008 - Wouter Verkerke


#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooCategory.h"
#include "TRandom3.h"
#include "TFile.h"
#include "fitF.cxx"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "fitFbkg.cxx"


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


using namespace RooFit ;
using namespace RooStats;

Double_t p_deadtime=0.0001;
Double_t p_timerange=10.0;
Double_t nsigma=10;


void getparms(Double_t* parms,Double_t* parmserr, Double_t* parmsmax, Double_t* parmsmin, Bool_t* isparmsfix,char* infile,Double_t nsig=10.)
{
    std::ifstream ifs(infile);
    Int_t rino;
    Double_t temp;
    Bool_t flagfix[knri][3];
    Double_t decayparms[knri][3];
    Double_t decayparms_err[knri][3];

    for (int i=0;i<knri;i++){
        ifs>>rino;
        for (int j=0;j<3;j++){
            ifs>>temp;
            if (temp>=0){
                flagfix[i][j]=true;
                decayparms[i][j]=temp;
            }else{
                flagfix[i][j]=false;
                decayparms[i][j]=(-temp);
            }
            ifs>>temp;decayparms_err[i][j]=temp;
        }
    }

    for (int i=0;i<knri;i++){
        for (int j=0;j<3;j++){
            if (j==0){//for half-life
                decayparms_err[i][j]=log(2)/decayparms[i][j]/decayparms[i][j]*decayparms_err[i][j];
                decayparms[i][j]=log(2)/decayparms[i][j];
            }
        }
    }


    //! calculate output
    for (int i=0;i<knri;i++){
        for (int j=0;j<2;j++){
            isparmsfix[j*knri+i]=flagfix[i][j];
            parms[j*knri+i]=decayparms[i][j];
            parmserr[j*knri+i]=decayparms_err[i][j];
            parmsmax[j*knri+i]=decayparms[i][j]+decayparms_err[i][j]*nsig;
            if ((decayparms[i][j]-decayparms_err[i][j]*nsig)>0)
                parmsmin[j*knri+i]=decayparms[i][j]-decayparms_err[i][j]*nsig;
            else
                parmsmin[j*knri+i]=log(2)/100000000000;
            //for pn
            if (j!=0) {parmsmin[j*knri+i]=0;parmsmax[j*knri+i]=1;}
        }
    }

    //! read p2n parrent
    parms[knri*2]=decayparms[0][2];
    parmserr[knri*2]=decayparms_err[0][2];
    parms[knri*2]=decayparms[0][2];
    parmsmin[knri*2]=0;parmsmax[knri*2]=1;
    isparmsfix[knri*2]=flagfix[0][2];

    //! read initial activity, background and random coincidence factor
    for (int i=0;i<5;i++){
        ifs>>parms[knri*2+i+1]>>parmsmin[knri*2+i+1]>>parmsmax[knri*2+i+1];
        parmserr[knri*2+i+1]=parms[knri*2+i+1]-parmsmin[knri*2+i+1];
        cout<<knri*2+i+1<<"\t"<<parms[knri*2+i+1]<<"\t"<<parmsmin[knri*2+i+1]<<"\t"<<parmsmax[knri*2+i+1]<<endl;
        isparmsfix[knri*2+i+1]=false;
        if (i==4) {
            parmsmin[knri*2+i+1]=0;parmsmax[knri*2+i+1]=1;
            if (parms[knri*2+i+1]>=0){
                isparmsfix[knri*2+i+1]=true;
            }else{
                parms[knri*2+i+1]=-parms[knri*2+i+1];
                isparmsfix[knri*2+i+1]=false;
            }
        }
    }
}


void getbackground(RooRealVar* nsig,RooRealVar* nbkg1,RooRealVar* nbkg2, char* infile)
{
    TFile *f=TFile::Open(infile);
    TTree* tree;
    f->GetObject("tree",tree);
    Double_t nall=tree->Draw("x",Form("x<%f&&x>%f",-p_deadtime,-p_timerange),"goff");
    Double_t n1ngate=tree->Draw("x",Form("x<%f&&x>%f&&y==1",-p_deadtime,-p_timerange),"goff");
    Double_t n2ngate=tree->Draw("x",Form("x<%f&&x>%f&&y==2",-p_deadtime,-p_timerange),"goff");
    nbkg1->setVal(n1ngate/nall);
    nbkg1->setMax(n1ngate/nall+n1ngate/nall*nsigma);
    nbkg1->setMin(n1ngate/nall-n1ngate/nall*nsigma);
    nbkg1->setError(n1ngate/nall*sqrt(1/n1ngate+1/nall));

    nbkg2->setVal(n2ngate/n1ngate);
    nbkg2->setMax(n2ngate/n1ngate+n2ngate/n1ngate*nsigma);
    nbkg2->setMin(n2ngate/n1ngate-n2ngate/n1ngate*nsigma);
    nbkg2->setError(n2ngate/n1ngate*sqrt(1/n2ngate+1/n1ngate));

    nsig->setVal(nall);
    nsig->setMax(nall+nall*nsigma);
    nsig->setMin(nall-nall*nsigma);
    nsig->setError(sqrt(nall));
    f->Close();
}
void getbackgroundcheck(char* infile){
    RooRealVar ntotalbkg("ntotalbkg","number of bkgs events",1e2,10,2e6);
    RooRealVar bkgratio1("bkg1","bkg1 nevt",0.5,0.,1.);
    RooRealVar bkgratio2("bkg2","bkg2 nevt",0.5,0.,1.);
    getbackground(&ntotalbkg,&bkgratio1,&bkgratio2,infile);
    cout<<"nsig = "<<ntotalbkg.getVal()<<" +/- "<<ntotalbkg.getError()<<endl;
    cout<<"bkg1ratio = "<<bkgratio1.getVal()<<" +/- "<<bkgratio1.getError()<<endl;
    cout<<"bkg2ratio = "<<bkgratio2.getVal()<<" +/- "<<bkgratio2.getError()<<endl;
}


void mlhfitv5(char* fitname, char* infile,char* parmsfile,char* outfile,Int_t fitoption=0)
{
    //! parameters
    Double_t parms[knri*2+6];
    Double_t parmserr[knri*2+6];
    Double_t parmsmax[knri*2+6];
    Double_t parmsmin[knri*2+6];
    Bool_t isparmsfix[knri*2+6];
    getparms(parms,parmserr,parmsmax,parmsmin,isparmsfix,parmsfile,nsigma);

    for (int i=0;i<knri*2+6;i++){
        cout<<"parm "<<i<<"\t"<<parms[i]<<"\t"<<parmserr[i]<<"\t"<<parmsmin[i]<<"\t"<<parmsmax[i]<<"\t"<<isparmsfix[i]<<endl;
    }
    cout<<"\n\n**************Setting parameters....\n"<<endl;

    //! define variables
    RooRealVar x("x","x",p_deadtime,p_timerange) ;

    RooRealVar* p[knri*2+6];
    for (int i=0;i<knri*2+6;i++){
        p[i]=new RooRealVar(Form("p%d",i),Form("p%d",i),parms[i],parmsmin[i],parmsmax[i]);
        if ((fitoption!=1&&fitoption!=2)||parmserr[i]==0){
            if (isparmsfix[i]){
                cout<<"fixed value p"<<i<<"\tval="<<parms[i]<<endl;
                p[i]->setConstant(kTRUE);
            }else{
                cout<<"variable value p"<<i<<"\tval "<<parms[i]<<"\tmin="<<parmsmin[i]<<"\tmax="<<parmsmax[i]<<endl;
            }
        }else{
            cout<<"variable value p"<<i<<"\tval "<<parms[i]<<"\tmin="<<parmsmin[i]<<"\tmax="<<parmsmax[i]<<endl;
        }
    }

    //! define discrete variable y
    RooCategory y("y","y");
    y.defineType("0neu",0);
    y.defineType("1neu",1);
    y.defineType("2neu",2);

    //! define variable for neutron detection efficinecy
    //RooRealVar neueff("neueff","neueff",0.613384,0.,1.) ;
    RooRealVar neueff("neueff","neueff",0.68,0.,1.) ;
    neueff.setConstant(kTRUE);

    //! define random coincidence factor
    /*
    p[knri*2+5]=new RooRealVar(Form("p%d",knri*2+5),Form("p%d",knri*2+5),parms[knri*2+5],0,1);
    if (parms[knri*2+5]>=0){
        p[knri*2+5]->setVal(parms[knri*2+5]);
        p[knri*2+5]->setConstant(kTRUE);
    }*/

    //! define roogaussian for error propagation
    RooGaussian* pconstr[19];

    RooArgSet constronlydecay;
    RooArgSet constronlydecaywbkg;

    for (int i=0;i<knri*2+1;i++){
        pconstr[i]=new RooGaussian(Form("p%dconstr",i),Form("p%dconstr",i),*p[i],RooConst(parms[i]),RooConst(parmserr[i]));
        if (isparmsfix[i]&&parmserr[i]>0&&(fitoption==1||fitoption==2)){
            cout<<"set constrain model for p"<<i<<"\tmean="<<parms[i]<<"\tstd="<<parmserr[i]<<endl;
            constronlydecay.add(*pconstr[i]);
            constronlydecaywbkg.add(*pconstr[i]);
        }
    }

    //! fix val for extended model
    p[19]->setMin(0);
    p[19]->setMax(2);
    p[19]->setVal(1);
    p[20]->setMin(-1);
    p[20]->setMax(1);
    p[20]->setVal(0);
    p[21]->setMin(-1);
    p[21]->setMax(1);
    p[21]->setVal(0);
    p[22]->setMin(-1);
    p[22]->setMax(1);
    p[22]->setVal(0);
    p[19]->setConstant(kTRUE);
    p[20]->setConstant(kTRUE);
    p[21]->setConstant(kTRUE);
    p[22]->setConstant(kTRUE);

    RooRealVar nsignal("nsignal","number of nsignal events",parms[19],parmsmin[19],parmsmax[19]);
    RooRealVar ntotalbkg("ntotalbkg","number of bkgs events",parms[20],parmsmin[20],parmsmax[20]);
    RooRealVar bkgratio1("bkg1","bkg1 nevt",0.5,0.,1.);
    RooRealVar bkgratio2("bkg2","bkg2 nevt",0.5,0.,1.);

    getbackground(&ntotalbkg,&bkgratio1,&bkgratio2,infile);
    //cout<<"BACKGROUND!"<<endl;


    /*
    //! add constrain here!
    RooGaussian* pconstrtotalbkg=new RooGaussian("pconstrtotalbkg","pconstrtotalbkg",ntotalbkg,RooConst(ntotalbkg.getVal()),RooConst(ntotalbkg.getError()));
    RooGaussian* pconstrbkgratio1=new RooGaussian("pconstrbkgratio1","pconstrbkgratio1",bkgratio1,RooConst(bkgratio1.getVal()),RooConst(bkgratio1.getError()));
    RooGaussian* pconstrbkgratio2=new RooGaussian("pconstrbkgratio2","pconstrbkgratio2",bkgratio2,RooConst(bkgratio2.getVal()),RooConst(bkgratio2.getError()));
    constronlydecaywbkg.add(*pconstrtotalbkg);
    constronlydecaywbkg.add(*pconstrbkgratio1);
    constronlydecaywbkg.add(*pconstrbkgratio2);
    */

    //! fix negative background

    ntotalbkg.setConstant(kTRUE);
    bkgratio1.setConstant(kTRUE);
    bkgratio2.setConstant(kTRUE);


    //! define p.d.f fit model
    fitF model("signal","signal",x,y,neueff,*p[0],*p[1],*p[2],*p[3],*p[4],*p[5],*p[6],*p[7],*p[8],*p[9],*p[10],*p[11],*p[12],*p[13],*p[14],*p[15],*p[16],*p[17],*p[18],*p[19],*p[20],*p[21],*p[22],*p[23]);

    fitFbkg bkgmodel("bkgmodel","bkgmodel",x,y,bkgratio1,bkgratio2);

    RooAddPdf final_pdf("final_pdf","final pdf",RooArgList(model,bkgmodel),RooArgList(nsignal,ntotalbkg));


    // Sample 100000 events in (x,y) from the model
    //RooDataSet* modelData = model.generate(RooArgSet(x,y),10000) ;

    // I m p o r t i n g   i n t e g e r   T T r e e   b r a n c h e s
    // ---------------------------------------------------------------

    // Import integer tree branch as RooCategory
    // will be imported as those are the only defined states
    TFile *f=TFile::Open(infile);
    TTree* tree;
    f->GetObject("tree",tree);
    RooDataSet* data=new RooDataSet("data","data",RooArgSet(x,y),Import(*tree)) ;
    data->Print() ;

    // P e r f o r m   f i t s   i n   i n d i v i d u a l   s i d e b a n d   r e g i o n s
    // -------------------------------------------------------------------------------------

    // Perform fit in SideBand1 region (RooAddPdf coefficients will be interpreted in full range)
    //model.fitTo(*modelData) ;


    RooFitResult* fitres;


    if (fitoption==1)
        fitres=final_pdf.fitTo(*data,Extended(),ExternalConstraints(constronlydecay),NumCPU(24),Save()) ;
    else if (fitoption==2)
        fitres=final_pdf.fitTo(*data,Extended(),ExternalConstraints(constronlydecaywbkg),NumCPU(24),Save()) ;
    else if (fitoption==0)
        fitres=final_pdf.fitTo(*data,Extended(),NumCPU(24),Save()) ;


     // Print results for comparison
    //r_sb1->Print() ;

    // Plot x distribution of data and projection of model on x = Int(dy) model(x,y)
    RooPlot* xframe0 = x.frame(Title("0 neutron fit")) ;
    //modelData->plotOn(xframe0,Cut("y==y::0neu")) ;
    data->plotOn(xframe0,Cut("y==y::0neu"),Binning(1000)) ;
    final_pdf.plotOn(xframe0,Slice(y,"0neu")) ;

    RooPlot* xframe1 = x.frame(Title("1 neutron fit")) ;
    //modelData->plotOn(xframe1,Cut("y==y::1neu")) ;
    data->plotOn(xframe1,Cut("y==y::1neu"),Binning(1000)) ;
    final_pdf.plotOn(xframe1,Slice(y,"1neu")) ;

    RooPlot* xframe2 = x.frame(Title("2 neutrons fit")) ;
    //modelData->plotOn(xframe2,Cut("y==y::2neu")) ;
    data->plotOn(xframe2,Cut("y==y::2neu"),Binning(1000)) ;
    final_pdf.plotOn(xframe2,Slice(y,"2neu")) ;

    RooHist* hresid0 = xframe0->residHist() ;
    RooPlot* rframe0 = x.frame(Title("0 neutron fit Residual Distribution"));
    rframe0->addPlotable(hresid0,"P");

    RooHist* hresid1 = xframe1->residHist() ;
    RooPlot* rframe1 = x.frame(Title("1 neutron fit Residual Distribution"));
    rframe1->addPlotable(hresid1,"P");

    RooHist* hresid2 = xframe2->residHist() ;
    RooPlot* rframe2 = x.frame(Title("2 neutrons fit Residual Distribution"));
    rframe2->addPlotable(hresid2,"P");

    //! book file and tree
    TFile* fout=new TFile(outfile,"recreate");
    fout->cd();
    // Make canvas and draw RooPlots

    char tempchar1[500];
    sprintf(tempchar1,"mlhfitresult%s",fitname);
    TCanvas *c = new TCanvas(tempchar1,tempchar1,1200, 800);
    c->Divide(2,3);
    c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe0->GetYaxis()->SetTitleOffset(1.4) ; xframe0->Draw() ;
    c->cd(3) ; gPad->SetLeftMargin(0.15) ; xframe1->GetYaxis()->SetTitleOffset(1.4) ; xframe1->Draw() ;
    c->cd(5) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.4) ; xframe2->Draw() ;
    c->cd(2) ; gPad->SetLeftMargin(0.15) ; rframe0->GetYaxis()->SetTitleOffset(1.4) ; rframe0->Draw() ;
    c->cd(4) ; gPad->SetLeftMargin(0.15) ; rframe1->GetYaxis()->SetTitleOffset(1.4) ; rframe1->Draw() ;
    c->cd(6) ; gPad->SetLeftMargin(0.15) ; rframe2->GetYaxis()->SetTitleOffset(1.4) ; rframe2->Draw() ;
    c->Write();
    fout->Close();
    c->Close();

    cout<<"\n********************SUMARRY***************"<<endl;
    cout<<"T1/2= "<<log(2)/p[0]->getVal()<<" +/- "<<log(2)/p[0]->getVal()/p[0]->getVal()*p[0]->getError()<<endl;
    cout<<"P1n= "<<p[9]->getVal()<<" +/- "<<p[9]->getError()<<endl;
    cout<<"P2n= "<<p[18]->getVal()<<" +/- "<<p[18]->getError()<<endl;
    cout<<endl;

    cout<<log(2)/p[0]->getVal()<<"\t"<<log(2)/p[0]->getVal()/p[0]->getVal()*p[0]->getError()<<endl;
    cout<<p[9]->getVal()<<"\t"<<p[9]->getError()<<endl;
    cout<<p[18]->getVal()<<"\t"<<p[18]->getError()<<endl;

    //! write to text file
    std::ofstream ofs("mlhfitresult.txt", std::ofstream::out | std::ofstream::app);
    ofs<<fitname<<",";
    ofs<<log(2)/p[0]->getVal()<<","<<log(2)/p[0]->getVal()/p[0]->getVal()*p[0]->getError()<<",";
    ofs<<p[9]->getVal()<<","<<p[9]->getError()<<",";
    ofs<<p[18]->getVal()<<","<<p[18]->getError()<<","<<p_deadtime<<","<<p_timerange<<endl;

    cout<<p[18]->getErrorLo()<<" -EEE- "<<p[18]->getErrorHi()<<endl;

}



