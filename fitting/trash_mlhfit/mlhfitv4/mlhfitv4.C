
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


using namespace RooFit ;
using namespace RooStats;


Double_t p_negtime=0.;
Double_t p_deadtime=0.;
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
        //cout<<knri*2+i+1<<"\t"<<parms[knri*2+i+1]<<"\t"<<parmsmin[knri*2+i+1]<<"\t"<<parmsmax[knri*2+i+1]<<endl;
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


void  mlhfitv4(char* fitname, char* infile,char* parmsfile,char* outfile)
{
    Int_t fitoption=0;
    //! Get parameters
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

    //! define variables with negative+positive
    RooRealVar x("x","x",p_negtime,p_timerange) ;
    //! define discrete variable y
    TFile *f = TFile::Open(infile);
    TTree* tree;
    f->GetObject("tree",tree);
    char tempchar1[500];

    //! Declare roorealvars
    RooRealVar* p[knri*2+4];

    // decay parameters
    for (int i=0;i<knri*2+1;i++){
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




    // form factor parameters
    RooRealVar nsignal("nsignal","number of nsignal events",100,0,10000);
    RooRealVar ntotalbkg("ntotalbkg","number of bkgs events",100,0,10000);


    //!models
    fitF model("signal","signal",x,*p[0],*p[1],*p[2],*p[3],*p[4],*p[5],*p[6],*p[7],*p[8],*p[9],*p[10],*p[11],*p[12],*p[13],*p[14],*p[15],*p[16],*p[17],*p[18],nsignal,ntotalbkg);

    //! data
    RooDataSet* data=new RooDataSet("data","data",x,Import(*tree)) ;
    data->Print() ;

    //! P e r f o r m   f i t s
    RooFitResult* fitres;
    fitres=model.fitTo(*data,NumCPU(24),Save()) ;


    //! Print and plot results for comparison
    fitres->Print() ;

    // Plot x distribution of data and projection of model on x = Int(dy) model(x,y)
    RooPlot* xframe0 = x.frame(Title("decay fit")) ;
    data->plotOn(xframe0,Binning(200)) ;
    model.plotOn(xframe0) ;

    RooHist* hresid0 = xframe0->residHist() ;
    RooPlot* rframe0 = x.frame(Title("decay fit Residual Distribution"));
    rframe0->addPlotable(hresid0,"P");

    //! book file and tree
    TFile* fout=new TFile(outfile,"recreate");
    fout->cd();
    // Make canvas and draw RooPlots

    sprintf(tempchar1,"mlhfitresult%s",fitname);
    TCanvas *c = new TCanvas(tempchar1,tempchar1,1200, 800);
    c->Divide(2,1);
    c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe0->GetYaxis()->SetTitleOffset(1.4) ; xframe0->Draw() ;
    c->cd(2) ; gPad->SetLeftMargin(0.15) ; rframe0->GetYaxis()->SetTitleOffset(1.4) ; rframe0->Draw() ;
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
}
