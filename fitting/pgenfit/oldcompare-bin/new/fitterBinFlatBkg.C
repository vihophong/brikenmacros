/*Fitter with Full decay chain
 *A Linear background is adoptted
 *Nov 6, 2018
 * Update fitting the negative background
 * Plot for residual and decay components
*/

#include "TFrame.h"
#include "TBox.h"

#include "fitFunctions.h"

Double_t rejectrange=0.075;//first 1 ms

Double_t neueff_mean=0.605;
Double_t neueff_err=0.053;

// Definition of shared parameter
// Decay parameters
int iparB[kmaxparms];

// Decay 1 neutron parameters
int iparSB[kmaxparms];

// Decay 2 neutrons parameters
int iparSB2[kmaxparms];

// Create the GlobalCHi2 structure

struct GlobalChi2 {
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2,
                ROOT::Math::IMultiGenFunction & f3) :
      fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}

   // parameter vector is first background (in common 1 and 2)
   // and then is signal (only in 2)
   double operator() (const double *par) const {
      double p1[knri*4+3];
      for (int i = 0; i < knri*4+3; ++i) p1[i] = par[iparB[i] ];

      double p2[knri*4+5];
      for (int i = 0; i < knri*4+5; ++i) p2[i] = par[iparSB[i] ];

      double p3[knri*4+6];
      for (int i = 0; i < knri*4+6; ++i) p3[i] = par[iparSB2[i] ];

      return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3);
   }
   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
   const  ROOT::Math::IMultiGenFunction * fChi2_3;
};



//! Monte carlo variation
void mc(TRandom3* rseed)
{
    for (int i=0;i<knri*3;i++){
        if (isparmsfix[i]==1&&parmserr[i]>0.){
            mcparms[i]=rseed->Gaus(parms[i],parmserr[i]);
        }
    }
    //! mc for neutron efficiency?
    if (neueff_err>0){// only if there is error
      //neueff=neueff_mean-neueff_err+rseed->Rndm()*neueff_err*2;//for an uniform distribution
      neueff=rseed->Gaus(neueff_mean,neueff_err);//for a gausian distribution
    }
}


void fitter(char* infilename,char* parmsfilename,char* outfilename,Int_t ninterations=0,Int_t rebin=1)
{    
    Double_t lowerlimit=-10;
    Double_t upperlimit=20;

    //! input decay parameters and make decay path
    makepath(parmsfilename);


    //! special for isomeric case (ground state population is 1- isomeric state population)
    if (flag_sum_isomer_ratio){
        for (Int_t i=knri*3;i<knri*4;i++){
            for (Int_t j=0;j<nisomers;j++){
                if ((i-knri*3) == groundstate[j]){
                    isparmsfix[i]=2;
                }
            }
        }
    }

    TRandom3* rseed=new TRandom3;
    //! construct params

    Double_t bkgactmaxmin=0.50; //100% of max min bkg or initial activity
    Double_t plotrange[]={0.075,5.,-2};

    //! input default value for other parameters
    parms[knri*4]=1000;//initial activity
    parms[knri*4+1]=2000;//background of no neutron gate curve
    parms[knri*4+2]=500;//background of 1n gate curve
    parms[knri*4+3]=100;//background of 2n gate curve

    parms[knri*4+4]=0.02;//random 1 neutron factor
    parms[knri*4+5]=0.022;//random gt0 neutron factor
    parms[knri*4+6]=0.005;//random 2 neutrons factor

    parms[knri*4+7]=0.;//background slope of no neutron gate curve
    parms[knri*4+8]=0.;//background slope of 1n gate curve
    parms[knri*4+9]=0.;//background slope of 2n gate curve

    parmsmin[knri*4]=0;parmsmax[knri*4]=10000;
    parmsmin[knri*4+1]=0;parmsmax[knri*4+1]=10000;
    parmsmin[knri*4+2]=0;parmsmax[knri*4+2]=10000;
    parmsmin[knri*4+3]=0;parmsmax[knri*4+3]=10000;

    parmsmin[knri*4+4]=0;parmsmax[knri*4+4]=1;
    parmsmin[knri*4+5]=0;parmsmax[knri*4+5]=1;
    parmsmin[knri*4+6]=0;parmsmax[knri*4+6]=1;

    parmsmin[knri*4+7]=-1000;parmsmax[knri*4+7]=1000;
    parmsmin[knri*4+8]=-1000;parmsmax[knri*4+8]=1000;
    parmsmin[knri*4+9]=-1000;parmsmax[knri*4+9]=1000;

    isparmsfix[knri*4]=0;//inital activity - vary
    isparmsfix[knri*4+1]=1;//background of no neutron gate curve ?fix
    isparmsfix[knri*4+2]=1;//background of 1n gate curve ?fix
    isparmsfix[knri*4+3]=1;//background of 2n gate curve ?fix

    isparmsfix[knri*4+4]=1;
    isparmsfix[knri*4+5]=1;
    isparmsfix[knri*4+6]=1;

    isparmsfix[knri*4+7]=1;//background slope of no neutron gate curve
    isparmsfix[knri*4+8]=1;//background slope of 1n gate curve
    isparmsfix[knri*4+9]=1;//background slope of 2n gate curve

    //!****************************GET HISTOGRAM FROM FILE********************
    //!
    //!
    TFile *f = TFile::Open(infilename);
    char tempchar1[1000];
    sprintf(tempchar1,"hdecay");
    TH1F* hdecay=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecay1n");
    TH1F* hdecay1n=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecay2n");
    TH1F* hdecay2n=(TH1F*) gDirectory->Get(tempchar1);

    sprintf(tempchar1,"hdecay1nbwd");
    TH1F* hdecay1nbwd=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecaygt0nbwd");
    TH1F* hdecaygt0nbwd=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecay2nbwd");
    TH1F* hdecay2nbwd=(TH1F*) gDirectory->Get(tempchar1);


    //! Calculate random coincidence paramters
    Double_t n1nbwd=(Double_t) hdecay1nbwd->GetEntries();
    Double_t gt0nbwd=(Double_t) hdecaygt0nbwd->GetEntries();
    Double_t n2nbwd=(Double_t) hdecay2nbwd->GetEntries();
    Double_t nball=(Double_t) hdecay->GetEntries();

    parms[knri*4+4]=n1nbwd/nball;
    parms[knri*4+5]=gt0nbwd/nball;
    parms[knri*4+6]=n2nbwd/nball;

    //! Binning
    Int_t binning=hdecay->GetNbinsX();
    //!Rebin
    hdecay->Rebin(rebin);
    hdecay1n->Rebin(rebin);
    hdecay2n->Rebin(rebin);
    TH1F * hB = (TH1F*) hdecay->Clone();
    TH1F * hSB = (TH1F*) hdecay1n->Clone();
    TH1F * hSB2 = (TH1F*) hdecay2n->Clone();

    //! background parameter estimation
    //! background as average of several first bin

    //background of decay no neutron gate

    hdecay->Fit("pol0","LQRE0","goff",lowerlimit,-rejectrange);

    parms[knri*4+1]=hdecay->GetFunction("pol0")->GetParameter(0);
    parmserr[knri*4+1]=hdecay->GetFunction("pol0")->GetParError(0);
    parmsmin[knri*4+1]=parms[knri*4+1]-parms[knri*4+1]*bkgactmaxmin;
    parmsmax[knri*4+1]=parms[knri*4+1]+parms[knri*4+1]*bkgactmaxmin;

    parms[knri*4+7]=0;
    parmserr[knri*4+7]=0;

    //activity
    parms[knri*4]=hdecay->GetBinContent(hdecay->GetXaxis()->FindBin(rejectrange))-parms[knri*4+1];
    parmserr[knri*4]=parms[knri*4]*bkgactmaxmin;
    parmsmin[knri*4]=parms[knri*4]-parms[knri*4]*bkgactmaxmin;
    parmsmax[knri*4]=parms[knri*4]*2+parms[knri*4]*bkgactmaxmin;

    //background of decay 1 neutron gate

    hdecay1n->Fit("pol0","LQRE0","goff",lowerlimit,-rejectrange);

    parms[knri*4+2]=hdecay1n->GetFunction("pol0")->GetParameter(0);
    parmserr[knri*4+2]=hdecay1n->GetFunction("pol0")->GetParError(0);
    parmsmin[knri*4+2]=parms[knri*4+2]-parms[knri*4+2]*bkgactmaxmin;
    parmsmax[knri*4+2]=parms[knri*4+2]+parms[knri*4+2]*bkgactmaxmin;

    parms[knri*4+8]=0;
    parmserr[knri*4+8]=0;

    //background of decay 2 neutrons gate

    hdecay2n->Fit("pol0","LQRE0","goff",lowerlimit,-rejectrange);

    parms[knri*4+3]=hdecay2n->GetFunction("pol0")->GetParameter(0);
    parmserr[knri*4+3]=hdecay2n->GetFunction("pol0")->GetParError(0);
    parmsmin[knri*4+3]=parms[knri*4+3]-parms[knri*4+3]*bkgactmaxmin;
    parmsmax[knri*4+3]=parms[knri*4+3]+parms[knri*4+3]*bkgactmaxmin;

    parms[knri*4+9]=0;
    parmserr[knri*4+9]=0;

    //! ********** Define FITTING FUNCTION
    //! Define function without neutron gate
    TF1* fB=new TF1("fB",fcn_decay,lowerlimit,upperlimit,knri*4+3);
    fB->SetNpx(2000);
    fB->SetLineWidth(2);
    fB->SetLineColor(8);

    //! initializing parameters
    for (Int_t i=0;i<knri*4;i++){
        fB->SetParameter(i,parms[i]);
        fB->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (!isparmsfix[i]) fB->FixParameter(i,parms[i]);
    }
    fB->SetParameter(knri*4,parms[knri*4]);//inital activity
    fB->SetParLimits(knri*4,parmsmin[knri*4],parmsmax[knri*4]);

    fB->SetParameter(knri*4+1,parms[knri*4+1]);//background
    fB->SetParLimits(knri*4+1,parmsmin[knri*4+1],parmsmax[knri*4+1]);

    fB->SetParameter(knri*4+2,parms[knri*4+7]);//background slope
    fB->SetParLimits(knri*4+2,parmsmin[knri*4+7],parmsmax[knri*4+7]);

    //! Define function with 1 neutron gate
    TF1* fSB=new TF1("fSB",fcn_1ndecay,lowerlimit,upperlimit,knri*4+5);
    fSB->SetNpx(2000);
    fSB->SetLineWidth(2);
    fSB->SetLineColor(8);

    // initializing parameters
    for (Int_t i=0;i<knri*4;i++){
        fSB->SetParameter(i,parms[i]);
        fSB->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (!isparmsfix[i]) fSB->FixParameter(i,parms[i]);
    }
    fSB->SetParameter(knri*4,parms[knri*4]);//inital activity
    fSB->SetParLimits(knri*4,parmsmin[knri*4],parmsmax[knri*4]);

    fSB->SetParameter(knri*4+1,parms[knri*4+2]);//background
    fSB->SetParLimits(knri*4+1,parmsmin[knri*4+2],parmsmax[knri*4+2]);

    fSB->SetParameter(knri*4+2,parms[knri*4+4]);//random 1 neutron factor
    fSB->SetParLimits(knri*4+2,parmsmin[knri*4+4],parmsmax[knri*4+4]);

    fSB->SetParameter(knri*4+3,parms[knri*4+5]);//random gt0 neutron factor
    fSB->SetParLimits(knri*4+3,parmsmin[knri*4+5],parmsmax[knri*4+5]);

    fSB->SetParameter(knri*4+4,parms[knri*4+8]);//background slope
    fSB->SetParLimits(knri*4+4,parmsmin[knri*4+8],parmsmax[knri*4+8]);

    //! Define function with 2 neutrons gate
    TF1* fSB2=new TF1("fSB2",fcn_2ndecay,lowerlimit,upperlimit,knri*4+6);
    fSB2->SetNpx(2000);
    fSB2->SetLineWidth(2);
    fSB2->SetLineColor(8);

    // initializing parameters
    for (Int_t i=0;i<knri*4;i++){
        fSB2->SetParameter(i,parms[i]);
        fSB2->SetParLimits(i,parmsmin[i],parmsmax[i]);
        if (!isparmsfix[i]) fSB2->FixParameter(i,parms[i]);
    }
    fSB2->SetParameter(knri*4,parms[knri*4]);//inital activity
    fSB2->SetParLimits(knri*4,parmsmin[knri*4],parmsmax[knri*4]);

    fSB2->SetParameter(knri*4+1,parms[knri*4+3]);//background
    fSB2->SetParLimits(knri*4+1,parmsmin[knri*4+3],parmsmax[knri*4+3]);

    fSB2->SetParameter(knri*4+2,parms[knri*4+4]);//random 1 neutron factor
    fSB2->SetParLimits(knri*4+2,parmsmin[knri*4+4],parmsmax[knri*4+4]);

    fSB2->SetParameter(knri*4+3,parms[knri*4+5]);//random gt0 neutron factor
    fSB2->SetParLimits(knri*4+3,parmsmin[knri*4+5],parmsmax[knri*4+5]);

    fSB2->SetParameter(knri*4+4,parms[knri*4+6]);//random 2 neutrons factor
    fSB2->SetParLimits(knri*4+4,parmsmin[knri*4+6],parmsmax[knri*4+6]);

    fSB2->SetParameter(knri*4+5,parms[knri*4+9]);//background slope
    fSB2->SetParLimits(knri*4+5,parmsmin[knri*4+9],parmsmax[knri*4+9]);

    cout<<fB->Eval(1.)<<endl;
    cout<<fSB->Eval(1.)<<endl;
    cout<<fSB2->Eval(1.)<<endl;


    //! Book output file and tree output for systematic half-life estimation
    TFile* outfile=new TFile(outfilename,"recreate");

    Double_t outparms[knri*4+10];
    Double_t outparmserr[knri*4+10];
    Int_t iparms[knri*4+10];
    Int_t isvary[knri*4+10];
    for (int i=0;i<knri*4+10;i++){
        outparms[i]=0;
        outparmserr[i]=0;
        iparms[i]=i;
        isvary[i]=isparmsfix[i];
    }
    sprintf(tempchar1,"treemc");
    TTree* treeout=new TTree(tempchar1,tempchar1);
    treeout->Branch("outparms",outparms,Form("outparms[%d]/D",knri*4+10));
    treeout->Branch("outparmserr",outparmserr,Form("outparmserr[%d]/D",knri*4+10));
    treeout->Branch("iparms",iparms,Form("iparms[%d]/I",knri*4+10));
    treeout->Branch("neueff",&neueff,"neueff/D");
    treeout->Branch("isvary",outparmserr,Form("isvary[%d]/I",knri*4+10));


    //! Define Simultaneous fitting function
    //! initialzing shared prameter matrix for simultaneous fitting
    for (int i = 0; i < knri*4; i++) {
        iparB[i]=i;
        iparSB[i]=i;
        iparSB2[i]=i;
    }
    iparB[knri*4]=knri*4;//initial activity
    iparB[knri*4+1]=knri*4+1;//background of no neutron gate curve
    iparB[knri*4+2]=knri*4+7;//background slope of no neutron gate curve

    iparSB[knri*4]=knri*4;//initial activity
    iparSB[knri*4+1]=knri*4+2;//background of 1n gate curve
    iparSB[knri*4+4]=knri*4+8;//background slope of 1n gate curve

    iparSB2[knri*4]=knri*4;//initial activity
    iparSB2[knri*4+1]=knri*4+3;//background of 1n gate curve
    iparSB2[knri*4+5]=knri*4+9;//background slope of 1n gate curve

    iparSB[knri*4+2]=knri*4+4;//random 1 neutron factor
    iparSB[knri*4+3]=knri*4+5;//random gt0 neutron factor

    iparSB2[knri*4+2]=knri*4+4;//random 1 neutron factor
    iparSB2[knri*4+3]=knri*4+5;//random gt0 neutron factor
    iparSB2[knri*4+4]=knri*4+6;//random 2 neutrons factor

    ROOT::Math::WrappedMultiTF1 wfB(*fB,1);
    ROOT::Math::WrappedMultiTF1 wfSB(*fSB,1);
    ROOT::Math::WrappedMultiTF1 wfSB2(*fSB2,1);
    ROOT::Fit::DataOptions opt;

    // limit within the fitting range
    opt.fUseRange  =true;
    // set the data range
    ROOT::Fit::DataRange rangeB;
    rangeB.SetRange(rejectrange,upperlimit);
    ROOT::Fit::BinData dataB(opt,rangeB);
    ROOT::Fit::FillData(dataB, hB);

    ROOT::Fit::DataRange rangeSB;
    rangeSB.SetRange(rejectrange,upperlimit);
    ROOT::Fit::BinData dataSB(opt,rangeSB);
    ROOT::Fit::FillData(dataSB, hSB);

    ROOT::Fit::DataRange rangeSB2;
    rangeSB2.SetRange(rejectrange,upperlimit);
    ROOT::Fit::BinData dataSB2(opt,rangeSB2);
    ROOT::Fit::FillData(dataSB2, hSB2);

    ROOT::Fit::PoissonLLFunction chi2_B(dataB, wfB);
    ROOT::Fit::PoissonLLFunction chi2_SB(dataSB, wfSB);
    ROOT::Fit::PoissonLLFunction chi2_SB2(dataSB2, wfSB2);
    GlobalChi2 globalChi2(chi2_B, chi2_SB, chi2_SB2);

    ROOT::Fit::Fitter fitter;

    //! ***********SETTING PARAMETERS************

    cout<<"\n\n*****SETTING PARAMETERS......\n"<<endl;
    //Double_t arr1[parms.size()];
    //std::copy(parms.begin(),parms.end(),arr1);
    fitter.Config().SetParamsSettings(knri*4+10,parms);
    for (int i=0;i<knri*4+10;i++){
        if (isparmsfix[i]){
           cout<<"fixed value p"<<i<<"\tval="<<parms[i]<<endl;
           fitter.Config().ParSettings(i).Fix();
        }else{
           cout<<"variable value p"<<i<<"\tval "<<parms[i]<<"\tmin="<<parmsmin[i]<<"\tmax="<<parmsmax[i]<<endl;
           fitter.Config().ParSettings(i).SetLimits(parmsmin[i],parmsmax[i]);
        }
    }

    //! Setting MINOS package
    fitter.Config().SetMinimizer("Minuit2","Migrad");
    //fitter.Config().SetMinosErrors();
    if (fitter.Config().MinosErrors()) cout<<"minos enabled"<<endl;

    //!Perform the fit
    fitter.FitFCN(knri*4+10,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),false);


    //! Printting Results
    ROOT::Fit::FitResult result = fitter.Result();
    const Double_t* resultpar=result.GetParams();
    const Double_t* resulterr=result.GetErrors();


    //! fix monte carlo parameter center
    for (int i=0;i<knri*4+10;i++) mcparms[i]=parms[i];
    for(int i=0;i<ninterations;i++) {
        cout<<"mc "<<i+1<<endl;
        mc(rseed);
//        Double_t arr2[mcparms.size()];
//        std::copy(mcparms.begin(),mcparms.end(),arr2);
        fitter.Config().SetParamsSettings(knri*4+10,mcparms);
        for (unsigned int j=0;j<knri*4+10;j++){
            fitter.Config().ParSettings(j).SetLimits(parmsmin[j],parmsmax[j]);
            if (isparmsfix[j]){
               fitter.Config().ParSettings(j).Fix();
            }
        }

        fitter.FitFCN(knri*4+10,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),false);

        const Double_t* resultparmc=fitter.Result().GetParams();
        const Double_t* resulterrmc=fitter.Result().GetErrors();
        for (unsigned int i=0;i<knri*4+10;i++){
            outparms[i]=resultparmc[i];
            outparmserr[i]=resulterrmc[i];
        }
        treeout->Fill();
    }

    //! Filling the output tree
    for (unsigned int i=0;i<knri*4+10;i++){
      outparms[i]=resultpar[i];
      outparmserr[i]=resulterr[i];
    }
    treeout->Fill();



    cout<<"\n*******PRINTING RESULT**********\n"<<endl;
    result.Print(std::cout);

    cout<<"\n*******PRINTING MAIN RESULTS**********\n"<<endl;
    cout<<"t1/2\terr\tp1n\terr\tp2n\terr"<<endl;
    cout<<log(2)/resultpar[0]<<"\t"<<log(2)/resultpar[0]/resultpar[0]*resulterr[0]<<"\t"<<resultpar[knri]<<"\t"<<resulterr[knri]<<"\t"<<resultpar[knri*2]<<"\t"<<resulterr[knri*2]<<endl;


    cout<<"\n*****************\n"<<endl;
    cout<<log(2)/resultpar[0]<<"\t"<<log(2)/resultpar[0]/resultpar[0]*resulterr[0]<<"\t";
    cout<<resultpar[knri]<<"\t"<<resulterr[knri]<<"\t";
    cout<<resultpar[knri*2]<<"\t"<<resulterr[knri*2]<<endl;
    cout<<log(2)/resultpar[0]<<"\t"<<log(2)/resultpar[0]/resultpar[0]*resulterr[0]<<"\t";
    cout<<resultpar[knri]<<"\t"<<resulterr[knri]<<"\t";
    cout<<resultpar[knri*2]<<"\t"<<resulterr[knri*2]<<endl;


    //! write to text file
    std::ofstream ofs("fitresults.txt", std::ofstream::out | std::ofstream::app);
    ofs<<parmsfilename<<"\t"<<log(2)/resultpar[0]<<"\t"<<log(2)/resultpar[0]/resultpar[0]*resulterr[0]<<"\t"<<resultpar[knri]<<"\t"<<resulterr[knri]<<"\t"<<resultpar[knri*2]<<"\t"<<resulterr[knri*2]<<endl;


    TH1F* histcomphB=new TH1F("residual_decay","Residual of decay curve",hB->GetNbinsX(),hB->GetXaxis()->GetXmin(),hB->GetXaxis()->GetXmax());
    for (Int_t i=0;i<hB->GetNbinsX();i++){
        double res;
        if (hB->GetBinError(i+1)<0.001)
          res=  0;
        else
        res=  (hB->GetBinContent(i+1)- fB->Eval( hB->GetBinCenter(i+1) ) )/hB->GetBinError(i+1);

        histcomphB->SetBinContent(i+1,res);
        histcomphB->SetBinError(i+1,0);
    }

    TH1F* histcomphSB=new TH1F("residual_decay1neu","Residual of decay curve with one neutron gate",hSB->GetNbinsX(),hSB->GetXaxis()->GetXmin(),hSB->GetXaxis()->GetXmax());
    for (Int_t i=0;i<hSB->GetNbinsX();i++){
        double res;
        if (hSB->GetBinError(i+1)<0.001)
         res =  0;
        else
         res =  (hSB->GetBinContent(i+1)- fSB->Eval( hSB->GetBinCenter(i+1) ) )/hSB->GetBinError(i+1);
        histcomphSB->SetBinContent(i+1,res);
        histcomphSB->SetBinError(i+1,0);
    }

    TH1F* histcomphSB2=new TH1F("residual_decay2neu","Residual of decay curve with two neutron gate",hSB2->GetNbinsX(),hSB2->GetXaxis()->GetXmin(),hSB2->GetXaxis()->GetXmax());
    for (Int_t i=0;i<hSB2->GetNbinsX();i++){
        double res;
        if (hSB2->GetBinError(i+1)<0.001)
         res =  0;
        else
         res =  (hSB2->GetBinContent(i+1)- fSB2->Eval( hSB2->GetBinCenter(i+1) ) )/hSB2->GetBinError(i+1);
        histcomphSB2->SetBinContent(i+1,res);
        histcomphSB2->SetBinError(i+1,0);
    }


    TCanvas* c1=new TCanvas("c1","c1",900,1200);

    gStyle->SetOptStat(11111);
    c1->Divide(1,3);
    c1->cd(1);
    fB->SetFitResult( result, iparB);
    fB->SetRange(plotrange[0], plotrange[1]);
    fB->SetLineColor(kRed);
    fB->SetNpx(binning);

    hB->SetTitle("Decay curve");
    hB->GetListOfFunctions()->Add(fB);
    hB->GetXaxis()->SetRangeUser(plotrange[2],plotrange[1]);
    hB->GetXaxis()->SetTitle("Time (s)");
    hB->GetYaxis()->SetTitle("Counts");
    hB->GetXaxis()->SetTitleSize(0.05);
    hB->GetYaxis()->SetTitleSize(0.05);
    hB->SetMarkerStyle(20);
    hB->SetMarkerSize(0.8);
    hB->GetXaxis()->SetLabelSize(0.05);
    hB->GetYaxis()->SetLabelSize(0.05);
    hB->Draw("P0 E");

    hdecay->GetFunction("pol0")->SetLineColor(kYellow);
    hdecay->GetFunction("pol0")->Draw("same");

    c1->cd(2);


    fSB->SetFitResult( result, iparSB);
    fSB->SetRange(plotrange[0], plotrange[1]);
    fSB->SetLineColor(kRed);
    fSB->SetNpx(binning);

    hSB->SetTitle("Decay curve with one neutron gate");
    hSB->GetListOfFunctions()->Add(fSB);
    hSB->GetXaxis()->SetRangeUser(plotrange[2],plotrange[1]);
    hSB->GetXaxis()->SetTitle("Time (s)");
    hSB->GetYaxis()->SetTitle("Counts");
    hSB->GetXaxis()->SetTitleSize(0.05);
    hSB->GetYaxis()->SetTitleSize(0.05);
    hSB->SetMarkerStyle(20);
    hSB->SetMarkerSize(0.8);
    hSB->GetXaxis()->SetLabelSize(0.05);
    hSB->GetYaxis()->SetLabelSize(0.05);
    hSB->Draw("P0 E");

    hdecay1n->GetFunction("pol0")->SetLineColor(kYellow);
    hdecay1n->GetFunction("pol0")->Draw("same");

    c1->cd(3);
    fSB2->SetFitResult( result, iparSB2);
    fSB2->SetRange(plotrange[0], plotrange[1]);
    fSB2->SetLineColor(kRed);
    fSB2->SetNpx(binning);
    hSB2->GetXaxis()->SetTitle("Time (s)");
    hSB2->GetYaxis()->SetTitle("Counts");
    hSB2->GetXaxis()->SetTitleSize(0.05);
    hSB2->GetYaxis()->SetTitleSize(0.05);
    hSB2->GetXaxis()->SetLabelSize(0.05);
    hSB2->GetYaxis()->SetLabelSize(0.05);
    hSB2->SetTitle("Decay curve with two neutron gate");
    hSB2->GetListOfFunctions()->Add(fSB2);
    hSB2->GetXaxis()->SetRangeUser(plotrange[2],plotrange[1]);

    hSB2->SetMarkerStyle(20);
    hSB2->SetMarkerSize(0.8);
    hSB2->GetXaxis()->SetLabelSize(0.05);
    hSB2->GetYaxis()->SetLabelSize(0.05);
    hSB2->Draw("P0 E");

    hdecay2n->GetFunction("pol0")->SetLineColor(kYellow);
    hdecay2n->GetFunction("pol0")->Draw("same");

    TCanvas* c2=new TCanvas("c2","c2",900,900);
    c2->Divide(1,3);
    c2->cd(1);
    histcomphB->Draw();
    c2->cd(2);
    histcomphSB->Draw();
    c2->cd(3);
    histcomphSB2->Draw();

    hB->Write();
    fB->Write();
    histcomphB->Write();

    hSB->Write();
    fSB->Write();
    histcomphSB->Write();

    hSB2->Write();
    fSB2->Write();
    histcomphSB2->Write();


    c2->Write();
    c1->Write();


    outfile->Close();
}

