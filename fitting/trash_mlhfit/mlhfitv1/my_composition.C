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


using namespace RooFit ;



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
    //! parameters
    Double_t p_deadtime=-10.0;
    Double_t p_timerange=-0.0001;

    RooRealVar* x=new RooRealVar("x","x",p_deadtime,p_timerange) ;
    RooCategory* y=new RooCategory("y","y");
    y->defineType("0neu",0);
    y->defineType("1neu",1);
    y->defineType("2neu",2);

    fitFbkg model("bkgfit","bkgfit",*x,*y,*nbkg1,*nbkg2);

    RooExtendPdf modelext("bkgext","bkgext p.d.f",model,*nsig) ;

    // Sample 100000 events in (x,y) from the model
    //RooDataSet* modelData = model.generate(RooArgSet(x,y),10000) ;

    // I m p o r t i n g   i n t e g e r   T T r e e   b r a n c h e s
    // ---------------------------------------------------------------

    // Import integer tree branch as RooCategory
    // will be imported as those are the only defined states
    TFile *f=TFile::Open(infile);
    TTree* tree;
    f->GetObject("tree",tree);
    RooDataSet* data=new RooDataSet("data","data",RooArgSet(*y,*x),Import(*tree)) ;
    data->Print() ;

    // P e r f o r m   f i t s   i n   i n d i v i d u a l   s i d e b a n d   r e g i o n s
    // -------------------------------------------------------------------------------------


    modelext.fitTo(*data,Extended(),NumCPU(8),Save()) ;
    delete data;
    delete x;
    delete y;
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



void my_composition(char* infile,char* parmsfile,Int_t fitoption=0,Int_t isbkgfrombackwardtime=0, Double_t randcointratio=0)
{


    //! parameters
    Double_t p_deadtime=0.;
    Double_t p_timerange=10.;

    Double_t nsigma=10;


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


    std::ifstream ifs(parmsfile);
    Int_t rino;
    Double_t temp;
    Bool_t flagfix[knri][3];
    Double_t decayparms[knri][3];
    Double_t decayparms_p[knri][3];
    Double_t decayparms_m[knri][3];
    Double_t decayparms_err[knri][3];


    for (int i=0;i<knri;i++){
        ifs>>rino;
        for (int j=0;j<3;j++){
            ifs>>temp;
            if (temp>=0){
                flagfix[i][j]=true;
                if (j==0) {
                    decayparms[i][j]=temp;
                    ifs>>temp;decayparms_m[i][j]=temp;decayparms_p[i][j]=temp;
                }else{
                    decayparms[i][j]=temp;
                    ifs>>temp;decayparms_p[i][j]=temp;decayparms_m[i][j]=temp;
                }
            }else{
                flagfix[i][j]=false;
                if (j==0) {
                    decayparms[i][j]=(-temp);
                    ifs>>temp;decayparms_p[i][j]=temp;
                    decayparms_m[i][j]=temp;
                }else{
                    decayparms[i][j]=-temp;
                    ifs>>temp;decayparms_p[i][j]=1;
                    decayparms_m[i][j]=0;
                }

            }
        }
    }
    Double_t bkg,bkg_p,bkg_m,init,init_p,init_m;
    Double_t bkg2,bkg2_p,bkg2_m;

    ifs>>init>>init_p>>init_m>>bkg>>bkg_m>>bkg_p>>bkg2>>bkg2_m>>bkg2_p;
    for (int i=0;i<knri;i++){
        cout<<"ri="<<i+1<<" : ";
        for (int j=0;j<3;j++){
            cout<<flagfix[i][j]<<"\t"<<decayparms[i][j]<<"\t"<<decayparms_p[i][j]<<"\t"<<decayparms_m[i][j]<<"\t";
            if (j==0){//for half-life
                decayparms_err[i][j]=log(2)/decayparms[i][j]/decayparms[i][j]*decayparms_m[i][j];
                decayparms_p[i][j]=log(2)/decayparms[i][j]+log(2)/decayparms[i][j]/decayparms[i][j]*decayparms_m[i][j];
                decayparms_m[i][j]=log(2)/decayparms[i][j]-log(2)/decayparms[i][j]/decayparms[i][j]*decayparms_m[i][j];
                decayparms[i][j]=log(2)/decayparms[i][j];
                cout<<decayparms[i][j]<<"-"<<decayparms_m[i][j]<<"-"<<decayparms_p[i][j]<<"-"<<decayparms[i][j]-decayparms_m[i][j]<<"-"<<decayparms_p[i][j]-decayparms[i][j]<<endl;
            }
        }
        cout<<endl;
    }
    cout<<"\t"<<bkg<<"\t"<<bkg_m<<"\t"<<bkg_p<<"\t"<<bkg2<<"\t"<<bkg2_m<<"\t"<<bkg2_p<<"\t"<<init<<"\t"<<init_p<<"\t"<<init_m<<endl;

    //! Fitting positive range
    RooRealVar x("x","x",p_deadtime,p_timerange) ;
    RooRealVar p0("p0","p0",log(2)/0.1,log(2)/1000000000000.,log(2)/0.00001);//l1
    RooRealVar p1("p1","p1",log(2)/0.1,log(2)/1000000000000.,log(2)/0.00001);//l2
    RooRealVar p2("p2","p2",log(2)/0.1,log(2)/1000000000000.,log(2)/0.00001);//l3
    RooRealVar p3("p3","p3",log(2)/0.1,log(2)/1000000000000.,log(2)/0.00001);//l4
    RooRealVar p4("p4","p4",log(2)/0.1,log(2)/1000000000000.,log(2)/0.00001);//l5
    RooRealVar p5("p5","p5",log(2)/0.1,log(2)/1000000000000.,log(2)/0.00001);//l7
    RooRealVar p6("p6","p6",log(2)/0.1,log(2)/1000000000000.,log(2)/0.00001);//l7
    RooRealVar p7("p7","p7",log(2)/0.1,log(2)/1000000000000.,log(2)/0.00001);//l8
    RooRealVar p8("p8","p8",log(2)/0.1,log(2)/1000000000000.,log(2)/0.00001);//l9
    RooRealVar p9("p9","p9",0.5,0,1);//pn1
    RooRealVar p10("p10","p10",0.5,0,1);//pn2
    RooRealVar p11("p11","p11",0.5,0,1);//pn3
    RooRealVar p12("p12","p12",0.5,0,1);//pn4
    RooRealVar p13("p13","p13",0.5,0,1);//pn5
    RooRealVar p14("p14","p14",0.5,0,1);//pn6
    RooRealVar p15("p15","p15",0.5,0,1);//pn7
    RooRealVar p16("p16","p16",0.5,0,1);//pn8
    RooRealVar p17("p17","p17",0.5,0,1);//pn9
    RooRealVar p18("p18","p18",0.5,0,1);//p2n1
    RooRealVar p19("p19","p19",100,0,1000000);//activity
    RooRealVar p20("p20","p20",100,0,1000000);//background
    RooRealVar p21("p21","p21",100,0,1000000);//background 2
    RooRealVar p22("p22","p22",100,0,1000000);//background 3
    RooRealVar p23("p23","p23",0.,0.,1.);//random coincidence

    //! Constrains for error propagation
    RooGaussian* pconstr[19];
    RooGaussian* pconstrbkg[3];

    RooArgSet constronlydecay;
    RooArgSet constronlydecaywbkg;

    for (int i=0;i<knri;i++){
        for (int j=0;j<2;j++){
            if (flagfix[i][j]){//fix value
                if (j*knri+i==0){
                    p0.setVal(decayparms[i][j]);
                    if (decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma>0) p0.setMin(decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma);
                    else p0.setMin(log(2)/1000000000000);
                    p0.setMax(decayparms[i][j]+(decayparms_p[i][j]-decayparms[i][j])*nsigma);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p0.setConstant(kTRUE);

                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[0]=new RooGaussian("p0constr","p0constr",p0,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[0]=new RooGaussian("p0constr","p0constr",p0,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        cout<<"set constrain model for p0"<<endl;
                        constronlydecay.add(*pconstr[0]);
                        constronlydecaywbkg.add(*pconstr[0]);
                    }

                }else if (j*knri+i==1){
                    p1.setVal(decayparms[i][j]);
                    if (decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma>0) p1.setMin(decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma);
                    else p1.setMin(log(2)/1000000000000);
                    p1.setMax(decayparms[i][j]+(decayparms_p[i][j]-decayparms[i][j])*nsigma);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p1.setConstant(kTRUE);

                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[1]=new RooGaussian("p1constr","p1constr",p1,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[1]=new RooGaussian("p1constr","p1constr",p1,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        cout<<"set constrain model for p1 "<<decayparms[i][j]<<"-"<<decayparms[i][j]-decayparms_m[i][j]<<"-"<<decayparms[i][j]-decayparms_m[i][j]<<endl;
                        constronlydecay.add(*pconstr[1]);
                        constronlydecaywbkg.add(*pconstr[1]);
                    }

                }else if (j*knri+i==2){
                    p2.setVal(decayparms[i][j]);
                    if (decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma>0) p2.setMin(decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma);
                    else p2.setMin(log(2)/1000000000000);
                    p2.setMax(decayparms[i][j]+(decayparms_p[i][j]-decayparms[i][j])*nsigma);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p2.setConstant(kTRUE);

                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[2]=new RooGaussian("p2constr","p2constr",p2,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[2]=new RooGaussian("p2constr","p2constr",p2,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        cout<<"set constrain model for p2"<<endl;
                        constronlydecay.add(*pconstr[2]);
                        constronlydecaywbkg.add(*pconstr[2]);
                    }
                }else if (j*knri+i==3){
                    p3.setVal(decayparms[i][j]);
                    if (decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma>0) p3.setMin(decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma);
                    else p3.setMin(log(2)/1000000000000);
                    p3.setMax(decayparms[i][j]+(decayparms_p[i][j]-decayparms[i][j])*nsigma);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p3.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[3]=new RooGaussian("p3constr","p3constr",p3,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[3]=new RooGaussian("p3constr","p3constr",p3,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        cout<<"set constrain model for p3"<<endl;
                        constronlydecay.add(*pconstr[3]);
                        constronlydecaywbkg.add(*pconstr[3]);
                    }
                }else if (j*knri+i==4){
                    p4.setVal(decayparms[i][j]);
                    if (decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma>0) p4.setMin(decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma);
                    else p4.setMin(log(2)/1000000000000);
                    p4.setMax(decayparms[i][j]+(decayparms_p[i][j]-decayparms[i][j])*nsigma);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p4.setConstant(kTRUE);

                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[4]=new RooGaussian("p4constr","p4constr",p4,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[4]=new RooGaussian("p4constr","p4constr",p4,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        cout<<"set constrain model for p4"<<endl;
                        constronlydecay.add(*pconstr[4]);
                        constronlydecaywbkg.add(*pconstr[4]);
                    }

                }else if (j*knri+i==5){
                    p5.setVal(decayparms[i][j]);
                    if (decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma>0) p5.setMin(decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma);
                    else p5.setMin(log(2)/1000000000000);
                    p5.setMax(decayparms[i][j]+(decayparms_p[i][j]-decayparms[i][j])*nsigma);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p5.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[5]=new RooGaussian("p5constr","p5constr",p5,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[5]=new RooGaussian("p5constr","p5constr",p5,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        cout<<"set constrain model for p5"<<endl;
                        constronlydecay.add(*pconstr[5]);
                        constronlydecaywbkg.add(*pconstr[5]);
                    }

                }else if (j*knri+i==6){
                    p6.setVal(decayparms[i][j]);
                    if (decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma>0) p6.setMin(decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma);
                    else p6.setMin(log(2)/1000000000000);
                    p6.setMax(decayparms[i][j]+(decayparms_p[i][j]-decayparms[i][j])*nsigma);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p6.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[6]=new RooGaussian("p6constr","p6constr",p6,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[6]=new RooGaussian("p6constr","p6constr",p6,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        cout<<"set constrain model for p6"<<endl;
                        constronlydecay.add(*pconstr[6]);
                        constronlydecaywbkg.add(*pconstr[6]);
                    }
                }else if (j*knri+i==7){
                    p7.setVal(decayparms[i][j]);
                    if (decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma>0) p7.setMin(decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma);
                    else p7.setMin(log(2)/1000000000000);
                    p7.setMax(decayparms[i][j]+(decayparms_p[i][j]-decayparms[i][j])*nsigma);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p7.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[7]=new RooGaussian("p7constr","p7constr",p7,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[7]=new RooGaussian("p7constr","p7constr",p7,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        cout<<"set constrain model for p7"<<endl;
                        constronlydecay.add(*pconstr[7]);
                        constronlydecaywbkg.add(*pconstr[7]);
                    }
                }else if (j*knri+i==8){
                    p8.setVal(decayparms[i][j]);
                    if (decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma>0) p8.setMin(decayparms[i][j]-(decayparms[i][j]-decayparms_m[i][j])*nsigma);
                    else p8.setMin(log(2)/1000000000000);
                    p8.setMax(decayparms[i][j]+(decayparms_p[i][j]-decayparms[i][j])*nsigma);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p8.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[8]=new RooGaussian("p8constr","p8constr",p8,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[8]=new RooGaussian("p8constr","p8constr",p8,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        cout<<"set constrain model for p8"<<endl;
                        constronlydecay.add(*pconstr[8]);
                        constronlydecaywbkg.add(*pconstr[8]);
                    }
                }else if (j*knri+i==9){
                    p9.setVal(decayparms[i][j]);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p9.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[9]=new RooGaussian("p9constr","p9constr",p9,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[9]=new RooGaussian("p9constr","p9constr",p9,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        cout<<"set constrain model for p9"<<endl;
                        constronlydecay.add(*pconstr[9]);
                        constronlydecaywbkg.add(*pconstr[9]);
                    }
                }else if (j*knri+i==10){
                    p10.setVal(decayparms[i][j]);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p10.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[10]=new RooGaussian("p10constr","p10constr",p10,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[10]=new RooGaussian("p10constr","p10constr",p10,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        cout<<"set constrain model for p10"<<endl;
                        constronlydecay.add(*pconstr[10]);
                        constronlydecaywbkg.add(*pconstr[10]);
                    }
                }else if (j*knri+i==11){
                    p11.setVal(decayparms[i][j]);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p11.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[11]=new RooGaussian("p11constr","p11constr",p11,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[11]=new RooGaussian("p11constr","p11constr",p11,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        cout<<"set constrain model for p11"<<endl;
                        constronlydecay.add(*pconstr[11]);
                        constronlydecaywbkg.add(*pconstr[11]);
                    }
                }else if (j*knri+i==12){
                    p12.setVal(decayparms[i][j]);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p12.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[12]=new RooGaussian("p12constr","p12constr",p12,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[12]=new RooGaussian("p12constr","p12constr",p12,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        constronlydecay.add(*pconstr[12]);
                        constronlydecaywbkg.add(*pconstr[12]);
                    }
                }else if (j*knri+i==13){
                    p13.setVal(decayparms[i][j]);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p13.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[13]=new RooGaussian("p13constr","p13constr",p13,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[13]=new RooGaussian("p13constr","p13constr",p13,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        constronlydecay.add(*pconstr[13]);
                        constronlydecaywbkg.add(*pconstr[13]);
                    }
                }else if (j*knri+i==14){
                    p14.setVal(decayparms[i][j]);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p14.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[14]=new RooGaussian("p14constr","p14constr",p14,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[14]=new RooGaussian("p14constr","p14constr",p14,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        constronlydecay.add(*pconstr[14]);
                        constronlydecaywbkg.add(*pconstr[14]);
                    }
                }else if (j*knri+i==15){
                    p15.setVal(decayparms[i][j]);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p15.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[15]=new RooGaussian("p15constr","p15constr",p15,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[15]=new RooGaussian("p15constr","p15constr",p15,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        constronlydecay.add(*pconstr[15]);
                        constronlydecaywbkg.add(*pconstr[15]);
                    }
                }else if (j*knri+i==16){
                    p16.setVal(decayparms[i][j]);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p16.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[16]=new RooGaussian("p16constr","p16constr",p16,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[16]=new RooGaussian("p16constr","p16constr",p16,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        constronlydecay.add(*pconstr[16]);
                        constronlydecaywbkg.add(*pconstr[16]);
                    }
                }else if (j*knri+i==17){
                    p17.setVal(decayparms[i][j]);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p17.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[17]=new RooGaussian("p17constr","p17constr",p17,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[17]=new RooGaussian("p17constr","p17constr",p17,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        constronlydecay.add(*pconstr[17]);
                        constronlydecaywbkg.add(*pconstr[17]);
                    }
                }else if (j*knri+i==18){
                    p18.setVal(decayparms[i][j]);
                    if ((fitoption!=1&&fitoption!=2)||(decayparms_m[i][j]==0||decayparms_m[i][j]==0))
                        p17.setConstant(kTRUE);
                    if (decayparms[i][j]-decayparms_m[i][j]>decayparms_p[i][j]-decayparms[i][j])
                        pconstr[18]=new RooGaussian("p18constr","p18constr",p18,RooConst(decayparms[i][j]),RooConst(decayparms[i][j]-decayparms_m[i][j]));
                    else
                        pconstr[18]=new RooGaussian("p18constr","p18constr",p18,RooConst(decayparms[i][j]),RooConst(decayparms_p[i][j]-decayparms[i][j]));
                    if (decayparms_m[i][j]!=0&&decayparms_m[i][j]!=0){
                        constronlydecay.add(*pconstr[18]);
                        constronlydecaywbkg.add(*pconstr[18]);
                    }
                }
                cout<<"Set fix parameter "<<j*knri+i<<" with val = "<<decayparms[i][j]<<endl;

            }else{//variable parameters
                if (j*knri+i==0){
                    p0.setVal(decayparms[i][j]);
                    p0.setMin(decayparms_m[i][j]);
                    p0.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==1){
                    p1.setVal(decayparms[i][j]);
                    p1.setMin(decayparms_m[i][j]);
                    p1.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==2){
                    p2.setVal(decayparms[i][j]);
                    p2.setMin(decayparms_m[i][j]);
                    p2.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==3){
                    p3.setVal(decayparms[i][j]);
                    p3.setMin(decayparms_m[i][j]);
                    p3.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==4){
                    p4.setVal(decayparms[i][j]);
                    p4.setMin(decayparms_m[i][j]);
                    p4.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==5){
                    p5.setVal(decayparms[i][j]);
                    p5.setMin(decayparms_m[i][j]);
                    p5.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==6){
                    p6.setVal(decayparms[i][j]);
                    p6.setMin(decayparms_m[i][j]);
                    p6.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==7){
                    p7.setVal(decayparms[i][j]);
                    p7.setMin(decayparms_m[i][j]);
                    p7.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==8){
                    p8.setVal(decayparms[i][j]);
                    p8.setMin(decayparms_m[i][j]);
                    p8.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==9){
                    p9.setVal(decayparms[i][j]);
                    p9.setMin(decayparms_m[i][j]);
                    p9.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==10){
                    p10.setVal(decayparms[i][j]);
                    p10.setMin(decayparms_m[i][j]);
                    p10.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==11){
                    p11.setVal(decayparms[i][j]);
                    p11.setMin(decayparms_m[i][j]);
                    p11.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==12){
                    p12.setVal(decayparms[i][j]);
                    p12.setMin(decayparms_m[i][j]);
                    p12.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==13){
                    p13.setVal(decayparms[i][j]);
                    p13.setMin(decayparms_m[i][j]);
                    p13.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==14){
                    p14.setVal(decayparms[i][j]);
                    p14.setMin(decayparms_m[i][j]);
                    p14.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==15){
                    p15.setVal(decayparms[i][j]);
                    p15.setMin(decayparms_m[i][j]);
                    p15.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==16){
                    p16.setVal(decayparms[i][j]);
                    p16.setMin(decayparms_m[i][j]);
                    p16.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==17){
                    p17.setVal(decayparms[i][j]);
                    p17.setMin(decayparms_m[i][j]);
                    p17.setMax(decayparms_p[i][j]);
                }else if (j*knri+i==18){
                    p18.setVal(decayparms[i][j]);
                    p18.setMin(decayparms_m[i][j]);
                    p18.setMax(decayparms_p[i][j]);
                }
                cout<<"Set valriable parameter "<<j*knri+i<<" with val = "<<decayparms[i][j]<<" varying from "<<decayparms_m[i][j]<<" to "<<decayparms_p[i][j]<<endl;
            }
        }
    }

    if (flagfix[0][2]){
        p18.setVal(decayparms[0][2]);
        if ((fitoption!=1&&fitoption!=2)||(decayparms_m[0][2]==0||decayparms_m[0][2]==0))
            p18.setConstant(kTRUE);
        cout<<"Set fix parameter 18(p2n) with val = "<<decayparms[0][2]<<endl;
    }else{
        p18.setVal(decayparms[0][2]);
        p18.setMin(decayparms_m[0][2]);
        p18.setMax(decayparms_p[0][2]);
        cout<<"Set valriable parameter 18(p2n) with val = "<<decayparms[0][2]<<" varying from "<<decayparms_m[0][2]<<" to "<<decayparms_p[0][2]<<endl;
    }

    //! set inital activity
    p19.setVal(init);
    p19.setMin(init_p);
    p19.setMax(init_m);

    p20.setVal(bkg);
    p20.setMin(bkg_m);
    p20.setMax(bkg_p);
    p21.setVal(bkg2);
    p21.setMin(bkg2_m);
    p21.setMax(bkg2_p);
    p22.setVal(bkg2);
    p22.setMin(bkg2_m);
    p22.setMax(bkg2_p);

    if (randcointratio>=0) {
        p23.setVal(randcointratio);
        p23.setConstant(kTRUE);
    } else {
        p23.setVal(0.5);
        p23.setMax(1.);
        p23.setMin(0.);
    }

    //! neutron detection efficinecy
   RooRealVar neueff("neueff","neueff",0.613384,0.,1.) ;
   neueff.setConstant(kTRUE);

   // C r e a t e   2 D   p d f   a n d   d a t a
   // -------------------------------------------

   // Define discrete variable y
   RooCategory y("y","y");
   y.defineType("0neu",0);
   y.defineType("1neu",1);
   y.defineType("2neu",2);



   //! fix val for extended model
   p19.setMin(0);
   p19.setMax(2);
   p19.setVal(1);
   p20.setMin(-1);
   p20.setMax(1);
   p20.setVal(0);
   p21.setMin(-1);
   p21.setMax(1);
   p21.setVal(0);
   p22.setMin(-1);
   p22.setMax(1);
   p22.setVal(0);
   p19.setConstant(kTRUE);
   p20.setConstant(kTRUE);
   p21.setConstant(kTRUE);
   p22.setConstant(kTRUE);

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


   ntotalbkg.setConstant(kTRUE);
   bkgratio1.setConstant(kTRUE);
   bkgratio2.setConstant(kTRUE);

   // pdf model
   fitF model("signal","signal",x,y,neueff,*p[0],*p[1],*p[2],*p[3],*p[4],*p[5],*p[6],*p[7],*p[8],*p[9],*p[10],*p[11],*p[12],*p[13],*p[14],*p[15],*p[16],*p[17],*p[18],p19,p20,p21,p22,p23);

   //fitF model("signal","signal",x,y,neueff,*p[0],*p[1],*p[2],*p[3],*p[4],*p[5],*p[6],*p[7],*p[8],*p[9],*p[10],*p[11],*p[12],*p[13],*p[14],*p[15],*p[16],*p[17],*p[18],*p[19],*p[20],*p[21],*p[22],*p[23]);

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
   RooDataSet* data=new RooDataSet("data","data",RooArgSet(y,x),Import(*tree)) ;
   data->Print() ;

   // P e r f o r m   f i t s   i n   i n d i v i d u a l   s i d e b a n d   r e g i o n s
   // -------------------------------------------------------------------------------------

   // Perform fit in SideBand1 region (RooAddPdf coefficients will be interpreted in full range)
   //model.fitTo(*modelData) ;

   if (fitoption==1)
       model.fitTo(*data,ExternalConstraints(constronlydecay)) ;
   else if (fitoption==2)
       model.fitTo(*data,ExternalConstraints(constronlydecaywbkg)) ;
   else
       final_pdf.fitTo(*data,NumCPU(16)) ;


    // Print results for comparison
   //r_sb1->Print() ;

   // Plot x distribution of data and projection of model on x = Int(dy) model(x,y)
   RooPlot* xframe0 = x.frame(Title("0 neutron fit")) ;
   //modelData->plotOn(xframe0,Cut("y==y::0neu")) ;
   data->plotOn(xframe0,Cut("y==y::0neu")) ;
   final_pdf.plotOn(xframe0,Slice(y,"0neu")) ;

   RooPlot* xframe1 = x.frame(Title("1 neutron fit")) ;
   //modelData->plotOn(xframe1,Cut("y==y::1neu")) ;
   data->plotOn(xframe1,Cut("y==y::1neu")) ;
   final_pdf.plotOn(xframe1,Slice(y,"1neu")) ;

   RooPlot* xframe2 = x.frame(Title("2 neutrons fit")) ;
   //modelData->plotOn(xframe2,Cut("y==y::2neu")) ;
   data->plotOn(xframe2,Cut("y==y::2neu")) ;
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

   // Make canvas and draw RooPlots
   TCanvas *c = new TCanvas("comp","comp",1200, 800);
   c->Divide(2,3);
   c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe0->GetYaxis()->SetTitleOffset(1.4) ; xframe0->Draw() ;
   c->cd(3) ; gPad->SetLeftMargin(0.15) ; xframe1->GetYaxis()->SetTitleOffset(1.4) ; xframe1->Draw() ;
   c->cd(5) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.4) ; xframe2->Draw() ;
   c->cd(2) ; gPad->SetLeftMargin(0.15) ; rframe0->GetYaxis()->SetTitleOffset(1.4) ; rframe0->Draw() ;
   c->cd(4) ; gPad->SetLeftMargin(0.15) ; rframe1->GetYaxis()->SetTitleOffset(1.4) ; rframe1->Draw() ;
   c->cd(6) ; gPad->SetLeftMargin(0.15) ; rframe2->GetYaxis()->SetTitleOffset(1.4) ; rframe2->Draw() ;

   cout<<"\n********************SUMARRY***************"<<endl;
   cout<<"T1/2= "<<log(2)/p0.getVal()<<" +/- "<<log(2)/p0.getVal()/p0.getVal()*p0.getError()<<endl;
   cout<<"P1n= "<<p9.getVal()<<" +/- "<<p9.getError()<<endl;
   cout<<"P2n= "<<p18.getVal()<<" +/- "<<p18.getError()<<endl;
   cout<<endl;


   cout<<log(2)/p0.getVal()<<"\t"<<log(2)/p0.getVal()/p0.getVal()*p0.getError()<<endl;
   cout<<p9.getVal()<<"\t"<<p9.getError()<<endl;
   cout<<p18.getVal()<<"\t"<<p18.getError()<<endl;

}


