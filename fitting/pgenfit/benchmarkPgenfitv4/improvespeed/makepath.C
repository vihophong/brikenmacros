#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <unistd.h>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"

#include "TLine.h"
#include "TArrow.h"


#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/PoissonLikelihoodFCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"

#include "TTree.h"
#include "TRandom3.h"
#include "TFile.h"
#include <fstream>


#include "TCutG.h"
#include "TLatex.h"
#include "TMath.h"
#include "TCanvas.h"

#include "TGraph.h"
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/GSLSimAnMinimizer.h"
#include "Math/Functor.h"

#include <list>

#include <sstream>
#include <string>

using namespace std;



//Global variable for core function
const Int_t kmaxndecay=200;
const Int_t kmaxpaths=200;
const Int_t kmaxparms=150;

Int_t npaths=100;
Int_t knri=200;

Int_t ndecay[kmaxpaths];
Int_t decaymap[kmaxpaths][kmaxndecay];
Int_t nneu[kmaxpaths][kmaxndecay];

TString riname[kmaxparms];
Double_t parms[kmaxparms];
Double_t parmserr[kmaxparms];
Double_t parmsmax[kmaxparms];
Double_t parmsmin[kmaxparms];
Int_t isparmsfix[kmaxparms];
Double_t mcparms[kmaxparms];

Int_t nisomers = 0;
Int_t groundstate[kmaxparms];
Int_t isomerstate[kmaxparms];

Int_t nparmsactive = 0;
Int_t indexparmsactive[kmaxparms];

void plainPlot(TCanvas* c1,Double_t xrange[],Double_t yrange[])
{

  Double_t minhalflife=0.0001;//100 ns

  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(12);
  gStyle->SetOptStat(0);

  Int_t NumRI;

  Int_t nprot;
  Int_t nneut;
  Int_t nmass;
  Double_t hlval;

  //TCanvas* c1=new TCanvas("c1","",900,700) ;


  NumRI = 5346;
  //  NumRI = 24;

  //  ofstream fout2;
  //  fout2.open("zzz.dat");

  ifstream fdat;
  fdat.open("FRDM-QRPA12-halflife.txt");

  cout << "Get data" << endl;

  TH2F *hchart = new TH2F("hist","",185,-0.5,184.5,127,-0.5,126.5);


  for (Int_t i=0; i<NumRI; i++) {
    fdat >> nprot >> nneut >> hlval;
    nmass = nneut+nprot;
    hchart->Fill(nneut,nprot,hlval);
  }

  c1->SetLogz(0);
  hchart->SetTitleSize(0.04);
  hchart->GetXaxis()->SetTitleOffset(1.0);
  hchart->GetYaxis()->SetTitleOffset(1.2);
  hchart->GetYaxis()->CenterTitle();
  hchart->GetXaxis()->SetLabelSize(0.03);
  hchart->GetYaxis()->SetLabelSize(0.03);
  //  hchart->GetXaxis()->SetTitleSize(1.1);
  hchart->GetYaxis()->SetTitle("N_{Proton}");
  hchart->GetXaxis()->SetTitle("N_{Neutron}");

  hchart->GetXaxis()->SetRangeUser(xrange[0],xrange[1]);
  hchart->GetYaxis()->SetRangeUser(yrange[0],yrange[1]);

  hchart->SetMinimum(minhalflife);


  c1->SetLogz();
  hchart->SetLineWidth(10);
  hchart->SetLineColor(1);
  hchart->Draw("COLZ");



  //! draw border of isotopes
  TBox b2;
  b2.SetFillStyle(0);
  b2.SetLineColor(2);
  b2.SetLineWidth(1);
  fdat.seekg(0, ios::beg);
  for (Int_t i=0; i<NumRI; i++) {
    fdat >> nprot >> nneut >> hlval;
    if(nprot>=yrange[0]&&nprot<=yrange[1]&&nneut>=xrange[0]&&nneut<=xrange[1]) b2.DrawBox(nneut-0.5,nprot-0.5,nneut+0.5,nprot+0.5);
  }
  fdat.close();



  //! Drawing magic number
  Double_t dd = 0.5;
  TLine a1;
  //  a1.SetLineWidth(1.5);
  a1.SetLineWidth(3.0);
  a1.SetLineColor(7);

  Int_t magicn[]={8,20,28,50,82,126};

  for (Int_t i=0;i<6;i++){
      a1.DrawLine(magicn[i]-dd,yrange[0]-dd,magicn[i]-dd,yrange[1]+dd); a1.DrawLine(magicn[i]+1-dd,yrange[0]-dd,magicn[i]+1-dd,yrange[1]+dd);
      a1.DrawLine(xrange[0]-dd,magicn[i]-dd,xrange[1]+dd,magicn[i]-dd); a1.DrawLine(xrange[0]-dd,magicn[i]+1-dd,xrange[1]+dd,magicn[i]+1-dd);
  }


  //! Drawing r-process path
  TLine a0;
  a0.SetLineWidth(4);
  a0.SetLineStyle(1);
  a0.SetLineColor(3);
  Double_t nn;
  Double_t pp1;
  Double_t pp2;
  ifstream rpathfile("r-process_path.txt");
   while (rpathfile.good()){
       rpathfile>>nn>>pp1>>pp2;
       if (nn>=xrange[0]&&nn<=xrange[1]){
           Bool_t isplot=true;
           if (pp1<yrange[0]) pp1=yrange[0];
           if (pp1>yrange[1]) isplot=false;
           if (pp2<yrange[0]) isplot=false;
           if (pp2>yrange[1]) pp2=yrange[1];
           if (isplot){
            a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5);a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5);
            a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5);a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);
           }
       }
   }

  ifstream fdat7;
  fdat7.open("stable.csv");
  cout << "Stable data " << endl;
  Int_t n7=287;
  Double_t xx7[300];
  Double_t yy7[300];

  for (Int_t i=0; i<n7; i++) {
    fdat7 >> yy7[i] >> xx7[i];
  }
  TGraph *gr7 = new TGraph(n7,xx7,yy7);
  gr7->SetMarkerStyle(21);
  gr7->SetMarkerColor(1);
  gr7->SetMarkerSize(1.8);
  gr7->Draw("PS");
}


typedef struct {
    // idendification
    Int_t id;
    Int_t z;
    Int_t n;
    TString name;

    // decay properies
    Double_t decay_hl;
    Double_t decay_lamda;

    Double_t population_ratio;//isomer population ratio
    Double_t population_ratioerr;
    Double_t population_ratioup;
    Double_t population_ratiolow;

    Double_t decay_p0n;//decay to several isomerics states or ground state
    Double_t decay_p1n;
    Double_t decay_p2n;

    Double_t decay_hlerr;
    Double_t decay_lamdaerr;
    Double_t decay_p0nerr;
    Double_t decay_p1nerr;
    Double_t decay_p2nerr;

    Double_t decay_hlup;
    Double_t decay_lamdaup;
    Double_t decay_p0nup;
    Double_t decay_p1nup;
    Double_t decay_p2nup;

    Double_t decay_hllow;
    Double_t decay_lamdalow;
    Double_t decay_p0nlow;
    Double_t decay_p1nlow;
    Double_t decay_p2nlow;

    // fit options 0-vary, 1-fix with error propagation, 2-fix without error propagation
    Int_t is_decay_hl_fix;
    Int_t is_decay_lamda_fix;
    Int_t is_decay_p0n_fix;
    Int_t is_decay_p1n_fix;
    Int_t is_decay_p2n_fix;

    Int_t is_population_ratio_fix;

    // paths to this ri
    vector< vector<Int_t> > path;
    vector< vector<Int_t> > nneupath;
    vector< vector<Double_t> > branching;
} MemberDef;

void ProcessMember(MemberDef* obj)
{
    //! further processing
    if (obj->decay_hl<0){
        obj->decay_hl=-obj->decay_hl;
        obj->is_decay_hl_fix=0;
    }else{//exclude decay half-life =0;
        obj->is_decay_hl_fix=1;
    }


    // convert half-life into activity
    obj->decay_lamda=log(2)/obj->decay_hl;
    obj->decay_lamdaerr=log(2)/obj->decay_hl/obj->decay_hl*obj->decay_hlerr;
    obj->decay_lamdalow=log(2)/obj->decay_hlup;
    obj->decay_lamdaup=log(2)/obj->decay_hllow;
    obj->is_decay_lamda_fix=obj->is_decay_hl_fix;

    if (obj->decay_p1n<0){
        obj->decay_p1n=-obj->decay_p1n;
        obj->is_decay_p1n_fix=0;
    }else if (obj->decay_p1n==0){//exclude decay with pn = 0;
        obj->is_decay_p1n_fix=2;
    }else{
        obj->is_decay_p1n_fix=1;
    }

    if (obj->decay_p2n<0){
        obj->decay_p2n=-obj->decay_p2n;
        obj->is_decay_p2n_fix=0;
    }else if (obj->decay_p2n==0){//exclude decay half-life =0;
        obj->is_decay_p2n_fix=2;
    }else{
        obj->is_decay_p2n_fix=1;
    }

    // convert pn in % to pn in 1
    obj->decay_p1n=obj->decay_p1n/100;
    obj->decay_p1nerr=obj->decay_p1nerr/100;
    obj->decay_p1nlow=obj->decay_p1nlow/100;
    obj->decay_p1nup=obj->decay_p1nup/100;

    obj->decay_p2n=obj->decay_p2n/100;
    obj->decay_p2nerr=obj->decay_p2nerr/100;
    obj->decay_p2nlow=obj->decay_p2nlow/100;
    obj->decay_p2nup=obj->decay_p2nup/100;
}

void CopyMember(MemberDef* source, MemberDef* destination)
{
    destination-> id = source->  id;
    destination-> z = source->  z;
    destination-> n = source->  n;
    destination-> name = source->  name;

    destination-> decay_hl = source->  decay_hl;
    destination-> decay_lamda = source->  decay_lamda;

    destination-> population_ratio = source->  population_ratio;
    destination-> population_ratioerr = source->  population_ratioerr;
    destination-> population_ratioup = source->  population_ratioup;
    destination-> population_ratiolow = source->  population_ratiolow;

    destination-> decay_p0n = source->  decay_p0n;
    destination-> decay_p1n = source->  decay_p1n;
    destination-> decay_p2n = source->  decay_p2n;

    destination-> decay_hlerr = source->  decay_hlerr;
    destination-> decay_lamdaerr = source->  decay_lamdaerr;
    destination-> decay_p0nerr = source->  decay_p0nerr;
    destination-> decay_p1nerr = source->  decay_p1nerr;
    destination-> decay_p2nerr = source->  decay_p2nerr;

    destination-> decay_hlup = source->  decay_hlup;
    destination-> decay_lamdaup = source->  decay_lamdaup;
    destination-> decay_p0nup = source->  decay_p0nup;
    destination-> decay_p1nup = source->  decay_p1nup;
    destination-> decay_p2nup = source->  decay_p2nup;

    destination-> decay_hllow = source->  decay_hllow;
    destination-> decay_lamdalow = source->  decay_lamdalow;
    destination-> decay_p0nlow = source->  decay_p0nlow;
    destination-> decay_p1nlow = source->  decay_p1nlow;
    destination-> decay_p2nlow = source->  decay_p2nlow;


    destination-> is_decay_hl_fix = source->  is_decay_hl_fix;
    destination-> is_decay_lamda_fix = source->  is_decay_lamda_fix;
    destination-> is_decay_p0n_fix = source->  is_decay_p0n_fix;
    destination-> is_decay_p1n_fix = source->  is_decay_p1n_fix;
    destination-> is_decay_p2n_fix = source->  is_decay_p2n_fix;
}


void readinput(char* infilename, list<MemberDef*>& list)
{
    Int_t id=0;
    nisomers = 0;
    std::string line;
    std::ifstream infile(infilename);
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (line[0]=='#') continue;
        Int_t a;
        // decay properies
        MemberDef* obj=new MemberDef();
        obj->id=id;

        //! read info without isomer
//        if (!(iss >> obj->name >> obj->z >> a >> obj->decay_hl >> obj->decay_hlerr >> obj->decay_hllow >> obj->decay_hlup >>
//              obj->decay_p1n >> obj->decay_p1nerr >> obj->decay_p1nlow >> obj->decay_p1nup >>
//              obj->decay_p2n >> obj->decay_p2nerr >> obj->decay_p2nlow >> obj->decay_p2nup)) break;

        obj->population_ratio = 1;
        obj->population_ratioerr = 0;
        obj->population_ratiolow = 0;
        obj->population_ratioup = 2;
        obj->is_population_ratio_fix = 2;

        //! read info with isomer
        if (!(iss >> obj->name >> obj->z >> a >> obj->decay_hl >> obj->decay_hlerr >> obj->decay_hllow >> obj->decay_hlup >>
              obj->decay_p1n >> obj->decay_p1nerr >> obj->decay_p1nlow >> obj->decay_p1nup >>
              obj->decay_p2n >> obj->decay_p2nerr >> obj->decay_p2nlow >> obj->decay_p2nup)) break;
        obj->n=a-obj->z;

        Bool_t flagisomer = false;
        std::string tempstring(obj->name.Data());
        MemberDef* objisomer;
        if (tempstring.back()=='*') {
            objisomer=new MemberDef();
            CopyMember(obj,objisomer);
            objisomer->id = obj->id + 1;
            cout<<"ISOMER of "<<obj->name<<endl;
            if (!(iss >> objisomer->population_ratio >> objisomer->population_ratioerr >> objisomer->population_ratiolow >> objisomer->population_ratioup >> objisomer->decay_hl >> objisomer->decay_hlerr >> objisomer->decay_hllow >> objisomer->decay_hlup >>
                  objisomer->decay_p1n >> objisomer->decay_p1nerr >> objisomer->decay_p1nlow >> objisomer->decay_p1nup >>
                  objisomer->decay_p2n >> objisomer->decay_p2nerr >> objisomer->decay_p2nlow >> objisomer->decay_p2nup)) break;
            obj->name = TString(tempstring.substr(0,tempstring.length()-1).data());

            if (objisomer->population_ratio>0){
                obj->is_population_ratio_fix = 1;
                objisomer->is_population_ratio_fix = 1;
            }else{
                objisomer->population_ratio = -objisomer->population_ratio;
                obj->is_population_ratio_fix = 0;
                objisomer->is_population_ratio_fix = 0;
            }

            obj->population_ratio = 1 - objisomer->population_ratio;
            obj->population_ratioerr = objisomer->population_ratioerr;
            obj->population_ratioup = 1 - objisomer->population_ratiolow;
            obj->population_ratiolow = 1 - objisomer->population_ratioup;



            flagisomer = true;
        }

        //! ground state
        ProcessMember(obj);
        list.emplace(list.end(),obj);
        id++;

        if (flagisomer){
            //! isomeric state
            ProcessMember(objisomer);
            list.emplace(list.end(),objisomer);
            id++;

            groundstate[nisomers] = obj->id;
            isomerstate[nisomers] = objisomer->id;
            nisomers++;
        }

    }

}

void appendvectors(vector< vector<Int_t> >& pathoriginal,vector< vector<Int_t> >& pathnew,Int_t newmember,vector< vector<Int_t> >& nneupathoriginal,vector< vector<Int_t> >& nneupathnew,Int_t nneunewmember)
{

    //! Vectors containing information about the node in the decay network
    // if original vector is empty (parent nuclei)
    if (pathoriginal.size()==0){
        vector<Int_t> row;
        row.push_back(0);
        row.push_back(newmember);
        pathnew.push_back(row);
    }

    // if original vector is not empty
    for (Size_t i=0;i<pathoriginal.size();i++)
    {
        // copy original vector
        vector<Int_t> row;
        for (Size_t j=0;j<pathoriginal[i].size();j++)
        {
            row.push_back(pathoriginal[i][j]);
        }
        // add one more element to the row
        row.push_back(newmember);
        pathnew.push_back(row);
    }

    //! Vectors containing information about how many neutron emitted in a decay
    // if original vector is empty (parent nuclei)
    if (nneupathoriginal.size()==0){
        vector<Int_t> nneurow;
        nneurow.push_back(-1); // no meaning for first row
        nneurow.push_back(nneunewmember);
        nneupathnew.push_back(nneurow);
    }

    // if original vector is not empty
    for (Size_t i=0;i<nneupathoriginal.size();i++)
    {
        // copy original vector
        vector<Int_t> nneurow;
        for (Size_t j=0;j<nneupathoriginal[i].size();j++)
        {
            nneurow.push_back(nneupathoriginal[i][j]);
        }
        // add one more element to the row
        nneurow.push_back(nneunewmember);
        nneupathnew.push_back(nneurow);
    }
}

void makepath(char* inputfile)
{
    list<MemberDef*> listofdecaymember;
    list<MemberDef*>::iterator listofdecaymember_it;
    list<MemberDef*>::iterator listofdecaymember_it2;
    readinput(inputfile,listofdecaymember);

    // display information of the decay members
    for (listofdecaymember_it = listofdecaymember.begin(); listofdecaymember_it != listofdecaymember.end(); listofdecaymember_it++)
    {
        cout<<(*listofdecaymember_it)->id+1<<"\t"<<(*listofdecaymember_it)->name<<"\t"<<(*listofdecaymember_it)->z<<"\t"<<
              (*listofdecaymember_it)->n<<"\t"<<(*listofdecaymember_it)->n+(*listofdecaymember_it)->z<<"\t"<<(*listofdecaymember_it)->decay_hl<<"\t"<<(*listofdecaymember_it)->decay_lamda<<"\t"<<
              (*listofdecaymember_it)->decay_p1n<<"\t"<<(*listofdecaymember_it)->decay_p2n<<"\t"<<
              (*listofdecaymember_it)->is_decay_hl_fix<<"\t"<<(*listofdecaymember_it)->is_decay_p1n_fix<<"\t"<<(*listofdecaymember_it)->is_decay_p2n_fix<<endl;
    }

    // main loop
    for (listofdecaymember_it = listofdecaymember.begin(); listofdecaymember_it != listofdecaymember.end(); listofdecaymember_it++)
    {
        for (listofdecaymember_it2 = listofdecaymember.begin(); listofdecaymember_it2 != listofdecaymember.end(); listofdecaymember_it2++)
        {
            if (((*listofdecaymember_it)->z-(*listofdecaymember_it2)->z)==1){
                if (((*listofdecaymember_it)->n-(*listofdecaymember_it2)->n)==-1){//p0n
                    appendvectors((*listofdecaymember_it2)->path,(*listofdecaymember_it)->path,(*listofdecaymember_it)->id,
                                  (*listofdecaymember_it2)->nneupath,(*listofdecaymember_it)->nneupath,0);
                }else if (((*listofdecaymember_it)->n-(*listofdecaymember_it2)->n)==-2){//p1n
                    appendvectors((*listofdecaymember_it2)->path,(*listofdecaymember_it)->path,(*listofdecaymember_it)->id,
                                  (*listofdecaymember_it2)->nneupath,(*listofdecaymember_it)->nneupath,1);
                }else if (((*listofdecaymember_it)->n-(*listofdecaymember_it2)->n)==-3){//p2n
                    appendvectors((*listofdecaymember_it2)->path,(*listofdecaymember_it)->path,(*listofdecaymember_it)->id,
                                  (*listofdecaymember_it2)->nneupath,(*listofdecaymember_it)->nneupath,2);
                }
            }
        }
    }

    // pass the listofdecaymember to the old arrays used for global function
    npaths=0;

    Double_t minz=100000;
    Double_t minn=100000;
    Double_t maxz=0;
    Double_t maxn=0;

    Double_t expandZ=3;
    Double_t expandN=3;

    for (listofdecaymember_it = listofdecaymember.begin(); listofdecaymember_it != listofdecaymember.end(); listofdecaymember_it++)
    {
        if ((*listofdecaymember_it)->z<minz) minz=(*listofdecaymember_it)->z;
        if ((*listofdecaymember_it)->n<minn) minn=(*listofdecaymember_it)->n;

        if ((*listofdecaymember_it)->z>maxz) maxz=(*listofdecaymember_it)->z;
        if ((*listofdecaymember_it)->n>maxn) maxn=(*listofdecaymember_it)->n;

        for (Size_t i=0;i<(*listofdecaymember_it)->path.size();i++){
            ndecay[npaths]=(*listofdecaymember_it)->path[i].size();
            for (Size_t j=0;j<(*listofdecaymember_it)->path[i].size();j++){
                decaymap[npaths][(Int_t)j]=(*listofdecaymember_it)->path[i][j];
            }

            for (Size_t j=1;j<(*listofdecaymember_it)->path[i].size();j++){
                nneu[npaths][(Int_t)j-1]=(*listofdecaymember_it)->nneupath[i][j];
            }

            npaths++;
        }
    }

    TCanvas* cc=new TCanvas("cc","",900,700) ;
    Double_t xrange[2]={minn-expandN,maxn+expandN};
    Double_t yrange[2]={minz-expandZ,maxz+expandZ};
    plainPlot(cc,xrange,yrange);

    TLatex latex;
    latex.SetTextAlign(12);
    latex.SetTextSize(0.025);

    TArrow arr;
    arr.SetLineColor(6);
    arr.SetFillColor(2);
    //arr.SetLineStyle(2);
    // Plot the flow
    npaths=0;

    Int_t isplotisoarr[kmaxparms];
    for (listofdecaymember_it = listofdecaymember.begin(); listofdecaymember_it != listofdecaymember.end(); listofdecaymember_it++)
    {
        cout<<"********* Go for Isotope "<<(*listofdecaymember_it)->id<<" : "<<(*listofdecaymember_it)->n+(*listofdecaymember_it)->z<<(*listofdecaymember_it)->name<<endl;
        Int_t isplotiso=0;
        for (Size_t i=0;i<(*listofdecaymember_it)->path.size();i++){
            Double_t prevz=0;
            Double_t prevn=0;
            Double_t previd=0;
            Double_t prevp1n=0;
            Double_t prevp2n=0;

            cout<<"row "<<i<<" = ";


            Bool_t isplot=true;

            for (Size_t j=0;j<(*listofdecaymember_it)->path[i].size();j++){
                cout<<(*listofdecaymember_it)->path[i][j]<<"-";
                //! draw arrow
                for (listofdecaymember_it2 = listofdecaymember.begin(); listofdecaymember_it2 != listofdecaymember.end(); listofdecaymember_it2++)
                {
                    if ((*listofdecaymember_it)->path[i][j]==(*listofdecaymember_it2)->id){

                        Double_t presz=(*listofdecaymember_it2)->z;
                        Double_t presn=(*listofdecaymember_it2)->n;
                        Double_t presid=(*listofdecaymember_it2)->id;


                        if (prevz+prevn==presz+presn){
                            if ((1-prevp1n+prevp2n)==0) {
                                isplot=false;
                            }
                        }else if(presz+presn==prevz+prevn-1){
                            if (prevp1n==0) {
                                isplot=false;
                            }
                        }else if((presz+presn==prevz+prevn-2)){
                            if (prevp2n==0) {
                                isplot=false;
                            }
                        }
                         cout<<presz+presn<<(*listofdecaymember_it2)->name<<"\t";

                        if (prevz!=0&&isplot){
                            arr.DrawArrow(prevn,prevz,presn,presz,0.01,">");
                            //! A trick for plotting, draw at the end of each track another arrow
                            arr.DrawArrow(presn,presz,presn-1,presz+1,0.01,">");
                            isplotiso++;
                            //latex.DrawLatex((*listofdecaymember_it2)->n-0.5,(*listofdecaymember_it2)->z,Form("%s",(*listofdecaymember_it2)->name.Data()));
                        }

                        prevz=presz;
                        prevn=presn;
                        prevp1n=(*listofdecaymember_it2)->decay_p1n;
                        prevp2n=(*listofdecaymember_it2)->decay_p2n;
                        previd=presid;
                    }
                }

            }
            cout<<endl;
            npaths++;
        }
        isplotisoarr[(*listofdecaymember_it)->id]=isplotiso;

    }
    for (listofdecaymember_it = listofdecaymember.begin(); listofdecaymember_it != listofdecaymember.end(); listofdecaymember_it++)
    {
        if ((*listofdecaymember_it)->id==0) latex.DrawLatex((*listofdecaymember_it)->n-0.5,(*listofdecaymember_it)->z,Form("%s",(*listofdecaymember_it)->name.Data()));
        else if (isplotisoarr[(*listofdecaymember_it)->id]>0) latex.DrawLatex((*listofdecaymember_it)->n-0.5,(*listofdecaymember_it)->z,Form("%s",(*listofdecaymember_it)->name.Data()));
    }

    cout<<endl;


    // number of ri
    knri=listofdecaymember.size();

    Int_t k=0;
    // put in parms variable
    for (listofdecaymember_it = listofdecaymember.begin(); listofdecaymember_it != listofdecaymember.end(); listofdecaymember_it++)
    {
        char tempa[100];
        sprintf(tempa,"%i",(*listofdecaymember_it)->z+(*listofdecaymember_it)->n);
        riname[k] = TString(tempa)+(*listofdecaymember_it)->name;
        parms[k]=(*listofdecaymember_it)->decay_lamda;
        parms[knri+k]=(*listofdecaymember_it)->decay_p1n;
        parms[knri*2+k]=(*listofdecaymember_it)->decay_p2n;
        parms[knri*3+k]=(*listofdecaymember_it)->population_ratio;


        parmserr[k]=(*listofdecaymember_it)->decay_lamdaerr;
        parmserr[knri+k]=(*listofdecaymember_it)->decay_p1nerr;
        parmserr[knri*2+k]=(*listofdecaymember_it)->decay_p2nerr;
        parmserr[knri*3+k]=(*listofdecaymember_it)->population_ratioerr;

        parmsmin[k]=(*listofdecaymember_it)->decay_lamdalow;
        parmsmin[knri+k]=(*listofdecaymember_it)->decay_p1nlow;
        parmsmin[knri*2+k]=(*listofdecaymember_it)->decay_p2nlow;
        parmsmin[knri*3+k]=(*listofdecaymember_it)->population_ratiolow;

        parmsmax[k]=(*listofdecaymember_it)->decay_lamdaup;
        parmsmax[knri+k]=(*listofdecaymember_it)->decay_p1nup;
        parmsmax[knri*2+k]=(*listofdecaymember_it)->decay_p2nup;
        parmsmax[knri*3+k]=(*listofdecaymember_it)->population_ratioup;

        isparmsfix[k]=(*listofdecaymember_it)->is_decay_lamda_fix;
        isparmsfix[knri+k]=(*listofdecaymember_it)->is_decay_p1n_fix;
        isparmsfix[knri*2+k]=(*listofdecaymember_it)->is_decay_p2n_fix;
        isparmsfix[knri*3+k]=(*listofdecaymember_it)->is_population_ratio_fix;
        k++;
    }

    k=0;
    for (Int_t i=0;i<knri*4;i++){
        cout<<"parms "<<i<<" : ";
        cout<<riname[i]<<"\t"<<parms[i]<<"\t"<<parmserr[i]<<"\t"<<parmsmin[i]<<"\t"<<parmsmax[i]<<"\t"<<isparmsfix[i]<<endl;
        if (isparmsfix[i]!=2){
            indexparmsactive[i]=k;
            k++;
        }
    }
    nparmsactive = k;

    //! delete all info
    for (listofdecaymember_it = listofdecaymember.begin(); listofdecaymember_it != listofdecaymember.end(); listofdecaymember_it++){
        delete (*listofdecaymember_it);
    }
    listofdecaymember.clear();




}
