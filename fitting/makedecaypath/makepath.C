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

#include "TCutG.h"
#include "TLatex.h"
#include "TMath.h"
#include "TCanvas.h"

#include "TGraph.h"
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/GSLSimAnMinimizer.h"
#include "Math/Functor.h"

#include <list>

#include <sstream>
#include <string>

using namespace std;



//Global variable for core function
const Int_t kmaxndecay=100;
const Int_t kmaxpaths=100;
const Int_t kmaxparms=500;

Int_t npaths=100;
Int_t ndecay[kmaxpaths];
Int_t decaymap[kmaxpaths][kmaxndecay];
Int_t nneu[kmaxpaths][kmaxndecay];
Int_t knri=200;

Double_t parms[kmaxparms];
Double_t parmserr[kmaxparms];
Double_t parmsmax[kmaxparms];
Double_t parmsmin[kmaxparms];
Int_t isparmsfix[kmaxparms];
Double_t mcparms[kmaxparms];

Bool_t reject=false;
Double_t rejectrange=0.05;//first 50 ms

typedef struct {
    // idendification
    Int_t id;
    Int_t z;
    Int_t n;
    TString name;

    // decay properies
    Double_t decay_hl;
    Double_t decay_lamda;
    Double_t decay_p1n;
    Double_t decay_p2n;

    Double_t decay_hlerr;
    Double_t decay_lamdaerr;
    Double_t decay_p1nerr;
    Double_t decay_p2nerr;

    Double_t decay_hlup;
    Double_t decay_lamdaup;
    Double_t decay_p1nup;
    Double_t decay_p2nup;

    Double_t decay_hllow;
    Double_t decay_lamdalow;
    Double_t decay_p1nlow;
    Double_t decay_p2nlow;

    // fit options 0-vary, 1-fix with error propagation, 2-fix without error propagation
    Int_t is_decay_hl_fix;
    Int_t is_decay_lamda_fix;
    Int_t is_decay_p1n_fix;
    Int_t is_decay_p2n_fix;

    // paths to this ri
    vector< vector<Int_t> > path;
    vector< vector<Int_t> > nneupath;
} MemberDef;


void readinput(char* infilename, list<MemberDef*>& list)
{
    Int_t id=0;
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
        if (!(iss >> obj->name >> obj->z >> a >> obj->decay_hl >> obj->decay_hlerr >> obj->decay_hllow >> obj->decay_hlup >>
              obj->decay_p1n >> obj->decay_p1nerr >> obj->decay_p1nlow >> obj->decay_p1nup >>
              obj->decay_p2n >> obj->decay_p2nerr >> obj->decay_p2nlow >> obj->decay_p2nup)) break;

        // further processing
        obj->n=a-obj->z;

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
        }else if (obj->decay_p1n==0){//exclude decay half-life =0;
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
        list.emplace(list.end(),obj);

        id++;
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
    for (Int_t i=0;i<pathoriginal.size();i++)
    {
        // copy original vector
        vector<Int_t> row;
        for (Int_t j=0;j<pathoriginal[i].size();j++)
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
    for (Int_t i=0;i<nneupathoriginal.size();i++)
    {
        // copy original vector
        vector<Int_t> nneurow;
        for (Int_t j=0;j<nneupathoriginal[i].size();j++)
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
    for (listofdecaymember_it = listofdecaymember.begin(); listofdecaymember_it != listofdecaymember.end(); listofdecaymember_it++)
    {
        for (Int_t i=0;i<(*listofdecaymember_it)->path.size();i++){
            ndecay[npaths]=(*listofdecaymember_it)->path[i].size();
            for (Int_t j=0;j<(*listofdecaymember_it)->path[i].size();j++){
                decaymap[npaths][j]=(*listofdecaymember_it)->path[i][j];
            }

            for (Int_t j=1;j<(*listofdecaymember_it)->path[i].size();j++){
                nneu[npaths][j-1]=(*listofdecaymember_it)->nneupath[i][j];
            }

            npaths++;
        }
    }
    // number of ri
    knri=listofdecaymember.size();

    Int_t k=0;
    // put in parms variable
    for (listofdecaymember_it = listofdecaymember.begin(); listofdecaymember_it != listofdecaymember.end(); listofdecaymember_it++)
    {
        parms[k]=(*listofdecaymember_it)->decay_lamda;
        parms[knri+k]=(*listofdecaymember_it)->decay_p1n;
        parms[knri*2+k]=(*listofdecaymember_it)->decay_p2n;


        parmserr[k]=(*listofdecaymember_it)->decay_lamdaerr;
        parmserr[knri+k]=(*listofdecaymember_it)->decay_p1nerr;
        parmserr[knri*2+k]=(*listofdecaymember_it)->decay_p2nerr;

        parmsmin[k]=(*listofdecaymember_it)->decay_lamdalow;
        parmsmin[knri+k]=(*listofdecaymember_it)->decay_p1nlow;
        parmsmin[knri*2+k]=(*listofdecaymember_it)->decay_p2nlow;

        parmsmax[k]=(*listofdecaymember_it)->decay_lamdaup;
        parmsmax[knri+k]=(*listofdecaymember_it)->decay_p1nup;
        parmsmax[knri*2+k]=(*listofdecaymember_it)->decay_p2nup;

        isparmsfix[k]=(*listofdecaymember_it)->is_decay_lamda_fix;
        isparmsfix[knri+k]=(*listofdecaymember_it)->is_decay_p1n_fix;
        isparmsfix[knri*2+k]=(*listofdecaymember_it)->is_decay_p2n_fix;
        k++;
    }

    for (Int_t i=0;i<knri*3;i++){
        cout<<"parms "<<i<<" : ";
        cout<<parms[i]<<"\t"<<parmserr[i]<<"\t"<<parmsmin[i]<<"\t"<<parmsmax[i]<<"\t"<<isparmsfix[i]<<endl;
    }

    /*
    for (Int_t i=0;i<npaths;i++){
        cout<<"path "<<i<<" : ";
        for (Int_t j=0;j<ndecay[i];j++){
            cout<<decaymap[i][j]<<"\t";
        }
        cout<<endl;
    }
    */




/*
    // display information of the decay members
    cout<<"\n\n"<<endl;
    for (listofdecaymember_it = listofdecaymember.begin(); listofdecaymember_it != listofdecaymember.end(); listofdecaymember_it++)
    {
        cout<<"ID = "<<(*listofdecaymember_it)->id<< "\tRI = "<<(*listofdecaymember_it)->name<<(*listofdecaymember_it)->z+(*listofdecaymember_it)->n<<" : "<<endl;
        cout<<"npath="<<(*listofdecaymember_it)->path.size()<<endl;
        cout<<"nneupath="<<(*listofdecaymember_it)->nneupath.size()<<endl;
        for (Int_t i=0;i<(*listofdecaymember_it)->path.size();i++){
            cout<<"row "<<i<<" = ";
            for (Int_t j=0;j<(*listofdecaymember_it)->path[i].size();j++){
                cout<<(*listofdecaymember_it)->path[i][j]<<"\t";
            }
            cout<<endl;

            cout<<"row "<<i<<" = ";
            for (Int_t j=0;j<(*listofdecaymember_it)->path[i].size();j++){
                cout<<(*listofdecaymember_it)->nneupath[i][j]<<"\t";
            }
            cout<<endl;
        }

        cout<<"***************************\n"<<endl;
    }
*/
}

//! Global Bateaman function
Double_t corefcn(Int_t ndecay,Int_t* decaymap,Int_t* nneu, Double_t* b1n,Double_t* b2n,Double_t* lamda,Double_t N0,Double_t t){
    Double_t fcnret=0;

    Double_t factor1=1.;
    //only parrent decay p2n
    for (int i=0;i<ndecay-1;i++){
        if (nneu[i]==0){
            factor1=factor1 * (1-b1n[decaymap[i]]-b2n[decaymap[i]])*lamda[decaymap[i]];
        }else if (nneu[i]==1){
            factor1=factor1 * b1n[decaymap[i]]*lamda[decaymap[i]];
        }else{
            factor1=factor1 * b2n[decaymap[i]]*lamda[decaymap[i]];
        }
    }

    Double_t factor2=0;
    for (int i=0;i<ndecay;i++){
        Double_t factor2i=exp(-lamda[decaymap[i]]*t);
        Double_t factor2ij=1;
        for (int j=0;j<ndecay;j++)
            if (j!=i) factor2ij=factor2ij*(lamda[decaymap[j]]-lamda[decaymap[i]]);
        factor2=factor2+factor2i/factor2ij;
    }

    fcnret=factor1*N0*factor2;
    return fcnret;
}

//! Global function
Double_t fcn_decay(Double_t *x, Double_t *par) {
    Double_t bkg=par[knri*3+1];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;

    Double_t* pn=&par[knri];
    Double_t* lamda=par;
    Double_t* p2n=&par[knri*2];
    Double_t N0=par[knri*3]/par[0];

    //! Parent nuclei
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0]);

    for (Int_t i=0;i<npaths;i++){
        returnval+=lamda[decaymap[i][ndecay[i]-1]]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejectrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

void fitter()
{

    //! construct params

    Double_t lowerlimit=-10;
    Double_t upperlimit=30;

    makepath("testinput.txt");
    //! define function
    TF1* fB=new TF1("fB",fcn_decay,lowerlimit,upperlimit,knri*3+2);

    fB->SetNpx(2000);
    fB->SetLineWidth(2);
    fB->SetLineColor(8);

    //! initializing parameters
    //fB->SetParameter(i,parms[i]);
    //fB->SetParLimits(i,parmsmin[i],parmsmax[i]);

}

