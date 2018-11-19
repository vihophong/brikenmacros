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
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/GSLSimAnMinimizer.h"
#include "Math/Functor.h"

#include <list>

#include <sstream>
#include <string>


#include "plainPlot.C"

using namespace std;



//Global variable for core function
const Int_t kmaxndecay=500;
const Int_t kmaxpaths=500;
const Int_t kmaxparms=1000;

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
    // Plot the flow
    npaths=0;
    for (listofdecaymember_it = listofdecaymember.begin(); listofdecaymember_it != listofdecaymember.end(); listofdecaymember_it++)
    {
        cout<<"********* Go for Isotope "<<(*listofdecaymember_it)->id<<" : "<<(*listofdecaymember_it)->n+(*listofdecaymember_it)->z<<(*listofdecaymember_it)->name<<endl;
        if ((*listofdecaymember_it)->id==0) latex.DrawLatex((*listofdecaymember_it)->n-0.5,(*listofdecaymember_it)->z,Form("^{%d}%s",(*listofdecaymember_it)->z+(*listofdecaymember_it)->n,(*listofdecaymember_it)->name.Data()));
        for (Int_t i=0;i<(*listofdecaymember_it)->path.size();i++){
            Double_t prevz=0;
            Double_t prevn=0;
            Double_t previd=0;
            Double_t prevp1n=0;
            Double_t prevp2n=0;

            cout<<"row "<<i<<" = ";


            Bool_t isplot=true;

            for (Int_t j=0;j<(*listofdecaymember_it)->path[i].size();j++){
                cout<<(*listofdecaymember_it)->path[i][j]<<"-";
                //! draw arrow
                for (listofdecaymember_it2 = listofdecaymember.begin(); listofdecaymember_it2 != listofdecaymember.end(); listofdecaymember_it2++)
                {
                    if ((*listofdecaymember_it)->path[i][j]==(*listofdecaymember_it2)->id){

                        Double_t presz=(*listofdecaymember_it2)->z;
                        Double_t presn=(*listofdecaymember_it2)->n;
                        Double_t presid=(*listofdecaymember_it2)->id;


                        if (prevz+prevn==presz+presn){
                            if ((1-prevp1n+prevp2n)==0) isplot=false;
                        }else if(presz+presn==prevz+prevn-1){
                            if (prevp1n==0) isplot=false;
                        }else if((presz+presn==prevz+prevn-2)){
                            if (prevp2n==0) isplot=false;
                        }
                         cout<<presz+presn<<(*listofdecaymember_it2)->name<<"\t";

                        if (prevz!=0&&isplot){
                            latex.DrawLatex((*listofdecaymember_it2)->n-0.5,(*listofdecaymember_it2)->z,Form("^{%d}%s",(*listofdecaymember_it2)->z+(*listofdecaymember_it2)->n,(*listofdecaymember_it2)->name.Data()));
                            arr.DrawArrow(prevn,prevz,presn,presz,0.01,">");
                            //! A trick for plotting, draw at the end of each track another arrow
                            arr.DrawArrow(presn,presz,presn-1,presz+1,0.01,">");
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
    }
    cout<<endl;


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
