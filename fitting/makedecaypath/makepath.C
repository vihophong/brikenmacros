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

using namespace std;

const Int_t kmaxndecay=50;
const Int_t kmaxpaths=100;

typedef struct {
    // idendification
    Int_t z;
    Int_t n;
    TString name;

    // decay properies
    Double_t decay_hl;
    Double_t decay_p1n;
    Double_t decay_p2n;

    // fit options
    Bool_t is_decay_hl_vary;
    Bool_t is_decay_p1n_vary;
    Bool_t is_decay_p2n_vary;

    // paths to this ri
    Int_t npaths;
    Int_t ndecay[kmaxpaths];
    Int_t decaymap[kmaxpaths][kmaxndecay];
    Int_t nneu[kmaxpaths][kmaxndecay];

} MemberDef;
void makepath()
{
    MemberDef* parent=new MemberDef();

    list<MemberDef*> listofdescendant;
    list<MemberDef*>::iterator listofdescendant_it;
    //! Update list of Parent and Descendants with identification and decay properties
    parent->z = 48;
    parent->n = 82;

    MemberDef* obj=new MemberDef();
    obj->npaths=3;
    listofdescendant.emplace(listofdescendant.begin(),obj);
    MemberDef* obj2=new MemberDef();
    obj2->npaths=5;
    listofdescendant.emplace(listofdescendant.begin(),obj2);
    cout<<listofdescendant.size()<<endl;
    for (listofdescendant_it = listofdescendant.begin(); listofdescendant_it != listofdescendant.end(); listofdescendant_it++)
        cout << (*listofdescendant_it)->npaths <<endl;
}
