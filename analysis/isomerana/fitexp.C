#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TCutG.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

#include "TF1.h"

void fitexp(TH1F* hin,Double_t rlow,Double_t rup)
{
    TString fdef("pol0(0)+expo(1)");
    hin->Fit("pol0","LQRE","",rlow,rup);
    Double_t p0=hin->GetFunction("pol0")->GetParameter(0);
    hin->Fit("expo","LQRE","",rlow,rup);
    Double_t p1=hin->GetFunction("expo")->GetParameter(0);
    Double_t p2=hin->GetFunction("expo")->GetParameter(1);
    Double_t p2err=hin->GetFunction("expo")->GetParError(1);
    cout<<-log(2)/p2<<"\t"<<log(2)/p2/p2*p2err<<endl;

    TF1* ffinal=new TF1("ffinal",(char*)fdef.Data(),rlow,rup);
    ffinal->SetParameter(0,p0);
    ffinal->SetParameter(1,p1);
    ffinal->SetParameter(2,p2);
    hin->Fit("ffinal","LQRE");
    p0=hin->GetFunction("ffinal")->GetParameter(0);
    p2=hin->GetFunction("ffinal")->GetParameter(2);
    p2err=hin->GetFunction("ffinal")->GetParError(2);
    cout<<p0<<endl;
    cout<<-log(2)/p2<<"\t"<<log(2)/p2/p2*p2err<<endl;
}
