#include "stdio.h"
#include "string.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TTree.h"
void offlinemlhfit()
{
    TString ri[100];
    Int_t binning[100];
    Int_t nri=17;

    ri[0]="Ag130";binning[0]=225;
    ri[1]="Ag131";binning[1]=150;

    ri[2]="Cd129";binning[2]=450;
    ri[3]="Cd130";binning[3]=900;
    ri[4]="Cd131";binning[4]=450;
    ri[5]="Cd132";binning[5]=225;
    ri[6]="Cd133";binning[6]=225;

    ri[7]="In131";binning[7]=450;
    ri[8]="In132";binning[8]=450;
    ri[9]="In133";binning[9]=225;
    ri[10]="In134";binning[10]=225;
    ri[11]="In135";binning[11]=225;

    ri[12]="Sn134";binning[12]=150;
    ri[13]="Sn135";binning[13]=150;
    ri[14]="Sn136";binning[14]=225;
    ri[15]="Sn137";binning[15]=225;
    ri[16]="Sn138";binning[16]=100;

    gROOT->ProcessLine(".L mlhfitv4.C");
    for (Int_t i=0;i<nri;i++){
        cout<<Form("mlhfitv4(\"%s\",\"../../fittrees/outree%s.root\",\"../../fitparms/parms%s.txt\",\"../../fitresultmlh/fitresult%s.root\")",(char*)ri[i].Data(),(char*)ri[i].Data(),(char*)ri[i].Data(),(char*)ri[i].Data())<<endl;
        gROOT->ProcessLine(Form("mlhfitv4(\"%s\",\"../../fittrees/outree%s.root\",\"../../fitparms/parms%s.txt\",\"../../fitresultmlh/fitresult%s.root\")",(char*)ri[i].Data(),(char*)ri[i].Data(),(char*)ri[i].Data(),(char*)ri[i].Data()));
        gSystem->Sleep(1000);
    }
}
