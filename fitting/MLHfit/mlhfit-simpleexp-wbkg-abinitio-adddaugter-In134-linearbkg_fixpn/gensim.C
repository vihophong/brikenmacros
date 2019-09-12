#include <TROOT.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TRandom3.h>
void gensim(){
    Double_t nentries=6000;

    Double_t rbkg=66.67;//%
    Double_t p1nbkg=10;
    Double_t p2nbkg=1;


    Double_t hl=0.128;
    //Double_t hl=0.5;
    Double_t p1n=50;
    Double_t p2n=25;
    //generate histo with background
    TH1F* hdecayc=new TH1F("hdecayc","hdecayc",2000,-10,10);

    TFile* fout=new TFile("outhist.root","RECREATE");
    Double_t tall=0;
    Int_t nmult=0;

    TTree* treeb=new TTree("tree","tree");
    treeb->Branch("x",&tall,"x/D");
    treeb->Branch("y",&nmult,"y/I");

    TRandom3* r3=new TRandom3(0);

    for (unsigned int jentry=0;jentry<nentries;jentry++){
        Double_t bkgratio=r3->Rndm()*100;
        if (bkgratio<rbkg){
            tall=r3->Rndm()*20-10;
            Double_t bkgmultratio=r3->Rndm()*100;
            if (bkgmultratio<p1nbkg) {
                nmult=1;
            }else if (bkgmultratio>=p1nbkg&&bkgmultratio<p1nbkg+p2nbkg) {
                nmult=2;
            }else{
                nmult=0;
            }
        }else{
            tall=r3->Exp(hl/log(2));
            hdecayc->Fill(tall);
            Double_t pratio=r3->Rndm()*100;

            if (pratio<p1n) {
                nmult=1;
            }else if (pratio>=p1n&&pratio<p1n+p2n) {
                nmult=2;
            }else{
                nmult=0;
            }
        }

        treeb->Fill();
    }
    hdecayc->Write();
    treeb->Write();
    fout->Close();
}

