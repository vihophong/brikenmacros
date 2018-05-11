#include <math.h>
#include <stdlib.h>
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AIDA.h"
#include "BELEN.h"
#include "Clover.h"
#include <map>
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class AIDAHit+;
#pragma link C++ class AIDACluster+;
#pragma link C++ class AIDA+;
#pragma link C++ class BELENHit+;
#pragma link C++ class BELEN+;
#pragma link C++ class CloverHit+;
#pragma link C++ class Clover+;
#endif



typedef struct {
    double T; 	 // Calibrated time
    double Tcorr; //correlated time
    double x,y,z;// number of pixel for AIDA, or number of tube for BELEN
    int type;
    int type2;
    int evt;
} datatype;

void convertsimdata(char* infile, char* outfile,Int_t opt=0)
{
    gROOT->ProcessLine(".L AIDA.h+");
    gROOT->ProcessLine(".L BELEN.h+");
    gROOT->ProcessLine(".L Clover.h+");


    std::multimap < double,std::pair< unsigned int, datatype> > fionbetaMap;
    std::multimap < double,std::pair< unsigned int, datatype> >::iterator fionbetaMap_it;

    std::multimap < double,std::pair< unsigned int, datatype> > fneuMap;
    std::multimap < double,std::pair< unsigned int, datatype> >::iterator fneuMap_it;


    datatype readbackion;
    datatype readbackbeta;
    datatype readbackneu;

    TTree* readbacktreebeta=0;
    TTree* readbacktreeion=0;
    TTree* readbacktreeneu=0;

    TFile* f = new TFile(infile);
    f->GetObject("ion",readbacktreeion);
    f->GetObject("beta",readbacktreebeta);
    f->GetObject("neutron",readbacktreeneu);

    readbacktreeion->SetBranchAddress("ion",&readbackion);
    readbacktreebeta->SetBranchAddress("beta",&readbackbeta);
    readbacktreeneu->SetBranchAddress("neutron",&readbackneu);

    Long64_t nentriesion=readbacktreeion->GetEntries();
    Long64_t nentriesbeta=readbacktreebeta->GetEntries();
    Long64_t nentriesneu=readbacktreeneu->GetEntries();


    if (opt==0){
        for (unsigned int jentry=0;jentry<nentriesion;jentry++){
            readbacktreeion->GetEntry(jentry);
            datatype aidadata;
            aidadata.x=readbackion.x;
            aidadata.y=readbackion.y;
            aidadata.z=readbackion.z;
            aidadata.T=readbackion.T;
            aidadata.type=readbackion.type;
            aidadata.type2=readbackion.type2;
            aidadata.evt=readbackion.evt;
            fionbetaMap.insert(make_pair(aidadata.T,make_pair(4,aidadata)));
        }

        for (unsigned int jentry=0;jentry<nentriesbeta;jentry++){
            readbacktreebeta->GetEntry(jentry);
            datatype aidadata;
            aidadata.x=readbackbeta.x;
            aidadata.y=readbackbeta.y;
            aidadata.z=readbackbeta.z;
            aidadata.T=readbackbeta.T;
            aidadata.type=readbackbeta.type;
            aidadata.type2=readbackbeta.type2;
            aidadata.evt=readbackbeta.evt;
            fionbetaMap.insert(make_pair(aidadata.T,make_pair(5,aidadata)));
        }
        cout<<"\ntotal number of implant-beta = "<<fionbetaMap.size()<<endl;
        cout<<"total number of neutron = "<<fneuMap.size()<<endl;

        //! making files
        TFile* fout=new TFile(outfile,"RECREATE");
        fout->cd();
        AIDASimpleStruct* aida=new AIDASimpleStruct;
        TTree* treeaida=new TTree("aida","aida simple tree");
        treeaida->Branch("aida",&aida);
        for (fionbetaMap_it=fionbetaMap.begin();fionbetaMap_it!=fionbetaMap.end();fionbetaMap_it++){
            datatype aidadata=fionbetaMap_it->second.second;
            Int_t id=(Int_t)fionbetaMap_it->second.first;
            aida->SetID(id);
            unsigned long long timedata=(unsigned long long)(aidadata.T*1e9);
            aida->SetTimestamp(timedata);
            aida->SetHitPosition(aidadata.x,aidadata.y,aidadata.z);
            aida->SetDtIon(999999999999);
            treeaida->Fill();
        }
        treeaida->Write();
        fout->Close();

    }else if (opt==1){
        TFile* fout=new TFile(outfile,"RECREATE");
        fout->cd();
        //! Book tree and histograms
        TTree* treeneutron=new TTree("neutron","tree neutron");
        TTree* treegamma=new TTree("gamma","tree gamma");
        TTree* treeanc=new TTree("anc","tree anc");

        BELENHit* neutron=new BELENHit;
        CloverHit* gamma=new CloverHit;
        BELENHit* anc=new BELENHit;
        treeneutron->Branch("neutron",&neutron);
        treegamma->Branch("gamma",&gamma);
        treeanc->Branch("anc",&anc);
        treeneutron->BranchRef();
        treegamma->BranchRef();
        treeanc->BranchRef();
        for (unsigned int jentry=0;jentry<nentriesneu;jentry++){
            readbacktreeneu->GetEntry(jentry);
            neutron->SetEnergy(600);
            neutron->SetTimestamp((unsigned long long)(readbackneu.T*1e9));
            treeneutron->Fill();
        }
        treeneutron->Write();
        treegamma->Write();
        treeanc->Write();
        fout->Close();
    }

}

void random3(){
    Double_t hl=0.128;
    //Double_t hl=0.5;
    //generate histo with background
    TH1F* hdecayc=new TH1F("hdecayc","hdecayc",2000,-10-1,10+1);

    TFile* fout=new TFile("outhist.root","RECREATE");
    Double_t tp1n=0;
    Double_t tp2n=0;
    Double_t tall=0;

    TTree* treeb=new TTree("treeb","treeb");
    treeb->Branch("x",&tall,"x/D");
    TTree* treebb=new TTree("treebb","treebb");
    treebb->Branch("x",&tall,"x/D");

    TTree* treep1n=new TTree("treep1n","treep1n");
    treep1n->Branch("x",&tp1n,"x/D");
    TTree* treep2n=new TTree("treep2n","treep2n");
    treep2n->Branch("x",&tp2n,"x/D");
    TRandom3* r3=new TRandom3();

    for (unsigned int jentry=0;jentry<150623;jentry++){
        tall=r3->Exp(hl/log(2));
        hdecayc->Fill(tall);
        Double_t pratio=r3->Rndm()*100;
        if (pratio<90) {
            tall=r3->Rndm()*22-11;
            treeb->Fill();
            if (tall<0) treebb->Fill();
        }
        treeb->Fill();
    }
    hdecayc->Write();
    treeb->Write();
    treebb->Write();
    treep1n->Write();
    treep2n->Write();
    fout->Close();
}


