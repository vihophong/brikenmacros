#include <iostream>
#include <iomanip>
#include <string>
#include <signal.h>
#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TCutG.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TStopwatch.h"
#include "TClonesArray.h"
#include "TVectorD.h"
#include "TSystem.h"
#include <fstream>

#include "DataStructNew.h"
#include "DataStructNewLinkDef.h"
#define MaxNRI 1000

void getdeadtime(char* listfile,Double_t* deadtimein,TH2F* h2){
    for (int i=0;i<6;i++) deadtimein[i]=0.;
    std::ifstream ifs(listfile);
    string filelist[1000];
    Int_t nfiles=0;
    while (!ifs.eof()){
        ifs>>filelist[nfiles];
        //cout<<filelist[nfiles]<<endl;
        nfiles++;
    }
    nfiles=nfiles-1;
    //cout<<"There are "<<nfiles<<" files in total!"<<endl;
    for (Int_t i=0;i<nfiles;i++){
        TFile *_file0 = TFile::Open(filelist[i].c_str());
        TVectorD* deadtimecontainer;
        deadtimecontainer=(TVectorD*) gDirectory->Get(Form("deadtime"));
        Double_t* deadtimearray=deadtimecontainer->GetMatrixArray();
        for (Int_t i=0;i<6;i++)
            deadtimein[i]+=deadtimearray[i];
        h2->Fill(i,(Double_t)deadtimearray[0]/(Double_t)deadtimearray[1]*100);
        _file0->Close();
    }
}


void chaineach(TChain* ch, char* listfile,char* pid){
  std::ifstream ifs(listfile);
  string filelist[1000];
  Int_t nfiles=0;
  while (!ifs.eof()){
      ifs>>filelist[nfiles];
      //cout<<filelist[nfiles]<<endl;
      nfiles++;
  }
  nfiles=nfiles-1;
  //cout<<"There are "<<nfiles<<" files in total!"<<endl;
  for (Int_t i=0;i<nfiles;i++){
      char tempchar2[1000];
      sprintf(tempchar2,"%s/tree%s",filelist[i].c_str(),pid);
      ch->Add(tempchar2);
  }
}

void makedecayx2hist(char* listfile, char* pidfile,char* outfile){
    gROOT->ProcessLine(".L BELEN.h+");
    gROOT->ProcessLine(".L Clover.h+");
    gROOT->ProcessLine(".L AIDA.h+");
    gROOT->ProcessLine(".L AIDA.h+");
    gROOT->ProcessLine(".L AIDA.h+");
    gROOT->ProcessLine(".L AIDA.h+");
    gROOT->ProcessLine(".L DataStructNew.h+");

    //! read pid file
    std::ifstream ifspid(pidfile);
    Int_t nri=0;
    Int_t nbinszet,nbinsaoq;
    Double_t zetrange[2];
    Double_t aoqrange[2];
    ifspid>>nbinsaoq>>aoqrange[0]>>aoqrange[1]>>nbinszet>>zetrange[0]>>zetrange[1];
    ifspid>>nri;

    Int_t enablepid[nri];
    Int_t enablepid2[nri];

    TString nameri[nri];
    TString latexnametri[nri];
    Double_t parmsri[nri][7];
    TString tempriname,tempria;

    Double_t halflife[nri];
    for (Int_t i=0;i<nri;i++){
        ifspid>>enablepid[i]>>enablepid2[i]>>tempria>>tempriname>>halflife[i];
        for(Int_t j=0;j<7;j++) ifspid>>parmsri[i][j];
        nameri[i]=tempriname+tempria;
        latexnametri[i]=TString("^{")+tempria+TString("}"+tempriname);
        for(Int_t j=0;j<7;j++) cout<<"\t"<<parmsri[i][j];
        //chaineach(listfile,(char*)nameri[i].Data());
    }

    //! data container
    Int_t idx=0;
    IonBetaMult* imp=new IonBetaMult;
    IonBeta* decay=new IonBeta;
    //! tree container
    TChain* chainall=0;
    TChain* chainRI[nri];
    for (Int_t i=0;i<nri;i++) chainRI[i]=0;
    //! get chain
    chainall = new TChain("tree");
    chaineach(chainall,listfile,(char*)TString("").Data());
    for (Int_t i=0;i<nri;i++) {
        char tempchar1[1000];
        sprintf(tempchar1,"tree%s",(char*)nameri[i].Data());
        chainRI[i]=new TChain(tempchar1);
        chaineach(chainRI[i],listfile,(char*)nameri[i].Data());
    }

    //! set branch
    chainall->SetBranchAddress("idx",&idx);
    chainall->SetBranchAddress("ion",&imp);
    chainall->SetBranchAddress("beta",&decay);

    for (Int_t i=0;i<nri;i++){
        chainRI[i]->SetBranchAddress("idx",&idx);
        chainRI[i]->SetBranchAddress("ion",&imp);
        chainRI[i]->SetBranchAddress("beta",&decay);
    }


    Long64_t nentriesall=chainall->GetEntries();
    Long64_t nentriesRI[nri];
    cout<<"\nFinished initializing! \n\nPrinting number of RIs entries:"<<endl;
    for (Int_t i=0;i<nri;i++){
        nentriesRI[i]=chainRI[i]->GetEntries();
        cout<<nameri[i]<<" : "<<nentriesRI[i]<<"\t";
    }
    cout<<endl;
    cout<<"\nAll RIs entries= " <<nentriesall<<endl;

    cout<<"\nPrinting first few timestamp"<<endl;
    for (Int_t i=0;i<nri;i++){
        for (Long64_t jentry=0;jentry<2;jentry++){
            chainRI[i]->GetEvent(jentry);
            cout<<"\t"<<nameri[i]<<"\t"<<decay->GetTimestamp();
        }

    }
    cout<<endl;
    cout<<endl;

    //!Read data here!
    //! Cut parameters
    Double_t dtioncut=600;
    Int_t sumexyrankcut=0;
    Int_t lightp_nzcut=6;
    Double_t neutronecut[2]={160,860};

    //! initilize parms for hisograming
    Int_t nbins[nri];
    Double_t maxval[nri];
    Double_t minval[nri];
    for (Int_t i=0;i<nri;i++){
        nbins[i]=2000;
        maxval[i]=10;
        minval[i]=-10;
    }

    //! initilize  hisograms
    TH1F* hdecay[nri];
    TH1F* hdecayp1n[nri];
    TH1F* hdecayp2n[nri];


    TH1F* hdecayp1nb[nri];
    TH1F* hdecayp2nb[nri];
    TH1F* hneuneuri[nri];
    TH1F* hneuneurib[nri];

    for (Int_t i=0;i<nri;i++){
        hdecay[i]=new TH1F(Form("hb%s",(char*)nameri[i].Data()),Form("hb%s",(char*)nameri[i].Data()),nbins[i],minval[i],maxval[i]);
        hdecayp1n[i]=new TH1F(Form("hb1n%s",(char*)nameri[i].Data()),Form("hb1n%s",(char*)nameri[i].Data()),nbins[i],minval[i],maxval[i]);
        hdecayp2n[i]=new TH1F(Form("hb2n%s",(char*)nameri[i].Data()),Form("hb2n%s",(char*)nameri[i].Data()),nbins[i],minval[i],maxval[i]);
        hdecayp1nb[i]=new TH1F(Form("hb1nb%s",(char*)nameri[i].Data()),Form("hb1nb%s",(char*)nameri[i].Data()),nbins[i],minval[i],maxval[i]);
        hdecayp2nb[i]=new TH1F(Form("hb2nb%s",(char*)nameri[i].Data()),Form("hb2nb%s",(char*)nameri[i].Data()),nbins[i],minval[i],maxval[i]);

        hneuneuri[i]=new TH1F(Form("hneuneuri%s",(char*)nameri[i].Data()),Form("hneuneuri%s",(char*)nameri[i].Data()),10000,-1000,1000);
        hneuneurib[i]=new TH1F(Form("hneuneurib%s",(char*)nameri[i].Data()),Form("hneuneurib%s",(char*)nameri[i].Data()),10000,-1000,1000);
    }

    //! 2 neutron correlation
    //!
    TH1F* hneuneu=new TH1F("hneuneu","tdiff",10000,-1000,1000);
    /*
    TH2F* h2neuneu1=new TH2F("h2neuneu1","e vs tdiff",10000,-1000,1000,500,0,1000);
    TH2F* h2neuneu2=new TH2F("h2neuneu2","pos vs tdiff",10000,-1000,1000,500,0,500);
    TH2F* h2neuneu3=new TH2F("h2neuneu3","ediff vs tdiff",10000,-1000,1000,1000,-1000,1000);
    TH2F* h2neuneu4=new TH2F("h2neuneu4","posdiff vs tdiff",10000,-1000,1000,500,0,1000);
    */


    TH1F* hneuneub=new TH1F("hneuneub","hneuneub",10000,-1000,1000);


    std::multimap < unsigned int,BELENHit* > neumap;
    std::multimap < unsigned int,BELENHit*> ::iterator neumap_it;
    std::multimap < unsigned int,BELENHit*> ::iterator neumap_it2;
    std::multimap < unsigned int,BELENHit* > neubmap;
    std::multimap < unsigned int,BELENHit*> ::iterator neubmap_it;
    std::multimap < unsigned int,BELENHit*> ::iterator neubmap_it2;

    //! fill histograms
    for (Int_t i=0;i<nri;i++){
        if (enablepid2[i]!=1) continue;
        for (Long64_t jentry=0;jentry<nentriesRI[i];jentry++){            
            //if (jentry>1) break;
            chainRI[i]->GetEvent(jentry);
            Double_t decaytime=((Double_t)((Long64_t)decay->GetTimestamp()-(Long64_t)imp->GetTimestamp()))/1e9;
            //! downstream veto cut
            Int_t ndownstreamveto=0;
            for (Int_t j=0;j<decay->GetAncMultipliticy();j++){
                if (decay->GetAncHit(j)->GetMyPrecious()==4) ndownstreamveto++;
            }
            //! beta cut
            if (decay->GetDtIon()<=dtioncut||decay->GetSumEXYRank()>sumexyrankcut||ndownstreamveto>0) continue;
            //if (decay->GetDtIon()<=dtioncut||ndownstreamveto>0) continue;
            //if (decay->GetDtIon()<=dtioncut) continue;
            //! neutron multipliciy counters and cut
            Int_t nneufreal=0;
            Int_t nneubreal=0;

            for (neumap_it=neumap.begin();neumap_it!=neumap.end();neumap_it++){
                BELENHit* hit1=(BELENHit*)neumap_it->second;
                delete hit1;
            }
            neumap.clear();
            for (neubmap_it=neubmap.begin();neubmap_it!=neubmap.end();neubmap_it++){
                BELENHit* hit1=(BELENHit*)neubmap_it->second;
                delete hit1;
            }
            neubmap.clear();

            for (Int_t j=0;j<decay->GetNeutronForwardMultipliticy();j++){
                if (decay->GetNeutronForwardHit(j)->GetEnergy()>neutronecut[0]&&decay->GetNeutronForwardHit(j)->GetEnergy()<neutronecut[1]&&decay->GetNeutronForwardHit(j)->GetFinalVetoTime()<0){
                    BELENHit* hitt=new BELENHit;
                    decay->GetNeutronForwardHit(j)->Copy(*hitt);
                    neumap.insert(make_pair(j,hitt));
                    nneufreal++;
                }
            }
            if (decaytime>0.06&&decaytime<2){
                for (neumap_it=neumap.begin();neumap_it!=neumap.end();neumap_it++){
                    unsigned int j1=neumap_it->first;
                    BELENHit* hit1=(BELENHit*)neumap_it->second;
                    for (neumap_it2=neumap.begin();neumap_it2!=neumap.end();neumap_it2++){
                        unsigned int j2=neumap_it2->first;
                        BELENHit* hit2=(BELENHit*)neumap_it2->second;
                        if (j2!=j1){
                            hneuneu->Fill(((Long64_t)hit2->GetTimestamp()-(Long64_t)hit1->GetTimestamp())/1e3);
                            hneuneuri[i]->Fill(((Long64_t)hit2->GetTimestamp()-(Long64_t)hit1->GetTimestamp())/1e3);
                            /*
                            h2neuneu1->Fill(((Long64_t)hit2->GetTimestamp()-(Long64_t)hit1->GetTimestamp())/1e3,hit1->GetEnergy());
                            h2neuneu2->Fill(((Long64_t)hit2->GetTimestamp()-(Long64_t)hit1->GetTimestamp())/1e3,sqrt(hit1->GetPosition().X()*hit1->GetPosition().X()+hit1->GetPosition().Y()*hit1->GetPosition().Y()));
                            h2neuneu3->Fill(((Long64_t)hit2->GetTimestamp()-(Long64_t)hit1->GetTimestamp())/1e3,hit2->GetEnergy()-hit1->GetEnergy());
                            h2neuneu4->Fill(((Long64_t)hit2->GetTimestamp()-(Long64_t)hit1->GetTimestamp())/1e3,sqrt((hit1->GetPosition().X()-hit2->GetPosition().X())*(hit1->GetPosition().X()-hit2->GetPosition().X())+(hit1->GetPosition().Y()-hit2->GetPosition().Y())*(hit1->GetPosition().Y()-hit2->GetPosition().Y())));
                            */
                        }
                    }
                }
            }else if (decaytime<-0.06&&decaytime>-2){
                for (neumap_it=neumap.begin();neumap_it!=neumap.end();neumap_it++){
                    unsigned int j1=neumap_it->first;
                    BELENHit* hit1=(BELENHit*)neumap_it->second;
                    for (neumap_it2=neumap.begin();neumap_it2!=neumap.end();neumap_it2++){
                        unsigned int j2=neumap_it2->first;
                        BELENHit* hit2=(BELENHit*)neumap_it2->second;
                        if (j2!=j1){
                            hneuneub->Fill(((Long64_t)hit2->GetTimestamp()-(Long64_t)hit1->GetTimestamp())/1e3);
                            hneuneurib[i]->Fill(((Long64_t)hit2->GetTimestamp()-(Long64_t)hit1->GetTimestamp())/1e3);
                            /*
                            h2neuneu1->Fill(((Long64_t)hit2->GetTimestamp()-(Long64_t)hit1->GetTimestamp())/1e3,hit1->GetEnergy());
                            h2neuneu2->Fill(((Long64_t)hit2->GetTimestamp()-(Long64_t)hit1->GetTimestamp())/1e3,sqrt(hit1->GetPosition().X()*hit1->GetPosition().X()+hit1->GetPosition().Y()*hit1->GetPosition().Y()));
                            h2neuneu3->Fill(((Long64_t)hit2->GetTimestamp()-(Long64_t)hit1->GetTimestamp())/1e3,hit2->GetEnergy()-hit1->GetEnergy());
                            h2neuneu4->Fill(((Long64_t)hit2->GetTimestamp()-(Long64_t)hit1->GetTimestamp())/1e3,sqrt((hit1->GetPosition().X()-hit2->GetPosition().X())*(hit1->GetPosition().X()-hit2->GetPosition().X())+(hit1->GetPosition().Y()-hit2->GetPosition().Y())*(hit1->GetPosition().Y()-hit2->GetPosition().Y())));
                            */
                        }
                    }
                }
            }


            for (Int_t j=0;j<decay->GetNeutronBackwardMultipliticy();j++){
                if (decay->GetNeutronBackwardHit(j)->GetEnergy()>neutronecut[0]&&decay->GetNeutronBackwardHit(j)->GetEnergy()<neutronecut[1]&&decay->GetNeutronBackwardHit(j)->GetFinalVetoTime()<0){
                    BELENHit* hitt=new BELENHit;
                    decay->GetNeutronBackwardHit(j)->Copy(*hitt);
                    neubmap.insert(make_pair(j,hitt));
                    nneubreal++;
                }
            }
            /*
            if (decaytime>0.06&&decaytime<1){
                for (neubmap_it=neubmap.begin();neubmap_it!=neubmap.end();neubmap_it++){
                    unsigned int j1=neubmap_it->first;
                    BELENHit* hit1=(BELENHit*)neubmap_it->second;
                    for (neubmap_it2=neubmap.begin();neubmap_it2!=neubmap.end();neubmap_it2++){
                        unsigned int j2=neubmap_it2->first;
                        BELENHit* hit2=(BELENHit*)neubmap_it2->second;
                        if (j2!=j1){
                            hneuneub->Fill(((Long64_t)hit2->GetTimestamp()-(Long64_t)hit1->GetTimestamp())/1e3);
                            hneuneurib[i]->Fill(((Long64_t)hit2->GetTimestamp()-(Long64_t)hit1->GetTimestamp())/1e3);
                        }
                    }
                }
            }
            */


            //! make histo
            hdecay[i]->Fill(decaytime);
            if (nneufreal==1) hdecayp1n[i]->Fill(decaytime);
            if (nneufreal==2) hdecayp2n[i]->Fill(decaytime);
            if (nneubreal==1) hdecayp1nb[i]->Fill(decaytime);
            if (nneubreal==2) hdecayp2nb[i]->Fill(decaytime);
        }
    }


    //! write hisograms to an output file
    TFile* fout=new TFile(outfile,"RECREATE");
    for (Int_t i=0;i<nri;i++){
        if (enablepid2[i]!=1) continue;
        hdecay[i]->Write();
        hdecayp1n[i]->Write();
        hdecayp2n[i]->Write();
        hdecayp1nb[i]->Write();
        hdecayp2nb[i]->Write();

        hneuneuri[i]->Write();
        hneuneurib[i]->Write();
    }
    hneuneu->Write();
    /*
    h2neuneu1->Write();
    h2neuneu2->Write();
    h2neuneu3->Write();
    h2neuneu4->Write();
    */
    hneuneub->Write();

    //! write deadtime
    TH2F* h2=new TH2F("deadtimeh2","deadtimeh2",chainall->GetNtrees(),0,chainall->GetNtrees(),500,0,100);
    Double_t deadtimec[6];
    getdeadtime(listfile,deadtimec,h2);
    TVectorD deadtimecontainer(7);
    for (Int_t i=0;i<6;i++) deadtimecontainer[i]=deadtimec[i];
    deadtimecontainer[6]=(Double_t)deadtimec[0]/(Double_t)deadtimec[1]*100;
    fout->cd();
    deadtimecontainer.Write("deadtime");
    h2->Write();
    fout->Close();
}


void makedecayx2histsim(char* infile, char* outfile){
    gROOT->ProcessLine(".L BELEN.h+");
    gROOT->ProcessLine(".L Clover.h+");
    gROOT->ProcessLine(".L AIDA.h+");
    gROOT->ProcessLine(".L AIDA.h+");
    gROOT->ProcessLine(".L AIDA.h+");
    gROOT->ProcessLine(".L AIDA.h+");
    gROOT->ProcessLine(".L DataStructNew.h+");
    TFile* f1=TFile::Open(infile);
    TTree* chainRI=(TTree*) f1->Get("tree");
    Int_t idx=0;
    IonBetaMult* imp=new IonBetaMult;
    IonBeta* decay=new IonBeta;
    //! set branch
    chainRI->SetBranchAddress("idx",&idx);
    chainRI->SetBranchAddress("ion",&imp);
    chainRI->SetBranchAddress("beta",&decay);
    Long64_t nentriesRI=chainRI->GetEntries();
    cout<<chainRI->GetEntries()<<endl;

    //! initilize  hisograms
    Int_t nbins=2000;
    Double_t maxval=10;
    Double_t minval=-10;
    TH1F* hdecay;
    TH1F* hdecayp1n;
    TH1F* hdecayp2n;
    TH1F* hdecayp1nb;
    TH1F* hdecayp2nb;
    hdecay=new TH1F("hbIn134","hbIn134",nbins,minval,maxval);
    hdecayp1n=new TH1F("hb1nIn134","hb1nIn134",nbins,minval,maxval);
    hdecayp2n=new TH1F("hb2nIn134","hb2nIn134",nbins,minval,maxval);
    hdecayp1nb=new TH1F("hb1nbIn134","hb1nbIn134",nbins,minval,maxval);
    hdecayp2nb=new TH1F("hb2nbIn134","hb2nbIn134",nbins,minval,maxval);

    for (Long64_t jentry=0;jentry<nentriesRI;jentry++){
        //if (jentry>1) break;
        chainRI->GetEvent(jentry);
        Int_t nneufreal=decay->GetNeutronForwardMultipliticy();
        Int_t nneubreal=decay->GetNeutronBackwardMultipliticy();
        //! make histo
        Double_t decaytime=((Double_t)((Long64_t)decay->GetTimestamp()-(Long64_t)imp->GetTimestamp()))/1e9;
        hdecay->Fill(decaytime);
        if (nneufreal==1) hdecayp1n->Fill(decaytime);
        if (nneufreal==2) hdecayp2n->Fill(decaytime);
        if (nneubreal==1) hdecayp1nb->Fill(decaytime);
        if (nneubreal==2) hdecayp2nb->Fill(decaytime);        
    }

    //! write hisograms to an output file
    TFile* fout=new TFile(outfile,"RECREATE");
    hdecay->Write();
    hdecayp1n->Write();
    hdecayp2n->Write();
    hdecayp1nb->Write();
    hdecayp2nb->Write();
    fout->Close();
}


void makedecaymlhtree(char* infile, char* outfile){
    gROOT->ProcessLine(".L BELEN.h+");
    gROOT->ProcessLine(".L Clover.h+");
    gROOT->ProcessLine(".L AIDA.h+");
    gROOT->ProcessLine(".L AIDA.h+");
    gROOT->ProcessLine(".L AIDA.h+");
    gROOT->ProcessLine(".L AIDA.h+");
    gROOT->ProcessLine(".L DataStructNew.h+");
    TFile* f1=TFile::Open(infile);
    TTree* chainRI=(TTree*) f1->Get("tree");
    Int_t idx=0;
    IonBetaMult* imp=new IonBetaMult;
    IonBeta* decay=new IonBeta;
    //! set branch
    chainRI->SetBranchAddress("idx",&idx);
    chainRI->SetBranchAddress("ion",&imp);
    chainRI->SetBranchAddress("beta",&decay);
    Long64_t nentriesRI=chainRI->GetEntries();
    cout<<chainRI->GetEntries()<<endl;

    //! initilize  hisograms
    Int_t nbins=2000;
    Double_t maxval=10;
    Double_t minval=-10;

    //! write hisograms to an output file
    TFile* fout=new TFile(outfile,"RECREATE");
    TTree* treebmlh = new TTree("treeb","treeb");
    TTree* treeb1nmlh = new TTree("treep1n","treep1n");
    TTree* treeb1nbmlh = new TTree("treep1nb","treep1nb");
    TTree* treeb2nmlh = new TTree("treep2n","treep2n");
    TTree* treeb2nbmlh = new TTree("treep2nb","treep2nb");

    Double_t bx;
    Double_t b1nx;
    Double_t b2nx;
    Double_t b1nbx;
    Double_t b2nbx;
    treebmlh->Branch("x",&bx,"x/D");
    treeb1nmlh->Branch("x",&b1nx,"x/D");
    treeb1nbmlh->Branch("x",&b1nbx,"x/D");
    treeb2nmlh->Branch("x",&b2nx,"x/D");
    treeb2nbmlh->Branch("x",&b2nbx,"x/D");

    TH1F* hdecay;
    TH1F* hdecayp1n;
    TH1F* hdecayp2n;
    TH1F* hdecayp1nb;
    TH1F* hdecayp2nb;
    hdecay=new TH1F("hbIn134","hbIn134",nbins,minval,maxval);
    hdecayp1n=new TH1F("hb1nIn134","hb1nIn134",nbins,minval,maxval);
    hdecayp2n=new TH1F("hb2nIn134","hb2nIn134",nbins,minval,maxval);
    hdecayp1nb=new TH1F("hb1nbIn134","hb1nbIn134",nbins,minval,maxval);
    hdecayp2nb=new TH1F("hb2nbIn134","hb2nbIn134",nbins,minval,maxval);

    for (Long64_t jentry=0;jentry<nentriesRI;jentry++){
        //if (jentry>1) break;
        chainRI->GetEvent(jentry);
        Int_t nneufreal=decay->GetNeutronForwardMultipliticy();
        Int_t nneubreal=decay->GetNeutronBackwardMultipliticy();
        //! make histo
        Double_t decaytime=((Double_t)((Long64_t)decay->GetTimestamp()-(Long64_t)imp->GetTimestamp()))/1e9;
        hdecay->Fill(decaytime);
        bx=decaytime;
        treebmlh->Fill();
        if (nneufreal==1) {
            hdecayp1n->Fill(decaytime);
            b1nx=decaytime;
            treeb1nmlh->Fill();
        }
        if (nneufreal==2) {
            hdecayp2n->Fill(decaytime);
            b2nx=decaytime;
            treeb2nmlh->Fill();
        }
        if (nneubreal==1) {
            hdecayp1nb->Fill(decaytime);
            b1nbx=decaytime;
            treeb1nbmlh->Fill();
        }
        if (nneubreal==2) {
            hdecayp2nb->Fill(decaytime);
            b2nbx=decaytime;
            treeb2nbmlh->Fill();
        }
    }

    treebmlh->Write();
    treeb1nmlh->Write();
    treeb2nmlh->Write();
    treeb1nbmlh->Write();
    treeb2nbmlh->Write();

    hdecay->Write();
    hdecayp1n->Write();
    hdecayp2n->Write();
    hdecayp1nb->Write();
    hdecayp2nb->Write();
    fout->Close();
}
