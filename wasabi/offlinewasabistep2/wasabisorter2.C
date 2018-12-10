#define tree_cxx
#include "tree.cpp"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <map>

#include "aidaclass.h"
#include "bigrips.cpp"



//! map bigrips
void TreeBranchBigRIPS(TTree* tree,bigrips* bigripsobj)
{
    tree->Branch("zet",&bigripsobj->zet,"zet/D");
    tree->Branch("aoq",&bigripsobj->aoq,"aoq/D");
    tree->Branch("beta",&bigripsobj->beta,"beta/D");
    tree->Branch("bripsts",&bigripsobj->ts,"bripsts/l");
}

void sorter2(char* infile,char* bripsfile, char* outfile, Int_t opt=0)
{
    //Long64_t fIonPidTWlow=10000;
    //Long64_t fIonPidTWup=10000;

    Long64_t fIonPidTWlow=10000;
    Long64_t fIonPidTWup=-7000;

//   In a ROOT session, you can do:
//      Root > .L tree.C
//      Root > tree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  //TFile* fin= TFile::Open(infile);


  //! input files
  TChain* group1ch = new TChain("group1");
  TChain* group2ch = new TChain("group2");

  std::ifstream ifs(infile);
  string filelist[1000];

  Int_t nfiles=0;
  while (!ifs.eof()){
      ifs>>filelist[nfiles];
      cout<<filelist[nfiles]<<endl;
      nfiles++;
  }
  nfiles=nfiles-1;
  cout<<"There are "<<nfiles<<" files in total!"<<endl;

  for (Int_t i=0;i<nfiles;i++){
      char tempchar[1000];
      sprintf(tempchar,"%s/group1",filelist[i].c_str());
      group1ch->Add(tempchar);

      sprintf(tempchar,"%s/group2",filelist[i].c_str());
      group2ch->Add(tempchar);
  }

  tree* inbeta=new tree(group1ch);
  inbeta->SetID(5);
  tree* inion=new tree(group2ch);
  inion->SetID(4);


  //! set configuration file
  inbeta->ReadConfigTable("dssdconfiglowe.txt");
  inion->ReadConfigTable("dssdconfighighe.txt");

   Long64_t nentriesbeta = inbeta->fChain->GetEntries();
   Long64_t nentriesion = inion->fChain->GetEntries();

   if (opt==1) nentriesbeta=0;
   cout<<nentriesbeta<<"-"<<nentriesion<<endl;



   //! read bigrips file
   std::multimap<Long64_t,Long64_t> fbigripsMap;
   std::multimap<Long64_t,Long64_t>::iterator fbigripsMap_it;

   bigrips* brips=new bigrips(bripsfile);
   brips->GetTsEntry(fbigripsMap);
   cout<<"read "<<fbigripsMap.size()<<" entries for BigRIPS TS table"<<endl;

   //! output
   //AIDAClass* wasabi=new AIDAClass;
   AIDAClass* wasabi=new AIDAClass;
   TFile* fout = new TFile(outfile,"RECREATE");
   TTree* tree = new TTree("wasabi","wasabi");
   tree->Branch("wasabi",&wasabi);

   //! special for implantation events
   tree->Branch("dssd_e",inion->dssd_e,Form("dssd_e[%d]/D",N_DSSD*(N_STRIP_X+N_STRIP_Y)));
   tree->Branch("dssd_bl",inion->dssd_bl,Form("dssd_bl[%d]/D",N_DSSD*(N_STRIP_X+N_STRIP_Y)));
   tree->Branch("dssd_ch",inion->dssd_ch,Form("dssd_ch[%d]/I",N_DSSD*(N_STRIP_X+N_STRIP_Y)));
   tree->Branch("dssd_EX",inion->dssd_EX,Form("dssd_EX[%d]/D",N_DSSD));
   tree->Branch("dssd_EY",inion->dssd_EY,Form("dssd_EY[%d]/D",N_DSSD));
   tree->Branch("dssd_NX",inion->dssd_NX,Form("dssd_NX[%d]/I",N_DSSD));
   tree->Branch("dssd_NY",inion->dssd_NY,Form("dssd_NY[%d]/I",N_DSSD));
   tree->Branch("dssd_No",inion->dssd_No,Form("dssd_No[%d]/I",N_DSSD));


   TreeBranchBigRIPS(tree,brips);

   std::multimap <unsigned long long,AIDAClass*> datamap; //! sort by timestamp
   std::multimap <unsigned long long,AIDAClass*>::iterator datamap_it; //! sort by timestamp

   Long64_t nbytes = 0, nb = 0;

   //! ion entries
   for (Long64_t jentry=0; jentry<nentriesion;jentry++) {
       Long64_t ientry = inion->LoadTree(jentry);
       if (ientry < 0) break;
       nb = inion->fChain->GetEntry(jentry);   nbytes += nb;
       inion->Reconstruction();
       //! get nz
       //Int_t nz[N_DSSD];
       //for (int ii=0;ii<N_DSSD;ii++) nz[ii]=0;
       Int_t nz=0;
       Int_t maxz=-1;
       for (unsigned short i=0;i<inion->GetNRecoData();i++){
           if (inion->GetRecoData(i)->EX>0&&inion->GetRecoData(i)->EY>0){
             //nz[(Int_t)inion->GetRecoData(i)->z]++;
               nz++;
               maxz=(Int_t)inion->GetRecoData(i)->z;
           }
       }
       for (unsigned short i=0;i<inion->GetNRecoData();i++){
           if (inion->GetRecoData(i)->EX>0&&inion->GetRecoData(i)->EY>0){
             if (maxz==(Int_t)inion->GetRecoData(i)->z){
                   AIDAClass* data=new AIDAClass;
                   syncrecodata(data,inion->GetRecoData(i));
                   //data->nz=nz[(Int_t)data->z];
                   data->nz=nz;
                   data->evt=jentry;
                   datamap.insert(std::make_pair(data->T,data));
             }
           }
       }
       inion->ClearRecoData();
   }

   //! beta entries
   for (Long64_t jentry=0; jentry<nentriesbeta;jentry++) {
       Long64_t ientry = inbeta->LoadTree(jentry);
       if (ientry < 0) break;
       nb = inbeta->fChain->GetEntry(jentry);   nbytes += nb;
       inbeta->Reconstruction();

       //! get nz
       //Int_t nz[N_DSSD];
       //for (int ii=0;ii<N_DSSD;ii++) nz[ii]=0;
       Int_t nz=0;
       for (unsigned short i=0;i<inbeta->GetNRecoData();i++){
           if (inbeta->GetRecoData(i)->EX>0&&inbeta->GetRecoData(i)->EY>0){
             //nz[(Int_t)inbeta->GetRecoData(i)->z]++;
             nz++;
             //cout<<nz<<endl;
           }
       }
       for (unsigned short i=0;i<inbeta->GetNRecoData();i++){
           if (inbeta->GetRecoData(i)->EX>0&&inbeta->GetRecoData(i)->EY>0){
               AIDAClass* data=new AIDAClass;
               syncrecodata(data,inbeta->GetRecoData(i));
               //data->nz=nz[(Int_t)data->z];
               data->nz=nz;
               datamap.insert(std::make_pair(data->T,data));
           }
       }
       inbeta->ClearRecoData();
   }

   //! fill time-ordered data to tree
   for(datamap_it=datamap.begin();datamap_it!=datamap.end();datamap_it++){
       AIDAClass* hit=(AIDAClass*)datamap_it->second;
       syncrecodata(wasabi,hit);
       Int_t ncorr=0;

       if (hit->ID==4){
           //! Correlate imp with bigrips
           Long64_t ts=(Long64_t)datamap_it->first;
           Long64_t ts1 = (Long64_t)ts - fIonPidTWlow;
           Long64_t ts2 = (Long64_t)ts + fIonPidTWup;
           Long64_t corrts = 0;
           Long64_t correntry = 0;
           Long64_t check_time = 0;
           fbigripsMap_it = fbigripsMap.lower_bound(ts1);
           while(fbigripsMap_it!=fbigripsMap.end()&&fbigripsMap_it->first<ts2){
               corrts = (Long64_t) fbigripsMap_it->first;
               correntry = (Long64_t) fbigripsMap_it->second;
               if (corrts!=check_time){
                   check_time=corrts;
                   //! fill data here
                   brips->GetEntry(correntry);
                   ncorr++;
                   break;
               }
               fbigripsMap_it++;
           }
       }

       if (!ncorr){
           brips->zet=-9999;
           brips->aoq=-9999;
           brips->beta=-9999;
           brips->ts=0;
       }

       //! reconstruction again
       inion->fChain->GetEntry(hit->evt);
       inion->Reconstruction();
       tree->Fill();
       inion->ClearRecoData();
   }
   tree->Write();
   fout->Close();

}
