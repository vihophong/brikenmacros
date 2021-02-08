#define makefinalhistos_cxx
#include "makefinalhistos.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <fstream>
#include "TMath.h"

double Median(const TH1F * h1) {

   int n = h1->GetXaxis()->GetNbins();
   double x[10000];
   double y[10000];
   for (int i=0;i<n;i++){
       x[i]=h1->GetBinCenter(i+1);
       y[i]=h1->GetBinContent(i+1);
   }
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, x, y);
}

void makefinalhistos::MakeFinalHisto(char* outfile, Double_t decaytmin=0.05, Double_t decaytmax=1, Int_t layer=-1, Double_t beta=0)
{


//   In a ROOT session, you can do:
//      root> .L makefinalhistos.C
//      root> makefinalhistos t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
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
   if (fChain == 0) return;

   //! stuff for hit disrtibution
   Double_t fHe3Id2posX[MaxID];
   Double_t fHe3Id2posY[MaxID];
   Double_t fHe3Id2posZ[MaxID];
   Double_t fHe3Id2diameter[MaxID];
   UShort_t fHe3Id2ring[MaxID];
   Double_t fHe3Id2length[MaxID];

   Double_t fCrystalId2posX[MaxIndex1][MaxIndex2];
   Double_t fCrystalId2posY[MaxIndex1][MaxIndex2];
   Double_t fCrystalId2posZ[MaxIndex1][MaxIndex2];

   //! read geo mapping
   std::ifstream inpf("He3_mapping_wecal.txt");
   if (inpf.fail()){
       cout<<"No BELEN Mapping file is given"<<endl;
       return;
   }
   cout<<"Start reading BELEN Mapping file: "<<endl;
   Int_t id,index1,index2;
   UShort_t ring;
   Double_t x,y,z;
   Double_t d,length;
   Double_t ftoadc,ecal1;
   Int_t mm=0;

   while (inpf.good()){
       inpf>>id>>index1>>index2>>d>>x>>y>>z>>ring>>length>>ftoadc>>ecal1;
       if (id<=500){//for he3
           fHe3Id2posX[id]=x;
           fHe3Id2posY[id]=y;
           fHe3Id2posZ[id]=z;
           fHe3Id2diameter[id]=d;
           fHe3Id2ring[id]=ring;
           fHe3Id2length[id]=length;
       }else if(id>500){ //for clover
           fCrystalId2posX[index1][index2]=x;
           fCrystalId2posY[index1][index2]=y;
           fCrystalId2posZ[index1][index2]=z;
       }
       //cout<<He3id<<"-"<<daqId<<"-"<<d<<"-"<<x<<"-"<<y<<"-"<<z<<endl;
       mm++;
   }
   cout<<"Read "<<mm<<" line"<<endl;
   inpf.close();


   //! output files
   TFile* ofile=new TFile(outfile,"recreate");
   ofile->cd();


   //! output histograms
   TH1F* hdecay=new TH1F("hdecay","hdecay",binning,lowlimit,uplimit);
   TH1F* hdecay1n=new TH1F("hdecay1n","hdecay1n",binning,lowlimit,uplimit);
   TH1F* hdecay2n=new TH1F("hdecay2n","hdecay2n",binning,lowlimit,uplimit);

   TH1F* hdecay1nbwd=new TH1F("hdecay1nbwd","hdecay1nbwd",binning,lowlimit,uplimit);
   TH1F* hdecaygt0nbwd=new TH1F("hdecaygt0nbwd","hdecaygt0nbwd",binning,lowlimit,uplimit);
   TH1F* hdecay2nbwd=new TH1F("hdecay2nbwd","hdecay2nbwd",binning,lowlimit,uplimit);



   TH1F* hibgfw=new TH1F("hibgfw","forward ionbeta gamma spec",8000,0,4000);
   TH1F* hibgbw=new TH1F("hibgbw","backward ionbeta gamma spec",8000,0,4000);

   TH1F* hibgfwnfw=new TH1F("hibgfwnfw","hibgfw forward ionbeta gamma spec with forward neutron ",8000,0,4000);
   TH1F* hibgbwnfw=new TH1F("hibgbwnbw","hibgfw backward ionbeta gamma spec with forward neutron",8000,0,4000);

   TH1F* hibgfwnbw=new TH1F("hibgfwnbw","hibgfw forward ionbeta gamma spec with backward neutron ",8000,0,4000);
   TH1F* hibgbwnbw=new TH1F("hibgbwnfw","hibgfw backward ionbeta gamma spec with backward neutron",8000,0,4000);

   TH1F* h1=new TH1F("hdistr","hit distribution",200,0,50);
   TH1F* h2=new TH1F("hibn","hibn",300,-10,20);
   TH1F* h3=new TH1F("hibngate","hibngate",300,-10,20);
   TH1F* h4=new TH1F("hdistrbwib","hit distribution backward time ion beta",200,0,50);
   TH1F* h5=new TH1F("hdistrbwbn","hit distribution backward time beta neutron",200,0,50);
   TH1F* h6=new TH1F("hdistrbwbnbwib","hit distribution backward time beta neutron of backward ionbeta",200,0,50);

   TH2F* h2z=new TH2F("h2z","Zdistributionvsdecay",binning,lowlimit,uplimit,6,0,6);

   TH1F* hbeta=new TH1F("hbeta","Beta distribution",500,0.6,0.7);

   //!

   Int_t ncountouterring=0;
   Int_t ncountinnerring=0;
   Int_t ncountouterringbwib=0;
   Int_t ncountinnerringbwib=0;
   Int_t ncountouterringbwbn=0;
   Int_t ncountinnerringbwbn=0;
   Int_t ncountouterringbwbnbwib=0;
   Int_t ncountinnerringbwbnbwib=0;


   Int_t ncountouterring1=0;
   Int_t ncountinnerring1=0;
   Int_t ncountouterringbwib1=0;
   Int_t ncountinnerringbwib1=0;
   Int_t ncountouterringbwbn1=0;
   Int_t ncountinnerringbwbn1=0;
   Int_t ncountouterringbwbnbwib1=0;
   Int_t ncountinnerringbwbnbwib1=0;

   Int_t ncountouterring2=0;
   Int_t ncountinnerring2=0;
   Int_t ncountouterringbwib2=0;
   Int_t ncountinnerringbwib2=0;
   Int_t ncountouterringbwbn2=0;
   Int_t ncountinnerringbwbn2=0;
   Int_t ncountouterringbwbnbwib2=0;
   Int_t ncountinnerringbwbnbwib2=0;



   //! output tree for MLH method
   TTree* otree = new TTree("tree","tree") ;
   Double_t* px = new Double_t ;
   Int_t* py = new Int_t ;
   Int_t* pz = new Int_t ;
   Int_t* pb = new Int_t ;
   otree->Branch("x",px,"x/D") ;
   otree->Branch("y",py,"y/I") ;


   //! output tree for backward neutrons
   TTree* otreebw = new TTree("treebw","treebw") ;
   Double_t* pxb = new Double_t ;
   Int_t* pyb = new Int_t ;
   otreebw->Branch("x",pxb,"x/D") ;
   otreebw->Branch("y",pyb,"y/I") ;

   Long64_t nentries = fChain->GetEntries();

   cout<<nentries<<endl;
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      hbeta->Fill(decay_beta);

      //! select layer if a non-negative value is specified
      if (layer>=0){
        if (decay_z!=layer) continue;
      }

      //! select beta cut
      if (beta!=0){
          //cout<<"beta"<<beta<<endl;
          if (beta<0){
              if (decay_beta>-beta) continue; //select velocity < input value
          }else{
              if (decay_beta<beta) continue; //select velocity > input value
          }
      }

      //! fill histos
      hdecay->Fill(decay_t);
      h2z->Fill(decay_t,decay_z);
      if (neu_hit==1) hdecay1n->Fill(decay_t);
      if (neu_hit==2) hdecay2n->Fill(decay_t);

      if (neub_hit==1) hdecay1nbwd->Fill(decay_t);
      if (neub_hit==2) hdecay2nbwd->Fill(decay_t);
      if (neub_hit>0) hdecaygt0nbwd->Fill(decay_t);


      //! fill tree
      //if (neu_hit<3){
      *px=decay_t;
      *pz=decay_z;
      *py=neu_hit;
      *pb=neub_hit;
      otree->Fill();


      *pxb=decay_t;
      *pyb=neub_hit;
      otreebw->Fill();
      //}
      //! get hit distributions
      //!selecting all multiplicity
      for (Int_t i=0;i<neu_hit;i++){
          h2->Fill(decay_t);
          if (decay_t>decaytmin&&decay_t<decaytmax) {
              Double_t xx=fHe3Id2posX[neu_ch[i]]/10;//in cm
              Double_t yy=fHe3Id2posY[neu_ch[i]]/10;//in cm
              h1->Fill(sqrt(xx*xx+yy*yy));
              h3->Fill(decay_t);
              if (sqrt(xx*xx+yy*yy)<outerinnersep){
                  ncountinnerring++;
                  if (neu_hit==1) ncountinnerring1++;
                  if (neu_hit==2) ncountinnerring2++;
              }else{
                  ncountouterring++;
                  if (neu_hit==1) ncountouterring1++;
                  if (neu_hit==2) ncountouterring2++;
              }
          }
      }


      //! backward ion beta hits
      for (Int_t i=0;i<neu_hit;i++){
          if (decay_t<-decaytmin&&decay_t>-decaytmax) {
              Double_t xx=fHe3Id2posX[neu_ch[i]]/10;//in cm
              Double_t yy=fHe3Id2posY[neu_ch[i]]/10;//in cm
              h4->Fill(sqrt(xx*xx+yy*yy));
              if (sqrt(xx*xx+yy*yy)<outerinnersep){
                  ncountinnerringbwib++;
                  if (neu_hit==1) ncountinnerringbwib1++;
                  if (neu_hit==2)ncountinnerringbwib2++;
              }else{
                  ncountouterringbwib++;
                  if (neu_hit==1) ncountouterringbwib1++;
                  if (neu_hit==2) ncountouterringbwib2++;
              }
          }
      }

      //! backward beta neutron hits
      for (Int_t i=0;i<neub_hit;i++){
          h2->Fill(decay_t);
          if (decay_t>decaytmin&&decay_t<decaytmax) {
              Double_t xx=fHe3Id2posX[neub_ch[i]]/10;//in cm
              Double_t yy=fHe3Id2posY[neub_ch[i]]/10;//in cm
              h5->Fill(sqrt(xx*xx+yy*yy));
              if (sqrt(xx*xx+yy*yy)<outerinnersep){
                  ncountinnerringbwbn++;
                  if (neub_hit==1) ncountinnerringbwbn1++;
                  if (neub_hit==2) ncountinnerringbwbn2++;
              }else{
                  ncountouterringbwbn++;
                  if (neub_hit==1) ncountouterringbwbn1++;
                  if (neub_hit==2) ncountouterringbwbn2++;
              }
          }
      }
      //! backward beta neutron hits of backward ion beta
      for (Int_t i=0;i<neub_hit;i++){
          h2->Fill(decay_t);
          if (decay_t<-decaytmin&&decay_t>-decaytmax) {
              Double_t xx=fHe3Id2posX[neub_ch[i]]/10;//in cm
              Double_t yy=fHe3Id2posY[neub_ch[i]]/10;//in cm
              h6->Fill(sqrt(xx*xx+yy*yy));
              if (sqrt(xx*xx+yy*yy)<outerinnersep){
                  ncountinnerringbwbnbwib++;
                  if (neub_hit==1) ncountinnerringbwbnbwib1++;
                  if (neub_hit==2) ncountinnerringbwbnbwib2++;
              }else{
                  ncountouterringbwbnbwib++;
                  if (neub_hit==1) ncountouterringbwbnbwib1++;
                  if (neub_hit==2) ncountouterringbwbnbwib2++;
              }
          }
      }

      //!forward ionbeta gamma spec
      if (decay_t>decaytmin&&decay_t<decaytmax) {
          for (Int_t i=0;i<gc_hit;i++){
                  hibgfw->Fill(gc_E[i]);
                  if (neu_hit>0) hibgfwnfw->Fill(gc_E[i]);
                  if (neub_hit>0) hibgfwnbw->Fill(gc_E[i]);
          }
      }
      //!backward ionbeta gamma spec
      if (decay_t<-decaytmin&&decay_t>-decaytmax) {
          for (Int_t i=0;i<gc_hit;i++){
                  hibgbw->Fill(gc_E[i]);
                  if (neu_hit>0) hibgbwnfw->Fill(gc_E[i]);
                  if (neub_hit>0) hibgbwnbw->Fill(gc_E[i]);
          }
      }
   }

   std::ofstream str("outimpandratio.txt",std::ios::app);
   ncountinnerring=ncountinnerring-ncountinnerringbwib-ncountinnerringbwbn+ncountinnerringbwbnbwib;
   ncountouterring=ncountouterring-ncountouterringbwib-ncountouterringbwbn+ncountouterringbwbnbwib;

   ncountinnerring1=ncountinnerring1-ncountinnerringbwib1-ncountinnerringbwbn1+ncountinnerringbwbnbwib1;
   ncountouterring1=ncountouterring1-ncountouterringbwib1-ncountouterringbwbn1+ncountouterringbwbnbwib1;

   ncountinnerring2=ncountinnerring2-ncountinnerringbwib2-ncountinnerringbwbn2+ncountinnerringbwbnbwib2;
   ncountouterring2=ncountouterring2-ncountouterringbwib2-ncountouterringbwbn2+ncountouterringbwbnbwib2;

   str<<nameri<<"\t"<<ncountinnerring<<"\t"<<ncountouterring<<"\t"<<ncountinnerring1<<"\t"<<ncountouterring1<<"\t"<<ncountinnerring2<<"\t"<<ncountouterring2<<"\t"<<fChainImp->GetEntries()<<endl;
   cout<<"inner = "<<ncountinnerring<<" outer = "<<ncountouterring<<endl;
   cout<<"inner 1n = "<<ncountinnerring1<<" outer 1n = "<<ncountouterring1<<endl;
   cout<<"inner 2n = "<<ncountinnerring2<<" outer 2n = "<<ncountouterring2<<endl;
   str.close();
   std::ofstream str2("outbetacut.txt",std::ios::app);
   str2<<nameri<<"\t"<<Median(hbeta)<<endl;

   //! write tree and histograms
   hdecay->Write();
   hdecay1n->Write();
   hdecay2n->Write();

   hdecay1nbwd->Write();
   hdecay2nbwd->Write();
   hdecaygt0nbwd->Write();

   h1->Write();
   h2->Write();
   h3->Write();
   h4->Write();
   h5->Write();
   h6->Write();

   TH1F* hibgreal=(TH1F*)hibgfw->Clone();
   hibgreal->SetName("hibgreal");
   hibgreal->SetTitle(" ionbeta gamma spec real(subtreacted)");
   hibgreal->Add(hibgbw,-1);

   TH1F* hibngreal=(TH1F*)hibgfwnfw->Clone();
   hibngreal->SetName("hibngreal");
   hibngreal->SetTitle(" ionbeta gamma-neutron spec real(subtreacted)");
   hibngreal->Add(hibgbwnfw,-1);
   hibngreal->Add(hibgfwnbw,-1);
   hibngreal->Add(hibgbwnbw);

   hibgreal->Write();
   hibngreal->Write();

   hibgfw->Write();
   hibgbw->Write();
   hibgfwnfw->Write();
   hibgbwnfw->Write();
   hibgfwnbw->Write();
   hibgbwnbw->Write();
   h2z->Write();
   hbeta->Write();

   otree->Write();
   otreebw->Write();

   ofile->Close();

   delete px;
   delete py;
   delete pz;
   delete pb;
}
