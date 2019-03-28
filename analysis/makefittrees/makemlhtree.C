#define makemlhtree_cxx
#include "makemlhtree.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void makemlhtree::Loop(char* outfile, Int_t binning, Int_t layer)
{
    Double_t lowlimit=-10;
    Double_t uplimit=20;

//   In a ROOT session, you can do:
//      root> .L makemlhtree.C
//      root> makemlhtree t
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

   TFile* ofile=new TFile(outfile,"recreate");
   ofile->cd();

   TH1F* hdecay=new TH1F("hdecay","hdecay",binning,lowlimit,uplimit);
   TH1F* hdecay1n=new TH1F("hdecay1n","hdecay1n",binning,lowlimit,uplimit);
   TH1F* hdecay2n=new TH1F("hdecay2n","hdecay2n",binning,lowlimit,uplimit);

   TH1F* hdecay1nbwd=new TH1F("hdecay1nbwd","hdecay1nbwd",binning,lowlimit,uplimit);
   TH1F* hdecaygt0nbwd=new TH1F("hdecaygt0nbwd","hdecaygt0nbwd",binning,lowlimit,uplimit);
   TH1F* hdecay2nbwd=new TH1F("hdecay2nbwd","hdecay2nbwd",binning,lowlimit,uplimit);

   TTree* otree = new TTree("tree","tree") ;
   Double_t* px = new Double_t ;
   Int_t* py = new Int_t ;
   Int_t* pz = new Int_t ;
   Int_t* pb = new Int_t ;
   otree->Branch("x",px,"x/D") ;
   otree->Branch("y",py,"y/I") ;
   otree->Branch("z",pz,"z/I") ;
   otree->Branch("b",pb,"b/I") ;

   Long64_t nentries = fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      *px=decay_t;
      *pz=decay_z;

      
      //cout<<"aaa"<<decay_z<<endl;
      //if (decay_ey>100){
      if (layer>=0){
	if (decay_z!=layer) continue;
      }
      hdecay->Fill(decay_t);
      if (neu_hit==1) hdecay1n->Fill(decay_t);
      if (neu_hit==2) hdecay2n->Fill(decay_t);

      if (neub_hit==1) hdecay1nbwd->Fill(decay_t);
      if (neub_hit==2) hdecay2nbwd->Fill(decay_t);
      if (neub_hit>0) hdecaygt0nbwd->Fill(decay_t);
      if (neu_hit<3){
          *py=neu_hit;
          *pb=neub_hit;
          otree->Fill();
      }

      //}
   }
   TTree* otreeimp = new TTree("treeimp","treeimp") ;
   otreeimp->Branch("z",pz,"z/I") ;
   Long64_t nentriesImp = fChainImp->GetEntries();
   cout<<nentriesImp<<endl;
   for (Long64_t jentry=0; jentry<nentriesImp;jentry++) {
       fChainImp->GetEntry(jentry);
       *pz=fz;
       otreeimp->Fill();
   }

   hdecay->Write();
   hdecay1n->Write();
   hdecay2n->Write();

   hdecay1nbwd->Write();
   hdecay2nbwd->Write();
   hdecaygt0nbwd->Write();


   otree->Write();
   otreeimp->Write();
   ofile->Close();
   delete px;
   delete py;
   delete pz;
   delete pb;
}
void makemlhtree::PlotNeuHitPattern(char* outfile,Double_t decaytmin,Double_t decaytmax)
{

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
    std::ifstream inpf("He3_mapping.txt");
    if (inpf.fail()){
        cout<<"No BELEN Mapping file is given"<<endl;
        return;
    }
    cout<<"Start reading BELEN Mapping file: "<<endl;
    Int_t id,index1,index2;
    UShort_t ring;
    Double_t x,y,z;
    Double_t d,length;
    Int_t mm=0;

    while (inpf.good()){
        inpf>>id>>index1>>index2>>d>>x>>y>>z>>ring>>length;
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

    if (fChain == 0) return;

    TFile* ofile=new TFile(outfile,"recreate");
    ofile->cd();

    TH1F* h1=new TH1F("h1","h1",200,0,2000);
    TH1F* h2=new TH1F("h2","h2",200,-10,10);
    TH1F* h3=new TH1F("h3","h3",200,-10,10);

    Long64_t nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       // if (Cut(ientry) < 0) continue;

       for (Int_t i=0;i<neu_hit;i++){
           if (neu_hit<5){
               h2->Fill(decay_t);
               Double_t xx=fHe3Id2posX[neu_ch[i]]/10;
               Double_t yy=fHe3Id2posY[neu_ch[i]]/10;
               //cout<<xx<<"\t"<<yy<<endl;
               if (decay_t>decaytmin&&decay_t<decaytmax) {
		 //cout<<xx<<"\t"<<yy<<endl;
		 h1->Fill(xx*xx+yy*yy);
                   h3->Fill(decay_t);
               }
           }
       }
    }

    h1->Write();
    h2->Write();
    h3->Write();
    ofile->Close();
}
