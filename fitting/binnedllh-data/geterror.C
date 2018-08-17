void geterror(Int_t parmsno, Int_t nbins,Double_t range,char* infile,char* riname,Int_t ishalflife=0){
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Sun Aug  5 17:12:46 2018 by ROOT version6.08/00)
//   from TTree treeSn135/treeSn135
//   found on file: fitresults/fitresultSn135lowin.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(infile);
   if (!f) {
      f = new TFile(infile);
   }
   TTree* tree;
   char treename[500];
   sprintf(treename,"tree%s",riname);
   f->GetObject(treename,tree);

//Declaration of leaves types
   Double_t        outparms[26];
   Double_t        outparmserr[26];
   Int_t           iparms[26];
   Double_t        neueff;
   Int_t           isvary[26];

   // Set branch addresses.
   tree->SetBranchAddress("outparms",outparms);
   tree->SetBranchAddress("outparmserr",outparmserr);
   tree->SetBranchAddress("iparms",iparms);
   tree->SetBranchAddress("neueff",&neueff);
   tree->SetBranchAddress("isvary",isvary);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// tree->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   Long64_t nentries = tree->GetEntries();
   tree->GetEntry(nentries-1);
   Double_t meanval;
   Double_t staterror;
   if (ishalflife==1) {
     meanval=log(2)/outparms[parmsno];
     staterror=log(2)/outparms[parmsno]/outparms[parmsno]*outparmserr[parmsno];
   }
   else {
     meanval=outparms[parmsno];   
     staterror=outparmserr[parmsno];
   }
   cout<<meanval<<endl;   
   TH1F* h1=new TH1F("h1","h1",nbins,meanval-range/2,meanval+range/2);

   Long64_t nbytes = 0;
   for (Long64_t i=0; i<nentries-1;i++) {
      nbytes += tree->GetEntry(i);
      if (ishalflife) h1->Fill(log(2)/outparms[parmsno]);
      else h1->Fill(outparms[parmsno]);
      //h1->Fill(outparms[parmsno]);
   }
   h1->Draw();
   Int_t binmean=h1->GetXaxis()->FindBin(meanval);
   TMarker meanmarker;   
   meanmarker.SetMarkerStyle(22);
   meanmarker.SetMarkerColor(2);
   meanmarker.SetMarkerSize(2);
   meanmarker.DrawMarker(meanval,h1->GetBinContent(binmean));

   Double_t nminus=0;
   for (Int_t i=binmean-1;i>0;i--) nminus+=h1->GetBinContent(i);
   Double_t nintegrate=0;
   Int_t binminus=0;
   Double_t minus=0;
   for (Int_t i=binmean-1;i>0;i--) {
     nintegrate+=h1->GetBinContent(i);
     if (nintegrate/nminus>0.682) {
       binminus=i;
       break;
     }
   }
   minus=h1->GetBinCenter(binminus);
   TMarker minusmarker;   
   minusmarker.SetMarkerStyle(23);
   minusmarker.SetMarkerColor(4);
   minusmarker.SetMarkerSize(2);
   minusmarker.DrawMarker(minus,h1->GetBinContent(binminus));
   Int_t nplus=0;
   for (Int_t i=binmean+1;i<h1->GetNbinsX();i++) nplus+=h1->GetBinContent(i);
   nintegrate=0;
   Int_t binplus=0;
   Double_t plus=0;
   for (Int_t i=binmean+1;i<h1->GetNbinsX();i++) {
     nintegrate+=h1->GetBinContent(i);
     if (nintegrate/nplus>0.682) {
       binplus=i;
       break;
     }
   }
   plus=h1->GetBinCenter(binplus);
   TMarker plusmarker;
   plusmarker.SetMarkerStyle(23);
   plusmarker.SetMarkerColor(4);
   plusmarker.SetMarkerSize(2);
   plusmarker.DrawMarker(plus,h1->GetBinContent(binplus));

   //! Gaussian fitting to the data
   h1->Fit("gaus","LR");
   TF1* f1=h1->GetFunction("gaus");
   cout<<f1->GetParameter(2)<<endl;
   std::ofstream ofs("resulttable.txt", std::ofstream::out | std::ofstream::app);
   ofs<<riname<<","<<parmsno<<","<<meanval<<","<<meanval-minus<<","<<plus-meanval<<","<<f1->GetParameter(2)<<","<<staterror<<endl;
}
