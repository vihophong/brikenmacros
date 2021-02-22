{
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Tue Feb 16 15:45:07 2021 by ROOT version6.08/00)
//   from TTree treeneuveto/treeneuveto
//   found on file: decaynew_lowin/decay_brips3107.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("decaynew_lowin/decay_brips3107.root");
   if (!f) {
      f = new TFile("decaynew_lowin/decay_brips3107.root");
   }
    f->GetObject("treeneuveto",tree);

//Declaration of leaves types
   Int_t           id;
   Long64_t        ts;

   // Set branch addresses.
   treeneuveto->SetBranchAddress("id",&id);
   treeneuveto->SetBranchAddress("ts",&ts);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// treeneuveto->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   Long64_t nentries = treeneuveto->GetEntries();

   Long64_t nbytes = 0;
//   for (Long64_t i=0; i<nentries;i++) {
//      nbytes += treeneuveto->GetEntry(i);
//   }
}
