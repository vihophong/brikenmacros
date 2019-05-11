void ana(){
  ofstream outfile("output.txt");
  TH1F *temphis;
	Int_t nRun=10;
	Int_t nGroup=6;
	Int_t k=0;	
   for (Int_t i=0;i<nRun;i++){
      TFile* h=new TFile(Form("Hist%d.root",i));
			TTree* BRIKEN=(TTree*) h->Get("BRIKEN");
//Get Energy 
			BRIKEN->Draw(Form("partEnergy>>histE%d(100,0,1000)",i),"E>180&&E<800","goff");
			temphistE=(TH1F*) h->Get(Form("histE%d",i));
			outfile<<temphistE->GetMean()<<" ";
//GetTotal fire on All
			BRIKEN->Draw(Form("E>>hist%d(100,0,1000)",i),"E>180&&E<800","goff");
			temphist=(TH1F*) h->Get(Form("hist%d",i));
			outfile<<temphist->GetEntries()<<" ";
//Get fire on each group
			for (Int_t j=0;j<nGroup;j++){
      BRIKEN->Draw(Form("E>>his%d(100,0,1000)",k),Form("E>180&&E<800&&detGroup==%d",j),"goff");
      temphis=(TH1F*) h->Get(Form("his%d",k));
      outfile<<temphis->GetEntries()<<" ";
			k++;
			}
			outfile<<"\n";
      h->Close("R");
   }
   outfile.close();
}
