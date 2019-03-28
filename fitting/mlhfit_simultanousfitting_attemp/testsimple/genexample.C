void genexample(){
    Double_t hl=0.128;
    //Double_t hl=0.5;
    Double_t p1n=50;
    Double_t p2n=25;
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

        pratio=r3->Rndm()*100;
        if (pratio<80) {
            tp1n=r3->Rndm()*22-11;
            treep1n->Fill();
        }
	pratio=r3->Rndm()*100;
        if (pratio<70) {
            tp2n=r3->Rndm()*22-11;
            treep2n->Fill();
        }

	pratio=r3->Rndm()*100;
        if (pratio<p1n) {
            tp1n=tall;
            treep1n->Fill();
        }else if (pratio>=p1n&&pratio<p1n+p2n){
	    tp2n=tall;
            treep2n->Fill();
	}

    }
    hdecayc->Write();
    treeb->Write();
    treebb->Write();
    treep1n->Write();
    treep2n->Write();
    fout->Close();
}

