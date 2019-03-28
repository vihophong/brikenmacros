
#include "TTree.h"
#include "TFile.h"

#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"

#include "TCut.h"
#include "TCutG.h"
#include "TLatex.h"
#include "TSpectrum.h"
#include "TF1.h"

#include <map>

#include "fstream"

#include "gammaclass.h"

typedef struct{
    Double_t gc_E;
    Double_t gc_T;//gamma time in ns
    Int_t gc_ch;
} gammahit;

typedef struct{
    Double_t ab_E;
    Double_t ab_T;//gamma time in ns
    Int_t ab_ch;//first hit channel
    Short_t ab_mult[4];//multiplicity
} gammaab;


void addbackcalibcheck(char* infile)
{


    std::multimap < unsigned long long, gammahit* >  gmap;

    gammaclass gin(infile);

    Long64_t nentries=gin.fChain->GetEntries();

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries/4;jentry++) {
        if (jentry%100000==0) cout<<jentry<<"\t/\t"<<nentries<<endl;
       Long64_t ientry = gin.LoadTree(jentry);
       if (ientry < 0) break;
       nb = gin.fChain->GetEntry(jentry);   nbytes += nb;

       gammahit* gc=new gammahit();
       gc->gc_T=gin.fts;
       gc->gc_ch=gin.fid;
       gc->gc_E=gin.fen;
       gmap.insert(make_pair(gin.fts,gc));
    }

    cout<<"finished step1"<<endl;



    TH1F* h1=new TH1F("h1","h1",2000,-10000,10000);

    std::multimap < unsigned long long, gammahit* > gmap2;

    std::multimap < unsigned long long, gammahit* >::iterator gmap_it;
    std::multimap < unsigned long long, gammahit* >::iterator gmap2_it;

    gmap2=gmap;

    Int_t k=0;
    Int_t nwithcorr=0;
    for (gmap2_it=gmap2.begin();gmap2_it!=gmap2.end();gmap2_it++){
        if (k%10000==0) cout<<k<<"/"<<gmap2.size()<<" ncorr="<<nwithcorr<<endl;

        unsigned long long ts=gmap2_it->first;
        Long64_t ts1 = (Long64_t)ts - (Long64_t)10000;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)10000;
        Long64_t corrts = 0;
        Int_t ncorr=0;
        Long64_t check_time = 0;
        gmap_it = gmap.lower_bound(ts1);
        while(gmap_it!=gmap.end()&&gmap_it->first<ts2){
            corrts = (Long64_t) gmap_it->first;
            //gammahit* gc = (gammahit*) gmap_it->second;
            //if (corrts!=check_time){
            check_time=corrts;
            if (corrts!=ts){

                h1->Fill(corrts-(Long64_t)ts);

                ncorr++;
            }
            gmap_it++;
        }

        if (ncorr>0) nwithcorr++;

        k++;
    }
    h1->SaveAs("test.root");

}

void addbackcalib(char* infile,char* outfile)
{
    //! gamma calibration
    Double_t fsep[8] = {4.912531E+05,4.998953E+05,5.248885E+05,7.622378E+05,1.108068E+06,1.147523E+06,1.164587E+06,9.553438E+05}; //separation points
    Double_t flow_offset[8] = {0.310056,0.602804	,0.594689,0.00575245,-0.409776,-0.739931,-0.231043,-0.410475};//low offset
    Double_t flow_gain[8] = {0.00158189,0.00154124,0.00153239,0.000668465,0.000707824,0.000687903,0.000672509,0.000700449};//low gain
    Double_t flow_se[8]= {2.5638600000E-12,5.1029500000E-12,8.0936100000E-12,-2.5440800000E-12,-3.2340000000E-12,-4.4067300000E-12,-3.9413800000E-12,-2.6372600000E-12};// low second order
    Double_t fhigh_offset[8] = {1.32609,4.65127,-17.4998,4.21305,-1.68140e+00,-1.13557,-2.79797,6.44019};//high offset
    Double_t fhigh_gain[8] = {0.00158113,0.00153239,0.00159751,0.000659638,7.05388e-04,0.000682542,0.000670503,0.000688148};//high gain
    Double_t fhigh_se[8] = {-9.9227100000E-14,6.6060100000E-12,-5.0294000000E-11,1.7949000000E-12,1.00018e-16,5.6552400000E-13,-3.2623400000E-13,2.7326500000E-12};//high second order
    Double_t fcgainold[8] =
    {0.0015836325,0.0015439776,0.0015329461,0.000664755,0.00070,0.000683305,0.0006685615,0.0006968099};
    Double_t fcoffsetold[8] =
    {-0.6017825233,-0.2629124127,0.2931096096,-0.0427639699,-0.5577366054,-0.6254548603,-0.1966707051,-0.2212697416};

    TH1F* habe=new TH1F("habe","habe",2000,0,2000);
    TH1F* hgce=new TH1F("hgce","hgce",2000,0,2000);
    TH1F* habt=new TH1F("habt","habt",5000,0,5000);

    std::multimap < unsigned long long, gammahit* >  gmap;
    std::multimap < unsigned long long, gammahit* >::iterator gmap_it;


    std::multimap < unsigned long long, gammahit* >  gmap2;
    std::multimap < unsigned long long, gammahit* >::iterator gmap2_it;


    gammaclass gin(infile);

    Long64_t nentries=gin.fChain->GetEntries();

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        if (jentry%100000==0) cout<<jentry<<"\t/\t"<<nentries<<endl;
       Long64_t ientry = gin.LoadTree(jentry);
       if (ientry < 0) break;
       nb = gin.fChain->GetEntry(jentry);   nbytes += nb;

       if (gin.fid<5){//clover 1
           gammahit* gc=new gammahit();
           gc->gc_T=gin.fts;
           gc->gc_ch=gin.fid;


           gc->gc_E=gin.fen;
           //convert back to adc
           gc->gc_E=(gc->gc_E-fcoffsetold[gin.fid-1])/fcgainold[gin.fid-1];
           //apply new calibration
           if (gc->gc_E<fsep[gin.fid-1]){//low energy calibration
               gc->gc_E=flow_offset[gin.fid-1]+flow_gain[gin.fid-1]*gc->gc_E+flow_se[gin.fid-1]*gc->gc_E*gc->gc_E;
           }else{//high energy calibration
               gc->gc_E=fhigh_offset[gin.fid-1]+fhigh_gain[gin.fid-1]*gc->gc_E+fhigh_se[gin.fid-1]*gc->gc_E*gc->gc_E;
           }

           hgce->Fill(gc->gc_E);
           gmap.insert(make_pair(gin.fts,gc));
       }else{//clover 2
           gammahit* gc=new gammahit();
           gc->gc_T=gin.fts;
           gc->gc_ch=gin.fid-4;
           gc->gc_E=gin.fen;
           //convert back to adc
           gc->gc_E=(gc->gc_E-fcoffsetold[gin.fid-1])/fcgainold[gin.fid-1];
           //apply new calibration
           if (gc->gc_E<fsep[gin.fid-1]){//low energy calibration
               gc->gc_E=flow_offset[gin.fid-1]+flow_gain[gin.fid-1]*gc->gc_E+flow_se[gin.fid-1]*gc->gc_E*gc->gc_E;
           }else{//high energy calibration
               gc->gc_E=fhigh_offset[gin.fid-1]+fhigh_gain[gin.fid-1]*gc->gc_E+fhigh_se[gin.fid-1]*gc->gc_E*gc->gc_E;
           }
           gmap2.insert(make_pair(gin.fts,gc));
       }

    }

    cout<<"finished step1"<<gmap.size()<<" / "<<gmap2.size()<<endl;


    gammaab* abdata=new gammaab();
    gammaab* abdata2=new gammaab();

    gammaab* gcdata=new gammaab();
    gammaab* gcdata2=new gammaab();

    TFile* f1=new TFile(outfile,"recreate");
    TTree* tree1 = new TTree("addback1","addback1");
    tree1->Branch("ab",&abdata);
    TTree* tree2 = new TTree("addback2","addback2");
    tree2->Branch("ab",&abdata2);

    TTree* tree3 = new TTree("gc1","gc1");
    tree3->Branch("gc",&gcdata);
    TTree* tree4 = new TTree("gc2","gc2");
    tree4->Branch("gc",&gcdata2);


    Bool_t ab_started=false;
    Long64_t ab_beg=0;
    Long64_t ab_end=0;
    Long64_t ab_window=1000;
    Long64_t ts_prev=0;

    Int_t k=0;
    Int_t ch_beg=0;
    Double_t esum=0;


    memset(abdata->ab_mult,0,sizeof(abdata->ab_mult));

    for (gmap_it=gmap.begin();gmap_it!=gmap.end();gmap_it++){
        if (k%10000==0) cout<<k<<"/"<<gmap.size()<<endl;
        Long64_t ts=gmap_it->first;
        gammahit* gc=gmap_it->second;


        gcdata->ab_ch=gc->gc_ch;

        gcdata->ab_E=gc->gc_E;
        gcdata->ab_T=gc->gc_T;
        tree3->Fill();


        if (!ab_started) {
            ab_beg=ts;
            ab_started=true;
        }

        if (ch_beg==0) ch_beg=gc->gc_ch;

        if (ts-ab_beg>ab_window&&ab_started){// end of event, next event start
            ab_end=ts_prev;
            //--------event operation here
            habt->Fill(ab_end-ab_beg);
            habe->Fill(esum);

            abdata->ab_ch=ch_beg;
            abdata->ab_E=esum;
            abdata->ab_T=ab_beg;
            tree1->Fill();

            ch_beg=gc->gc_ch;
            ab_beg=ts;// start new event
            esum=0;
            memset(abdata->ab_mult,0,sizeof(abdata->ab_mult));
        }


        //" colleting hits to an event here
        esum+=gc->gc_E;
        if (gc->gc_ch==1) abdata->ab_mult[0]++;
        else if (gc->gc_ch==2) abdata->ab_mult[1]++;
        else if (gc->gc_ch==3) abdata->ab_mult[2]++;
        else if (gc->gc_ch==4) abdata->ab_mult[3]++;

        ts_prev=ts;
        k++;
    }


    ab_started=false;
    ab_beg=0;
    ab_end=0;
    ab_window=1000;
    ts_prev=0;

    k=0;
    ch_beg=0;
    esum=0;


    memset(abdata2->ab_mult,0,sizeof(abdata2->ab_mult));

    for (gmap2_it=gmap2.begin();gmap2_it!=gmap2.end();gmap2_it++){
        if (k%10000==0) cout<<k<<"/"<<gmap2.size()<<endl;
        Long64_t ts=gmap2_it->first;
        gammahit* gc=gmap2_it->second;

        gcdata2->ab_ch=gc->gc_ch;
        gcdata2->ab_E=gc->gc_E;
        gcdata2->ab_T=gc->gc_T;
        tree4->Fill();


        if (!ab_started) {
            ab_beg=ts;
            ab_started=true;
        }

        if (ch_beg==0) ch_beg=gc->gc_ch;

        if (ts-ab_beg>ab_window&&ab_started){// end of event, next event start
            ab_end=ts_prev;
            //--------event operation here

            abdata2->ab_ch=ch_beg;
            abdata2->ab_E=esum;
            abdata2->ab_T=ab_beg;
            tree2->Fill();

            ch_beg=gc->gc_ch;
            ab_beg=ts;// start new event
            esum=0;
            memset(abdata2->ab_mult,0,sizeof(abdata2->ab_mult));
        }


        //" colleting hits to an event here
        esum+=gc->gc_E;
        if (gc->gc_ch==1) abdata2->ab_mult[0]++;
        else if (gc->gc_ch==2) abdata2->ab_mult[1]++;
        else if (gc->gc_ch==3) abdata2->ab_mult[2]++;
        else if (gc->gc_ch==4) abdata2->ab_mult[3]++;

        ts_prev=ts;
        k++;
    }

    tree1->Write();
    tree2->Write();
    tree3->Write();
    tree4->Write();
    habt->Write();
    habe->Write();
    hgce->Write();
    f1->Close();
}
