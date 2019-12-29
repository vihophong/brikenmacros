#define MERGER_H

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

#include "fstream"
#include "map"

#include "treebr.h"
#include "treem0_wasabi.h"
#include "treem1_wasabi.h"
#include "treem2_wasabi.h"
#include "treem3_wasabi.h"


//DSSDs constants
#define N_DSSD 2
#define N_STRIP_X 16
#define N_STRIP_Y 16

//DE constants
#define N_DE 18

//CLOVER constants
#define N_GE_CRYSTAL 8

using namespace std;

treebr* tree_bigrips;
treem0_wasabi* treew0;
treem1_wasabi* treew1;
treem2_wasabi* treew2;
treem3_wasabi* treew3;

const Long64_t timewindow_plus=10000;
const Long64_t timewindow_minus=10000;

const Long64_t tswidth_offsetsearch = 600;
const Long64_t print_low_offsetsearch = 32.17685e12;
const Long64_t print_high_offsetsearch = 32.17686e12;

Long64_t tsoffset_b_minus_w=32.17685e12+4.29222e+06;


void mergem2(char* infilebrips, char* infilewas,char* outfilestr){
    tree_bigrips=new treebr(infilebrips);
    treew0=new treem0_wasabi(infilewas);
    treew1=new treem1_wasabi(infilewas);
    treew2=new treem2_wasabi(infilewas);
    treew3=new treem3_wasabi(infilewas);

    //! output file
    TFile* outfile=new TFile(outfilestr,"RECREATE");
    TH1F* htscorr=new TH1F("htscorr","tscorr",500,-timewindow_minus,timewindow_plus);
    TTree* mergedtree = new TTree("tree","tree");
    mergedtree->Branch("bigrips",tree_bigrips);
    mergedtree->Branch("wasabi",treew2);

    //! fill merged map
    std::multimap < unsigned long long, unsigned int> fbigripsMap;
    std::multimap < unsigned long long, unsigned int>::iterator fbigripsMap_it;

    std::multimap < unsigned long long, unsigned int> fwasabiMap;
    std::multimap < unsigned long long, unsigned int>::iterator fwasabiMap_it;

    cout<<"BIGRIPS entries = "<<tree_bigrips->fChain->GetEntries()<<endl;
    cout<<"WASABI entries = "<<treew2->fChain->GetEntries()<<endl;

    for (Long64_t jentry=0; jentry<tree_bigrips->fChain->GetEntries();jentry++) {
        tree_bigrips->GetEntry(jentry);
        fbigripsMap.insert(make_pair(tree_bigrips->ts-tsoffset_b_minus_w,jentry));

    }
    cout<<tree_bigrips->ts<<endl;

    for (Long64_t jentry=0; jentry<treew2->fChain->GetEntries();jentry++) {
        treew2->GetEntry(jentry);
        fwasabiMap.insert(make_pair(treew2->ts*8,jentry));

    }
    cout<<treew2->ts*8<<endl;

    for (fbigripsMap_it=fbigripsMap.begin();fbigripsMap_it!=fbigripsMap.end();fbigripsMap_it++){
        Long64_t ts=fbigripsMap_it->first;
        Long64_t entry=fbigripsMap_it->second;
        tree_bigrips->GetEntry(entry);


        //! Correlate f11 with neutron
        Long64_t ts1 = (Long64_t)ts - timewindow_minus;
        Long64_t ts2 = (Long64_t)ts + timewindow_plus;
        Long64_t corrts = 0;
        Long64_t correntry = 0;
        Int_t ncorr=0;
        fwasabiMap_it = fwasabiMap.lower_bound(ts1);
        while(fwasabiMap_it!=fwasabiMap.end()&&fwasabiMap_it->first<ts2){
            corrts = (Long64_t) fwasabiMap_it->first;
            correntry = (Long64_t) fwasabiMap_it->second;
            htscorr->Fill(corrts - ts);
            treew2->GetEntry(correntry);
            mergedtree->Fill();
            ncorr++;
            fwasabiMap_it++;
        }

    }

    htscorr->Write();
    mergedtree->Write();
    outfile->Close();    

}

void mergem1m2(char* infilebrips, char* infilewas,char* outfilestr){
    tree_bigrips=new treebr(infilebrips);
    treew0=new treem0_wasabi(infilewas);
    treew1=new treem1_wasabi(infilewas);
    treew2=new treem2_wasabi(infilewas);
    treew3=new treem3_wasabi(infilewas);

    //! output file
    TFile* outfile=new TFile(outfilestr,"RECREATE");
    TH1F* htscorr=new TH1F("htscorr","tscorr",500,-timewindow_minus,timewindow_plus);
    TTree* mergedtree = new TTree("tree","tree");
    mergedtree->Branch("bigrips",tree_bigrips);
    mergedtree->Branch("wasabim1",treew1);
    mergedtree->Branch("wasabim2",treew2);

    //! fill merged map
    std::multimap < unsigned long long, unsigned int> fbigripsMap;
    std::multimap < unsigned long long, unsigned int>::iterator fbigripsMap_it;

    std::multimap < unsigned long long, unsigned int> fwasabiMap;
    std::multimap < unsigned long long, unsigned int>::iterator fwasabiMap_it;

    cout<<"BIGRIPS entries = "<<tree_bigrips->fChain->GetEntries()<<endl;
    cout<<"WASABI M1 entries = "<<treew1->fChain->GetEntries()<<endl;
    cout<<"WASABI M2 entries = "<<treew2->fChain->GetEntries()<<endl;


    for (Long64_t jentry=0; jentry<tree_bigrips->fChain->GetEntries();jentry++) {
        tree_bigrips->GetEntry(jentry);
        if (jentry==0) cout<<"tsbegin_br = "<<tree_bigrips->ts<<endl;
        fbigripsMap.insert(make_pair(tree_bigrips->ts,jentry));
    }
    cout<<"tsend_br = "<<tree_bigrips->ts<<endl;

    for (Long64_t jentry=0; jentry<treew1->fChain->GetEntries();jentry++) {
        treew1->GetEntry(jentry);
        if (jentry==0) cout<<"tsbegin_w1 = "<<treew1->ts*8<<endl;
    }
    cout<<"tsend_w1 = "<<treew1->ts*8<<endl;

    for (Long64_t jentry=0; jentry<treew2->fChain->GetEntries();jentry++) {
        treew2->GetEntry(jentry);
        if (jentry==0) cout<<"tsbegin_w2 = "<<treew2->ts*8<<endl;
        fwasabiMap.insert(make_pair(treew2->ts*8,jentry));
    }
    cout<<"tsend_w2 = "<<treew2->ts*8<<endl;


    for (fbigripsMap_it=fbigripsMap.begin();fbigripsMap_it!=fbigripsMap.end();fbigripsMap_it++){
        Long64_t ts=fbigripsMap_it->first;
        Long64_t entry=fbigripsMap_it->second;
        tree_bigrips->GetEntry(entry);


        //! Correlate f11 with neutron
        Long64_t ts1 = (Long64_t)ts - timewindow_minus;
        Long64_t ts2 = (Long64_t)ts + timewindow_plus;
        Long64_t corrts = 0;
        Long64_t correntry = 0;
        Int_t ncorr=0;
        fwasabiMap_it = fwasabiMap.lower_bound(ts1);
        while(fwasabiMap_it!=fwasabiMap.end()&&fwasabiMap_it->first<ts2){
            corrts = (Long64_t) fwasabiMap_it->first;
            correntry = (Long64_t) fwasabiMap_it->second;
            htscorr->Fill(corrts - ts);
            treew2->GetEntry(correntry);
            mergedtree->Fill();
            ncorr++;
            fwasabiMap_it++;
        }

    }

    htscorr->Write();
    mergedtree->Write();
    outfile->Close();

}

void findtimestampoffset(char* infilebrips, char* infilewas,char* outfilestr){
    tree_bigrips=new treebr(infilebrips);
    treew2=new treem2_wasabi(infilewas);

    std::multimap < unsigned long long, pair<Long64_t,unsigned long long> > fbigripsMap;
    std::multimap < unsigned long long, pair<Long64_t,unsigned long long> >::iterator fbigripsMap_it;

    std::multimap < unsigned long long, pair<Long64_t,unsigned long long> > fwasabiMap;
    std::multimap < unsigned long long, pair<Long64_t,unsigned long long> >::iterator fwasabiMap_it;

    TFile* outfile=new TFile(outfilestr,"RECREATE");
    TH1F* hcorr = new TH1F("hcorr","hcorr",200,-tswidth_offsetsearch-tswidth_offsetsearch*0.1,tswidth_offsetsearch+tswidth_offsetsearch*0.1);

    long long w2_runtime = 0;
    for (Long64_t jentry=0; jentry<treew2->fChain->GetEntries();jentry++) {
        treew2->GetEntry(jentry);
        if (jentry==0) {
            cout<<"tsbegin_w2 = "<<treew2->ts*8<<endl;
            w2_runtime = treew2->ts*8;
        }
    }
    w2_runtime = treew2->ts*8-w2_runtime;
    cout<<"tsend_w2 = "<<treew2->ts*8<<endl;
    cout<<"runtime_w2 = "<<(double)w2_runtime/1e9<<" s"<<endl;


    long long br_runtime = 0;

    long long br_begints = 0;
    Long64_t nentriesbr = 0;
    for (Long64_t jentry=0; jentry<tree_bigrips->fChain->GetEntries();jentry++) {
        tree_bigrips->GetEntry(jentry);

        nentriesbr = jentry;
        if (jentry==0) {
            cout<<"tsbegin_br = "<<tree_bigrips->ts<<endl;
            br_begints = tree_bigrips->ts;
            br_runtime = tree_bigrips->ts;
        }else{
            if ((tree_bigrips->ts-br_begints)>w2_runtime) break;
        }
    }
    cout<<"tsend_br = "<<tree_bigrips->ts<<endl;
    br_runtime = tree_bigrips->ts - br_runtime;
    cout<<"runtime_br = "<<(double)br_runtime/1e9<<" s"<<endl;
    cout<<"last entry = "<<nentriesbr<<endl;

    //! FILL SEARCHING MAP

    long long tsw_prev=0;
    for (Long64_t jentry=0; jentry<treew2->fChain->GetEntries();jentry++) {
        treew2->GetEntry(jentry);
        long long tsw_now=treew2->ts*8;

        if (tsw_prev>0){
            long long tsw_distance = tsw_now-tsw_prev;
            fwasabiMap.insert(make_pair(tsw_distance,make_pair(jentry,tsw_now)));
        }
        tsw_prev=tsw_now;
    }
    long long tsb_prev = 0;
    for (Long64_t jentry=0; jentry<nentriesbr;jentry++) {
        tree_bigrips->GetEntry(jentry);
        long long tsb_now = tree_bigrips->ts;
        if (tsb_prev>0){
            long long tsb_distance = tsb_now - tsb_prev;
            fbigripsMap.insert(make_pair(tsb_distance,make_pair(jentry,tsb_now)));
        }
        tsb_prev = tsb_now;
    }

    for (fbigripsMap_it=fbigripsMap.begin();fbigripsMap_it!=fbigripsMap.end();fbigripsMap_it++){
        Long64_t tsb_distance=fbigripsMap_it->first;
        Long64_t jentryb=fbigripsMap_it->second.first;
        Long64_t tsb_now=fbigripsMap_it->second.second;
        //cout<<tsb_distance<<"\t"<<jentryb<<"\t"<<tsb_now<<"\t"<<endl;
    }

    for (fwasabiMap_it=fwasabiMap.begin();fwasabiMap_it!=fwasabiMap.end();fwasabiMap_it++){
        Long64_t tsw_distance=fwasabiMap_it->first;
        Long64_t jentryw=fwasabiMap_it->second.first;
        Long64_t tsw_now=fwasabiMap_it->second.second;
        //cout<<tsw_distance<<"\t"<<jentryw<<"\t"<<tsw_now<<"\t"<<endl;
    }

    //! SEARCH for OFFSET
    //!
    TTree* tree=new TTree("tree","tree");
    Long64_t offsetval = 0;
    Long64_t widthval = 0;
    tree->Branch("o",&offsetval,"o/L");
    tree->Branch("w",&widthval,"w/L");

    for (fbigripsMap_it=fbigripsMap.begin();fbigripsMap_it!=fbigripsMap.end();fbigripsMap_it++){
        Long64_t tsb_distance=fbigripsMap_it->first;
        Long64_t jentryb=fbigripsMap_it->second.first;
        Long64_t tsb_now=fbigripsMap_it->second.second;

        //! Correlate f11 with neutron
        Long64_t ts1 = (Long64_t)tsb_distance - tswidth_offsetsearch;
        Long64_t ts2 = (Long64_t)tsb_distance + tswidth_offsetsearch;
        Long64_t tsw_distance=0;
        Long64_t jentryw=0;
        Long64_t tsw_now=0;

        Int_t ncorr=0;
        fwasabiMap_it = fwasabiMap.lower_bound(ts1);
        while(fwasabiMap_it!=fwasabiMap.end()&&fwasabiMap_it->first<ts2){
            tsw_distance=fwasabiMap_it->first;
            jentryw=fwasabiMap_it->second.first;
            tsw_now=fwasabiMap_it->second.second;
            widthval = tsw_distance - tsb_distance;
            offsetval = tsb_now - tsw_now;
            hcorr->Fill(widthval);
            tree->Fill();

            //! print val
            if (offsetval>print_low_offsetsearch&&offsetval<print_high_offsetsearch){
                cout<<offsetval<<endl;
            }
            ncorr++;
            fwasabiMap_it++;
        }
    }

    hcorr->Write();
    tree->Write();
    outfile->Close();
}

void merger(char* infilebrips, char* infilewas,char* outfilestr)
{
    mergem2(infilebrips,infilewas,outfilestr);
    //findtimestampoffset(infilebrips,infilewas,outfilestr);
}
