
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TThread.h"

#include "dataStruct.h"


TCondition* bigripsFinished;
TCondition* implantFinished;
TCondition* belenFinished;

TThread *t[4];


AIDASimpleStruct* faida = NULL;
BELENHit* fneutron = NULL;
CloverHit* fclover = NULL;
BELENHit* fanc = NULL;
TreeData* fbigrips = NULL;

TTree* ftrAIDA = NULL;
TTree* ftrBigrips = NULL;
TTree* ftrNeutron = NULL;
TTree* ftrGamma = NULL;
TTree* ftrAnc = NULL;

Long64_t fnentriesAIDA = 0;
Long64_t fnentriesNeutron = 0;
Long64_t fnentriesGamma = 0;
Long64_t fnentriesAnc = 0;
Long64_t fnentriesBigrips = 0;


Long64_t fTsNow;
Long64_t fTsMax[4];

std::multimap < unsigned long long, TreeData* > fbigripsMap;
std::multimap < unsigned long long, TreeData* >::iterator fbigripsMap_it;

std::multimap < unsigned long long, std::pair<TreeData*, AIDASimpleStruct* > > fimplantMap;
std::multimap < unsigned long long, std::pair<TreeData*, AIDASimpleStruct* > >::iterator fimplantMap_it;

std::multimap < unsigned long long, BELENHit* > fneutronMap;
std::multimap < unsigned long long, BELENHit* >::iterator fneutronMap_it;
std::multimap < unsigned long long, BELENHit* > fancMap;
std::multimap < unsigned long long, BELENHit* >::iterator fancMap_it;
std::multimap < unsigned long long, CloverHit* > fcloverMap;
std::multimap < unsigned long long, CloverHit* >::iterator fcloverMap_it;



void *readBigrips(void *ptr){
    cout<<"Reading "<<fnentriesBigrips<<" bigrips items"<<endl;
    cout<<"Printing first few timestamp:"<<endl;
    for (Long64_t i=0;i<fnentriesBigrips;i++){
        ftrBigrips->GetEvent(i);
        TreeData* data = new TreeData;
        data->ts = fbigrips->ts;
        data->sts = fbigrips->sts;
        data->tof = fbigrips->ts;
        data->zet = fbigrips->ts;
        data->aoq = fbigrips->ts;
        data->f5x = fbigrips->ts;
        data->f11x = fbigrips->ts;
        data->f11y = fbigrips->ts;
        data->f11dt = fbigrips->ts;
        data->beta = fbigrips->ts;
        fbigripsMap.insert(make_pair(data->ts,data));
        if (i<10) cout<<"bigripsts "<<fbigrips->ts<<endl;
    }
    cout<<"Finished reading Bigrips"<<endl;
    bigripsFinished->Broadcast();
    return 0;
}

void *readAidaImplant(void *ptr){
    bigripsFinished->Wait();
    Long64_t fIonPidTWup = 0;
    Long64_t fIonPidTWlow = 20000;
    cout<<"Reading "<<fnentriesAIDA<<" enetries in AIDA tree"<<endl;
    cout<<"Printing first few timestamp:"<<endl;
    cout<<"Bigrip size = "<<fbigripsMap.size()<<endl;
    for (Long64_t i=0;i<fnentriesAIDA;i++){
        ftrAIDA->GetEvent(i);
        if (faida->GetID()==4){
            //! Correlate imp with bigrips
            Long64_t ts1 = (Long64_t)faida->GetTimestamp() - (Long64_t)fIonPidTWlow;
            Long64_t ts2 = (Long64_t)faida->GetTimestamp() + (Long64_t)fIonPidTWup;
            Long64_t corrts = 0;
            Int_t ncorr=0;
            Long64_t check_time = 0;
            fbigripsMap_it = fbigripsMap.lower_bound(ts1);
            while(fbigripsMap_it!=fbigripsMap.end()&&fbigripsMap_it->first<ts2){
                corrts = (Long64_t) fbigripsMap_it->first;
                if (corrts!=check_time){
                    check_time=corrts;
                    TreeData* correntry = (TreeData*) fbigripsMap_it->second;
                    AIDASimpleStruct* data =new AIDASimpleStruct();
                    faida->Copy(*data);
                    fimplantMap.insert(make_pair(faida->GetTimestamp(),make_pair(correntry,data)));
                    ncorr++;
                    break;
                }
                fbigripsMap_it++;
            }
        }
    }
    cout<<"Finished reading Implant"<<endl;
    return 0;
}

void *readBriken(void *ptr){
    cout<<"Reading "<<fnentriesNeutron<<" neutrons, "
       <<fnentriesGamma<<" gammas and "<<fnentriesAnc<<" anc hits in BELEN tree"<<endl;
    for (Long64_t i=0;i<fnentriesNeutron;i++){
        ftrNeutron->GetEvent(i);
        BELENHit* data =new BELENHit();
        fneutron->Copy(*data);
        fneutronMap.insert(make_pair(fneutron->GetTimestamp(),data));
    }
    for (Long64_t i=0;i<fnentriesGamma;i++){
        ftrGamma->GetEvent(i);
        CloverHit* data =new CloverHit();
        fclover->Copy(*data);
        fcloverMap.insert(make_pair(fclover->GetTimestamp(),data));
    }
    for (Long64_t i=0;i<fnentriesAnc;i++){
        ftrAnc->GetEvent(i);
        BELENHit* data =new BELENHit();
        fanc->Copy(*data);
        fancMap.insert(make_pair(fanc->GetTimestamp(),data));
    }
    cout<<"Finished reading BELEN"<<endl;
    return 0;
}

void *doMerge(void *ptr){
    cout<<"Merging data:"<<endl;
    for (Long64_t i=0;i<fnentriesAIDA;i++){
        ftrAIDA->GetEvent(i);
        
        if (faida->GetID()==5){
            Long64_t ts = faida->GetTimestamp();
            short betaz=faida->GetHitPositionZ();//! correct dZ
            int betaminx=faida->GetMinHitPositionX();
            int betamaxx=faida->GetMaxHitPositionX();
            int betaminy=faida->GetMinHitPositionY();
            int betamaxy=faida->GetMaxHitPositionY();
            fTsNow = faida->GetTimestamp();
            //cout<<"ftsnow = "<<fTsNow<<endl;
            
            //! correlation start
            Long64_t ts1 = (Long64_t)ts - 10e9;
            Long64_t ts2 = (Long64_t)ts + 10e9;
            Long64_t corrts = 0;
            Long64_t check_time = 0;
            
            fimplantMap_it = fimplantMap.lower_bound(ts1);

            while(fimplantMap_it!=fimplantMap.end()&&fimplantMap_it->first<ts2){
                corrts = (Long64_t) fimplantMap_it->first;
                AIDASimpleStruct* imp = (AIDASimpleStruct*) fimplantMap_it->second.second;
                
                short impz= imp->GetHitPositionZ()+imp->GetDZ(); // corrected dZ
                
                int impminx=imp->GetMinHitPositionX();
                int impmaxx=imp->GetMaxHitPositionX();
                int impminy=imp->GetMinHitPositionY();
                int impmaxy=imp->GetMaxHitPositionY();
                
                if (corrts!=check_time&&betaz==impz){// avoid multiple filling corrts!=check_time
                    //! implant-beta spatial correlation
                    //! new using overlaping area
                    if (!( ((impminx<=betamaxx&&impminx>=betaminx)||(impmaxx<=betamaxx&&impmaxx>=betaminx))&&
                        ((impminy<=betamaxy&&impminy>=betaminy)||(impmaxy<=betamaxy&&impmaxy>=betaminy)) )) {
                        fimplantMap_it++;
                        continue;
                    }
                    
                    //! Fill correlation
                    //cout<<"merge "<< corrts << "\t" << faida->GetTimestamp()<<endl;
                    
                    
                    check_time=corrts;
                }
                
                fimplantMap_it++;
            }
            //! correlation end
            
        }
    }
    return 0;
}

void mergertest(char* finputAida, char* finputBriken, char* finputBigrips)
{
    //gROOT->ProcessLine(".L dataStruct.h+");
    //! Initialization and Read First Entries
    
    fTsNow = 0;
    for (Int_t i=0;i<4;i++) fTsMax[i] = 0;
    
    if (finputAida!=NULL){
        TFile* fAidaFile = new TFile(finputAida);
        fAidaFile->GetObject("aida",ftrAIDA);
        ftrAIDA->SetBranchAddress("aida",&faida);
        fnentriesAIDA = ftrAIDA->GetEntries();
    }

    if (finputBriken!=NULL){
        //! init briken
        TFile*  fBrikenFile = new TFile(finputBriken);
        fBrikenFile->GetObject("neutron",ftrNeutron);
        fBrikenFile->GetObject("gamma",ftrGamma);
        fBrikenFile->GetObject("anc",ftrAnc);
        ftrNeutron->SetBranchAddress("neutron",&fneutron);
        ftrGamma->SetBranchAddress("gamma",&fclover);
        ftrAnc->SetBranchAddress("anc",&fanc);
        fnentriesNeutron = ftrNeutron->GetEntries();
        fnentriesGamma = ftrGamma->GetEntries();
        fnentriesAnc = ftrAnc->GetEntries();
    }

    if (finputBigrips!=NULL){
        //! init bigrips
        TFile* fBigripsFile = new TFile(finputBigrips);
        fBigripsFile->GetObject("tree",ftrBigrips);
        ftrBigrips->SetBranchAddress("bigrips",&fbigrips);
        fnentriesBigrips = ftrBigrips->GetEntries();
    }
    
    bigripsFinished = new TCondition(0);
    implantFinished = new TCondition(0);
    belenFinished = new TCondition(0);

    t[0] = new TThread("t0",readAidaImplant,(void*) 0);
    t[1] = new TThread("t1",readBigrips,(void*) 1);
    t[2] = new TThread("t2",readBriken,(void*) 2);
    
    t[0]->Run();
    t[1]->Run();
    t[2]->Run();
    
    t[0]->Join();
    t[1]->Join();
    t[2]->Join();
    
    cout<<"Merging data:"<<endl;
    for (Long64_t i=0;i<fnentriesAIDA;i++){
        ftrAIDA->GetEvent(i);
        
        if (faida->GetID()==5){
            Long64_t ts = faida->GetTimestamp();
            short betaz=faida->GetHitPositionZ();//! correct dZ
            int betaminx=faida->GetMinHitPositionX();
            int betamaxx=faida->GetMaxHitPositionX();
            int betaminy=faida->GetMinHitPositionY();
            int betamaxy=faida->GetMaxHitPositionY();
            fTsNow = faida->GetTimestamp();
            //cout<<"ftsnow = "<<fTsNow<<endl;
            
            //! correlation start
            Long64_t ts1 = (Long64_t)ts - 10e9;
            Long64_t ts2 = (Long64_t)ts + 10e9;
            Long64_t corrts = 0;
            Long64_t check_time = 0;
            
            fimplantMap_it = fimplantMap.lower_bound(ts1);

            while(fimplantMap_it!=fimplantMap.end()&&fimplantMap_it->first<ts2){
                corrts = (Long64_t) fimplantMap_it->first;
                AIDASimpleStruct* imp = (AIDASimpleStruct*) fimplantMap_it->second.second;
                
                short impz= imp->GetHitPositionZ()+imp->GetDZ(); // corrected dZ
                
                int impminx=imp->GetMinHitPositionX();
                int impmaxx=imp->GetMaxHitPositionX();
                int impminy=imp->GetMinHitPositionY();
                int impmaxy=imp->GetMaxHitPositionY();
                
                if (corrts!=check_time&&betaz==impz){// avoid multiple filling corrts!=check_time
                    //! implant-beta spatial correlation
                    //! new using overlaping area
                    if (!( ((impminx<=betamaxx&&impminx>=betaminx)||(impmaxx<=betamaxx&&impmaxx>=betaminx))&&
                        ((impminy<=betamaxy&&impminy>=betaminy)||(impmaxy<=betamaxy&&impmaxy>=betaminy)) )) {
                        fimplantMap_it++;
                        continue;
                    }
                    
                    //! Fill correlation
                    //cout<<"merge "<< corrts << "\t" << faida->GetTimestamp()<<endl;
                    
                    
                    check_time=corrts;
                }
                
                fimplantMap_it++;
            }
            //! correlation end
            
        }
    }

    cout<<"finished"<<endl;
}

mergertest((char*)"aidarootfiles/aida_R6_427to446.root",(char*)"belenrootfiles/belen42.root",(char*)"bigripsrootfiles/run3039_CORRECTED.root");
