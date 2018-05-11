#include "AIDA.h"

bool AIDA::BetaGetPos(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[])
{
    int maxmult=1;
    Double_t dssdH_E_X[NumDSSD][NumStrX][maxmult];
    Double_t dssdH_E_Y[NumDSSD][NumStrY][maxmult];
    Long64_t dssdH_T_X[NumDSSD][NumStrX][maxmult];
    Long64_t dssdH_T_Y[NumDSSD][NumStrY][maxmult];

    Long64_t dssdF_T_X[NumDSSD][NumStrX][maxmult];
    Long64_t dssdF_T_Y[NumDSSD][NumStrY][maxmult];

    int mult_strip_x[NumDSSD][NumStrX];
    int mult_strip_y[NumDSSD][NumStrX];
    for (int i=0;i<NumDSSD;i++){
        for (int j=0;j<NumStrX;j++){
            mult_strip_x[i][j]=0;
            mult_strip_y[i][j]=0;
            for (int k=0;k<maxmult;k++){
                dssdH_E_X[i][j][k]=0;
                dssdH_E_Y[i][j][k]=0;
                dssdH_T_X[i][j][k]=0;
                dssdF_T_Y[i][j][k]=0;
                dssdF_T_X[i][j][k]=0;
                dssdF_T_Y[i][j][k]=0;
            }
        }
    }
    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        int xy = hit->GetXY();
        double energy = hit->GetEnergy();
        unsigned long long time = hit->GetTimestamp();
        unsigned long long fastime = hit->GetFastTimestamp();
        if (xy<128){
            if (mult_strip_x[z][xy]<maxmult) {
                dssdH_E_X[z][xy][mult_strip_x[z][xy]] = energy;
                dssdH_T_X[z][xy][mult_strip_x[z][xy]] = time;
                dssdF_T_X[z][xy][mult_strip_x[z][xy]] = fastime;
            }
            mult_strip_x[z][xy]++;
        }else{
            if (mult_strip_y[z][xy-128]<maxmult) {
                dssdH_E_Y[z][xy-128][mult_strip_y[z][xy-128]] = energy;
                dssdH_T_Y[z][xy-128][mult_strip_y[z][xy-128]] = time;
                dssdF_T_Y[z][xy-128][mult_strip_y[z][xy-128]] = fastime;
            }
            mult_strip_y[z][xy-128]++;
        }
    }

    //! Old code

    Int_t maxStrInCluster = 1000;
    Int_t maxNposBeta = 1000;
    Int_t maxNCluster = 1000;

    for(Int_t z=0; z<NumDSSD; z++){
        Int_t mult_z = 0;
        vector<pair<Double_t,pair<pair<Double_t,pair<Double_t,Int_t> >,pair<Double_t,pair<Double_t,Int_t> > >  > > beta_dssd_pre;
        //vector<pair<E-corr,pair<pair<PosX,pair<EX,NClusterX> >,pair<PosY,pair<EY,NClusterY> > >  > >
        vector<pair<Double_t,pair<pair<Double_t,pair<Double_t,Int_t> >,pair<Double_t,pair<Double_t,Int_t> > >  > >::iterator ibeta_dssd_pre;

        Double_t posX=-1;
        Double_t posY=-1;

        Int_t nClusterX=0;
        Int_t nClusterY=0;
        Int_t nStripsX=0;
        Double_t E_X=0;
        Double_t E_X_ch=0;

        for(Int_t x=0; x<NumStrX; x++){
            //cout << "x-strip found" << endl;
            if (dssdH_E_X[z][x][0]>0) {
                E_X+=dssdH_E_X[z][x][0];
                E_X_ch+=x*dssdH_E_X[z][x][0];
                nStripsX++;
            }
            if ((dssdH_E_X[z][x][0]<=0||x==NumStrX-1)&&nStripsX>0) {//end of X cluster
                if (nStripsX<=maxStrInCluster){ //if cluster of less than 3 strips
                    posX=E_X_ch/E_X;
                    //another loop for Y strips
                    nClusterY=0;
                    Int_t nStripsY=0;
                    Double_t E_Y=0;
                    Double_t E_Y_ch=0;
                    for(Int_t y=0; y<NumStrY; y++){
                        //cout << "x-strip found" << endl;
                        if (dssdH_E_Y[z][y][0]>0) {
                            E_Y+=(Double_t) dssdH_E_Y[z][y][0];
                            E_Y_ch+=(Double_t) y*dssdH_E_Y[z][y][0];
                            nStripsY++;
                        }
                        if ((dssdH_E_Y[z][y][0]<=0||y==NumStrY-1)&&nStripsY>0) {//end of X cluster
                            if (nStripsY<=maxStrInCluster){ //if cluster of less than 3 strips
                                posY=E_Y_ch/E_Y;
                                //cout<<"x"<<posX <<"y"<<posY<<"e" <<(E_Y/E_X-1)*(E_Y/E_X-1)<<endl;
                                beta_dssd_pre.push_back(make_pair((E_Y/E_X-1)*(E_Y/E_X-1),make_pair(make_pair(posX,make_pair(E_X,nStripsX)),make_pair(posY,make_pair(E_Y,nStripsY)))));
                            }
                            E_Y=0;
                            E_Y_ch=0;
                            nStripsY=0;
                            nClusterY++;
                        }
                    }
                    //ok!
                }
                E_X=0;
                E_X_ch=0;
                nStripsX=0;
                nClusterX++;
            }
        }

        //Additional step to filter all possible combination:
        if (nClusterX>0&&nClusterX<maxNCluster&&nClusterY>0&&nClusterY<maxNCluster){
            //Record number of cluster here
            vector <Double_t> xindex;
            vector <Double_t> yindex;
            sort(beta_dssd_pre.begin(),beta_dssd_pre.end());
            for (ibeta_dssd_pre=beta_dssd_pre.begin();ibeta_dssd_pre<beta_dssd_pre.end();++ibeta_dssd_pre){
                if (((corr_cut<=0)&&(find(xindex.begin(),xindex.end(),ibeta_dssd_pre->second.first.first)==xindex.end())&&(find(yindex.begin(),yindex.end(),ibeta_dssd_pre->second.second.first)==yindex.end()))||
                        ((corr_cut>0)&&(ibeta_dssd_pre->first<corr_cut*corr_cut)&&(find(xindex.begin(),xindex.end(),ibeta_dssd_pre->second.first.first)==xindex.end())&&(find(yindex.begin(),yindex.end(),ibeta_dssd_pre->second.second.first)==yindex.end())))
                {
                    //beta_dssd.push_back(make_pair(make_pair(ibeta_dssd_pre->second.first.first,ibeta_dssd_pre->second.first.second),make_pair(ibeta_dssd_pre->second.second.first,ibeta_dssd_pre->second.second.second)));
                    if (mult_z<maxNposBeta && ibeta_dssd_pre->second.first.second.first>sumexcut[z] && ibeta_dssd_pre->second.second.second.first>sumeycut[z]){
                        AIDACluster *cluster=new AIDACluster;
                        cluster->SetHitPosition(ibeta_dssd_pre->second.first.first,ibeta_dssd_pre->second.second.first,z);
                        cluster->SetXEnergy(ibeta_dssd_pre->second.first.second.first);
                        cluster->SetYEnergy(ibeta_dssd_pre->second.second.second.first);
                        cluster->SetXMult(ibeta_dssd_pre->second.first.second.second);
                        cluster->SetYMult(ibeta_dssd_pre->second.second.second.second);
                        //if (ibeta_dssd_pre->second.first.second.second!=ibeta_dssd_pre->second.second.second.second) cout<<"eurica!"<<cluster->GetXMultiplicity()<<"-"<<cluster->GetYMultiplicity()<<endl;

                        //Timing
                        int hitx = (int) round(ibeta_dssd_pre->second.first.first);
                        int hity = (int) round(ibeta_dssd_pre->second.second.first);
                        if (dssdH_T_X[z][hitx][0]>0&&dssdH_T_Y[z][hity][0]>0){

                            //! Take the ealiest time stamp
                            if (dssdH_T_X[z][hitx][0]>dssdH_T_Y[z][hity][0]) {
                                cluster->SetTimestamp(dssdH_T_Y[z][hity][0]);
                                cluster->SetFastTimestamp(dssdF_T_Y[z][hity][0]);
                            }else {
                                cluster->SetTimestamp(dssdH_T_X[z][hitx][0]);
                                cluster->SetFastTimestamp(dssdF_T_X[z][hitx][0]);
                            }
                            /*
                            //! if there is available fast time stamp
                            if (dssdF_T_X[z][hitx][0]==0||dssdF_T_Y[z][hity][0]==0){
                                cluster->SetFastTimestamp(dssdF_T_X[z][hitx][0] + dssdF_T_Y[z][hity][0]);
                            }
                            */

                        }else{
                            cout<<__PRETTY_FUNCTION__<<" something wrong!"<<endl;
                        }
                        this->AddCluster(cluster);
                    }
                    xindex.push_back(ibeta_dssd_pre->second.first.first);
                    yindex.push_back(ibeta_dssd_pre->second.second.first);
                    mult_z++;
                }
            }
        }//end of additional step...

    }//end of loop on all dssd
    if(this->GetNClusters()>0) return true;
    return false;
}

bool AIDA::BetaGetPosNew(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[])
{
    //!CAUTION: condition if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
    //! is apply so we should assume "good" calibrationn


    //! A  tiny event builder by dimension!
    //! So far the method can only handle the first hit in each strip!
    Int_t maxStrInCluster = 1000;
    Int_t maxNposBeta = 1000;
    Int_t maxNCluster = 1000;

    Int_t nClusterX[NumDSSD];
    Int_t nClusterY[NumDSSD];
    memset(nClusterX,0,sizeof(nClusterX));
    memset(nClusterY,0,sizeof(nClusterY));

    //! Prepare maps of spatial distribution of hits
    //! NOTICE: remember to dallocate AIDAHit* !
    //! NOTICE: we use map: just select the first hit in each channel
    vector< map<short, AIDAHit* > > xhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator xhitmaps_it;
    vector< map<short, AIDAHit* > > yhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator yhitmaps_it;

    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        short xy = hit->GetXY();
        if (xy<128){
            xhitmaps[z].insert(std::make_pair(xy, hit));
            //if (z==0&&hit->GetEnergy()>0) cout<<xy<<"b-"<<hit->GetEnergy()<<endl;
        }else{
            yhitmaps[z].insert(std::make_pair(xy-128, hit));
        }
    }

    //cluster map sort by EX EY
    vector< map<Double_t,AIDACluster*> > cluster_map(NumDSSD);
    //multimap<E-corr, pair<Xhit,Yhit> >
    map<Double_t,AIDACluster*>::iterator cluster_map_it;

    Double_t posX=-11.;
    Double_t posY=-11.;
    //if (thereis) cout<<"----"<<endl;

    //! make them clusters!
    for (int z = 0;z < NumDSSD;z++){
        short x_prev = -11;
        short nStripInClusterX = 0;
        double E_X = 0;
        double E_X_ch = 0;
        //! loop through the x axis

        for(xhitmaps_it = xhitmaps[z].begin(); xhitmaps_it != xhitmaps[z].end(); xhitmaps_it++){
            short x = xhitmaps_it->first;
            double xen = xhitmaps_it->second->GetEnergy();
            //if (this->GetHit(0)->GetTimestamp()==23112644336) cout<<z<<"eee-"<<x<<"-"<<xen<<endl;
            //if (z==0) cout<<x<<"cc-"<<xen<<endl;

            //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"*"<<endl;

            //!out of cluster
            if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;
                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;
                    nClusterY[z] = 0;//ok...
                    //! loop through the  y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){

                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                //cout<<"x"<<posX <<"y"<<posY<<"e" <<(E_Y/E_X-1)*(E_Y/E_X-1)<<endl;

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        y_prev = y;
                    }//y
                     nClusterX[z]++;
                }//copy to the last cluster handling in X

                //****************************************//

                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;
            }
            //!still in a same cluster
            if ( xen>0){
                //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"xprev"<<x_prev<<endl;
                E_X += xen;
                E_X_ch += (double) x * xen;
                nStripInClusterX++;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<x<<"dd-"<<xen<<"-"<<nStripInClusterX<<endl;
            }
            //! handle last cluster!
            if (xhitmaps_it == --xhitmaps[z].end() && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;

                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;
                    nClusterY[z] = 0;//ok...
                    //! loop through the y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }
                        y_prev = y;
                    }//y

                    nClusterX[z]++;
                }//copy to the last cluster handling in X
                //****************************************//
                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;

            }
            x_prev = x;
        }//x
    }//z

    Int_t maxZ=-1;

    for (int z = 0;z < NumDSSD;z++){
        Int_t mult_z = 0;
        //Record number of cluster here
        vector <Double_t> xindex;
        vector <Double_t> yindex;
        if (nClusterX[z]>0&&nClusterX[z]<maxNCluster&&nClusterY[z]>0&&nClusterY[z]<maxNCluster){
            for(cluster_map_it = cluster_map[z].begin(); cluster_map_it != cluster_map[z].end(); cluster_map_it++){
                AIDACluster* cluster = cluster_map_it->second;
                //! to fix bug on .9999999 number
                double xx=cluster->GetHitPositionX();
                double yy=cluster->GetHitPositionY();
                if(cluster->GetXMultiplicity()==1) xx=round(xx);
                if(cluster->GetYMultiplicity()==1) yy=round(yy);
                cluster->SetHitPosition(xx,yy,cluster->GetHitPositionZ());

                if (((corr_cut<=0)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end()))||
                        ((corr_cut>0)&&(cluster_map_it->first<corr_cut)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end())))
                {
                    //cout<<cluster_map_it->first<<"-"<<cluster->GetXEnergy()<<"-"<<cluster->GetYEnergy()<<endl;
                    if (mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        //! Deallocating memory
                        delete cluster;
                    }
                    xindex.push_back(cluster->GetHitPositionX());
                    yindex.push_back(cluster->GetHitPositionY());
                    mult_z++;
                }else{
                    //! Deallocating memory
                    delete cluster;
                }
            }//loop all cluster map
        }//if condition
    }//all z


    if(this->GetNClusters()>0) {
        this->SetMaxZ(maxZ);
        return true;
    }
    return false;

}


bool AIDA::BetaGetPosAllNew(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[])
{
    //!CAUTION: condition if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
    //! is apply so we should assume "good" calibrationn


    //! A  tiny event builder by dimension!
    //! So far the method can only handle the first hit in each strip!
    Int_t maxStrInCluster = 1000;
    Int_t maxNposBeta = 1000;
    Int_t maxNCluster = 1000;

    Int_t nClusterX[NumDSSD];
    Int_t nClusterY[NumDSSD];
    memset(nClusterX,0,sizeof(nClusterX));
    memset(nClusterY,0,sizeof(nClusterY));

    //! Prepare maps of spatial distribution of hits
    //! NOTICE: remember to dallocate AIDAHit* !
    //! NOTICE: we use map: just select the first hit in each channel
    vector< map<short, AIDAHit* > > xhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator xhitmaps_it;
    vector< map<short, AIDAHit* > > yhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator yhitmaps_it;

    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        short xy = hit->GetXY();
        if (xy<128){
            xhitmaps[z].insert(std::make_pair(xy, hit));
            //if (z==0&&hit->GetEnergy()>0) cout<<xy<<"b-"<<hit->GetEnergy()<<endl;
        }else{
            yhitmaps[z].insert(std::make_pair(xy-128, hit));
        }
    }

    //cluster map sort by EX EY
    vector< map<Double_t,AIDACluster*> > cluster_map(NumDSSD);
    //multimap<E-corr, pair<Xhit,Yhit> >
    map<Double_t,AIDACluster*>::iterator cluster_map_it;

    Double_t posX=-11.;
    Double_t posY=-11.;
    //if (thereis) cout<<"----"<<endl;

    //! make them clusters!
    for (int z = 0;z < NumDSSD;z++){
        short x_prev = -11;
        short nStripInClusterX = 0;
        double E_X = 0;
        double E_X_ch = 0;
        //! loop through the x axis

        for(xhitmaps_it = xhitmaps[z].begin(); xhitmaps_it != xhitmaps[z].end(); xhitmaps_it++){
            short x = xhitmaps_it->first;
            double xen = xhitmaps_it->second->GetEnergy();
            //if (this->GetHit(0)->GetTimestamp()==23112644336) cout<<z<<"eee-"<<x<<"-"<<xen<<endl;
            //if (z==0) cout<<x<<"cc-"<<xen<<endl;

            //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"*"<<endl;

            //!out of cluster
            if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;
                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;
                    nClusterY[z] = 0;//ok...
                    //! loop through the  y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){

                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                //cout<<"x"<<posX <<"y"<<posY<<"e" <<(E_Y/E_X-1)*(E_Y/E_X-1)<<endl;

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        y_prev = y;
                    }//y
                     nClusterX[z]++;
                }//copy to the last cluster handling in X

                //****************************************//

                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;
            }
            //!still in a same cluster
            if ( xen>0){
                //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"xprev"<<x_prev<<endl;
                E_X += xen;
                E_X_ch += (double) x * xen;
                nStripInClusterX++;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<x<<"dd-"<<xen<<"-"<<nStripInClusterX<<endl;
            }
            //! handle last cluster!
            if (xhitmaps_it == --xhitmaps[z].end() && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;

                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;
                    nClusterY[z] = 0;//ok...
                    //! loop through the y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;                                
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());


                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }
                        y_prev = y;
                    }//y

                    nClusterX[z]++;
                }//copy to the last cluster handling in X
                //****************************************//
                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;

            }
            x_prev = x;
        }//x
    }//z

    Int_t maxZ=-1;

    for (int z = 0;z < NumDSSD;z++){
        Int_t mult_z = 0;
        //Record number of cluster here
        vector <Double_t> xindex;
        vector <Double_t> yindex;
        if (nClusterX[z]>0&&nClusterX[z]<maxNCluster&&nClusterY[z]>0&&nClusterY[z]<maxNCluster){
            for(cluster_map_it = cluster_map[z].begin(); cluster_map_it != cluster_map[z].end(); cluster_map_it++){
                AIDACluster* cluster = cluster_map_it->second;

                //! to fix bug on .9999999 number
                double xx=cluster->GetHitPositionX();
                double yy=cluster->GetHitPositionY();
                if(cluster->GetXMultiplicity()==1) xx=round(xx);
                if(cluster->GetYMultiplicity()==1) yy=round(yy);
                cluster->SetHitPosition(xx,yy,cluster->GetHitPositionZ());



                /*
                if (((corr_cut<=0)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end()))||
                        ((corr_cut>0)&&(cluster_map_it->first<corr_cut)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end())))
                {
                    //cout<<cluster_map_it->first<<"-"<<cluster->GetXEnergy()<<"-"<<cluster->GetYEnergy()<<endl;
                    if (mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        //! Deallocating memory
                        delete cluster;
                    }
                    xindex.push_back(cluster->GetHitPositionX());
                    yindex.push_back(cluster->GetHitPositionY());
                    mult_z++;

                }else{
                    //! Deallocating memory
                    delete cluster;
                }
                */
                if ((corr_cut<=0)||(corr_cut>0)&&(cluster_map_it->first<corr_cut)&&mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                    this->AddCluster(cluster);
                    maxZ = z;
                }
                else delete cluster;
            }//loop all cluster map
        }//if condition
    }//all z


    if(this->GetNClusters()>0) {
        this->SetMaxZ(maxZ);
        return true;
    }
    return false;
}

bool AIDA::IonGetPos()
{
    int maxmult=1;
    Double_t dssdL_E_X[NumDSSD][NumStrX][maxmult];
    Double_t dssdL_E_Y[NumDSSD][NumStrY][maxmult];

    unsigned long long dssdL_T_X[NumDSSD][NumStrX][maxmult];
    unsigned long long dssdL_T_Y[NumDSSD][NumStrX][maxmult];

    int mult_strip_x[NumDSSD][NumStrX];
    int mult_strip_y[NumDSSD][NumStrX];
    for (int i=0;i<NumDSSD;i++){
        for (int j=0;j<NumStrX;j++){
            mult_strip_x[i][j]=0;
            mult_strip_y[i][j]=0;
            for (int k=0;k<maxmult;k++){
                dssdL_E_X[i][j][k]=0;
                dssdL_E_Y[i][j][k]=0;
                dssdL_T_X[i][j][k]=0;
                dssdL_T_Y[i][j][k]=0;
            }
        }
    }
    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);

        int z = hit->GetZ();
        int xy = hit->GetXY();
        double energy = hit->GetEnergy();
        unsigned long long timestamp = hit->GetTimestamp();
        if (xy<128){
            if (mult_strip_x[z][xy]<maxmult) {
                dssdL_E_X[z][xy][mult_strip_x[z][xy]] = energy;
                dssdL_T_X[z][xy][mult_strip_x[z][xy]]=timestamp;
            }
            mult_strip_x[z][xy]++;
        }else{
            if (mult_strip_y[z][xy-128]<maxmult) {
                dssdL_E_Y[z][xy-128][mult_strip_y[z][xy-128]] = energy;
                dssdL_T_Y[z][xy-128][mult_strip_y[z][xy-128]]=timestamp;
            }
            mult_strip_y[z][xy-128]++;
        }
    }
    //!Old code
    //get pos by weighting method
    Double_t Ion_posX;
    Double_t Ion_posY;
    Bool_t is_X_hit=false;
    Bool_t is_Y_hit=false;
    Int_t my_ion_mult_x,my_ion_mult_y;
    Int_t maxZ=-1;

    for(Int_t z=0; z<NumDSSD; z++){
        is_X_hit=false;
        is_Y_hit=false;
        my_ion_mult_x=0;
        my_ion_mult_y=0;
        Double_t sum_ex=0;
        Double_t sum_weightx=0;
        Double_t sum_ey=0;
        Double_t sum_weighty=0;
        //! newly added
        Double_t maxEX = 0;
        Int_t maxX = 0;
        Double_t maxEY = 0;
        Int_t maxY = 0;
        for(Int_t x=0; x<NumStrX; x++){
            if(dssdL_E_X[z][x][0]>0){
                is_X_hit=true;
                sum_ex+=dssdL_E_X[z][x][0];
                sum_weightx+=x*dssdL_E_X[z][x][0];
                my_ion_mult_x++;
                if (dssdL_E_X[z][x][0]>maxEX) {
                    maxX = x;
                    maxEX = dssdL_E_X[z][x][0];
                }
            }
        }
        if (is_X_hit){
            for(Int_t y=0; y<NumStrY; y++){
                if(dssdL_E_Y[z][y][0]>0){
                   is_Y_hit=true;
                   sum_ey+=dssdL_E_Y[z][y][0];
                   sum_weighty+=y*dssdL_E_Y[z][y][0];
                   my_ion_mult_y++;
                   if (dssdL_E_Y[z][y][0]>maxEY) {
                       maxY = y;
                       maxEY = dssdL_E_Y[z][y][0];
                   }
               }
            }
            if(is_Y_hit){
                Ion_posX=sum_weightx/sum_ex;
                Ion_posY=sum_weighty/sum_ey;
                //! only one cluster
                AIDACluster *cluster=new AIDACluster;
                cluster->SetXEnergy(sum_ex);
                cluster->SetYEnergy(sum_ey);
                cluster->SetHitPosition(Ion_posX,Ion_posY,z);
                cluster->SetXMult(my_ion_mult_x);
                cluster->SetYMult(my_ion_mult_y);
                //! take ealiest time stamp

                //cluster->SetTimestamp(fhits.at(0)->GetTimestamp());
                //! newly added from the Jose email
                if (maxEX>maxEY) cluster->SetTimestamp(dssdL_T_X[z][maxX][0]);
                else cluster->SetTimestamp(dssdL_T_Y[z][maxY][0]);

                this->AddCluster(cluster);
                maxZ=z;
            }
        }
    }
    this->SetMaxZ(maxZ);
    if (maxZ!=-1) return true;
    return false;
}


bool AIDA::IonGetPosNew()//Clustering algorithm applied for Ion events
{
    Double_t corr_cut=-1;
    int maxmult=1;
    Double_t dssdH_E_X[NumDSSD][NumStrX][maxmult];
    Double_t dssdH_E_Y[NumDSSD][NumStrY][maxmult];
    Long64_t dssdH_T_X[NumDSSD][NumStrX][maxmult];
    Long64_t dssdH_T_Y[NumDSSD][NumStrY][maxmult];

    Long64_t dssdF_T_X[NumDSSD][NumStrX][maxmult];
    Long64_t dssdF_T_Y[NumDSSD][NumStrY][maxmult];

    int mult_strip_x[NumDSSD][NumStrX];
    int mult_strip_y[NumDSSD][NumStrX];
    for (int i=0;i<NumDSSD;i++){
        for (int j=0;j<NumStrX;j++){
            mult_strip_x[i][j]=0;
            mult_strip_y[i][j]=0;
            for (int k=0;k<maxmult;k++){
                dssdH_E_X[i][j][k]=0;
                dssdH_E_Y[i][j][k]=0;
                dssdH_T_X[i][j][k]=0;
                dssdF_T_Y[i][j][k]=0;
                dssdF_T_X[i][j][k]=0;
                dssdF_T_Y[i][j][k]=0;
            }
        }
    }
    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        int xy = hit->GetXY();
        double energy = hit->GetEnergy();
        unsigned long long time = hit->GetTimestamp();
        unsigned long long fastime = hit->GetFastTimestamp();
        if (xy<128){
            if (mult_strip_x[z][xy]<maxmult) {
                dssdH_E_X[z][xy][mult_strip_x[z][xy]] = energy;
                dssdH_T_X[z][xy][mult_strip_x[z][xy]] = time;
                dssdF_T_X[z][xy][mult_strip_x[z][xy]] = fastime;
            }
            mult_strip_x[z][xy]++;
        }else{
            if (mult_strip_y[z][xy-128]<maxmult) {
                dssdH_E_Y[z][xy-128][mult_strip_y[z][xy-128]] = energy;
                dssdH_T_Y[z][xy-128][mult_strip_y[z][xy-128]] = time;
                dssdF_T_Y[z][xy-128][mult_strip_y[z][xy-128]] = fastime;
            }
            mult_strip_y[z][xy-128]++;
        }
    }

    //! Old code

    Int_t maxStrInCluster = 1000;
    Int_t maxNposBeta = 1000;
    Int_t maxNCluster = 1000;
    Int_t maxZ=-1;

    for(Int_t z=0; z<NumDSSD; z++){
        Int_t mult_z = 0;
        vector<pair<Double_t,pair<pair<Double_t,pair<Double_t,Int_t> >,pair<Double_t,pair<Double_t,Int_t> > >  > > beta_dssd_pre;
        //vector<pair<E-corr,pair<pair<PosX,pair<EX,NClusterX> >,pair<PosY,pair<EY,NClusterY> > >  > >
        vector<pair<Double_t,pair<pair<Double_t,pair<Double_t,Int_t> >,pair<Double_t,pair<Double_t,Int_t> > >  > >::iterator ibeta_dssd_pre;

        Double_t posX=-1;
        Double_t posY=-1;

        Int_t nClusterX=0;
        Int_t nClusterY=0;
        Int_t nStripsX=0;
        Double_t E_X=0;
        Double_t E_X_ch=0;

        for(Int_t x=0; x<NumStrX; x++){
            //cout << "x-strip found" << endl;
            if (dssdH_E_X[z][x][0]>0) {
                E_X+=dssdH_E_X[z][x][0];
                E_X_ch+=(Double_t)x*(Double_t)dssdH_E_X[z][x][0];
                nStripsX++;
            }
            if ((dssdH_E_X[z][x][0]<=0||x==NumStrX-1)&&nStripsX>0) {//end of X cluster
                if (nStripsX<=maxStrInCluster){ //if cluster of less than 3 strips
                    posX=E_X_ch/E_X;
                    //another loop for Y strips
                    nClusterY=0;
                    Int_t nStripsY=0;
                    Double_t E_Y=0;
                    Double_t E_Y_ch=0;
                    for(Int_t y=0; y<NumStrY; y++){
                        //cout << "x-strip found" << endl;
                        if (dssdH_E_Y[z][y][0]>0) {
                            E_Y+=(Double_t)dssdH_E_Y[z][y][0];
                            E_Y_ch+=(Double_t)y*(Double_t)dssdH_E_Y[z][y][0];
                            nStripsY++;
                        }
                        if ((dssdH_E_Y[z][y][0]<=0||y==NumStrY-1)&&nStripsY>0) {//end of X cluster
                            if (nStripsY<=maxStrInCluster){ //if cluster of less than 3 strips
                                posY=E_Y_ch/E_Y;
                                //cout<<"x"<<posX <<"y"<<posY<<"e" <<(E_Y/E_X-1)*(E_Y/E_X-1)<<endl;

                                //if (nStripsX==1) posX=round(posX);
                                //if (nStripsY==1) posY=round(posY);
                                if (E_X>E_Y)
                                beta_dssd_pre.push_back(make_pair(1-E_Y/E_X,make_pair(make_pair(posX,make_pair(E_X,nStripsX)),make_pair(posY,make_pair(E_Y,nStripsY)))));
                                else
                                beta_dssd_pre.push_back(make_pair(1-E_X/E_Y,make_pair(make_pair(posX,make_pair(E_X,nStripsX)),make_pair(posY,make_pair(E_Y,nStripsY)))));
                            }
                            E_Y=0;
                            E_Y_ch=0;
                            nStripsY=0;
                            nClusterY++;
                        }
                    }
                    //ok!
                }
                E_X=0;
                E_X_ch=0;
                nStripsX=0;
                nClusterX++;
            }
        }

        //Additional step to filter all possible combination:
        if (nClusterX>0&&nClusterX<maxNCluster&&nClusterY>0&&nClusterY<maxNCluster){
            //Record number of cluster here
            vector <Double_t> xindex;
            vector <Double_t> yindex;
            sort(beta_dssd_pre.begin(),beta_dssd_pre.end());
            for (ibeta_dssd_pre=beta_dssd_pre.begin();ibeta_dssd_pre<beta_dssd_pre.end();++ibeta_dssd_pre){
                if (((corr_cut<=0)&&(find(xindex.begin(),xindex.end(),ibeta_dssd_pre->second.first.first)==xindex.end())&&(find(yindex.begin(),yindex.end(),ibeta_dssd_pre->second.second.first)==yindex.end()))||
                        ((corr_cut>0)&&(ibeta_dssd_pre->first<corr_cut)&&(find(xindex.begin(),xindex.end(),ibeta_dssd_pre->second.first.first)==xindex.end())&&(find(yindex.begin(),yindex.end(),ibeta_dssd_pre->second.second.first)==yindex.end())))
                {
                    //beta_dssd.push_back(make_pair(make_pair(ibeta_dssd_pre->second.first.first,ibeta_dssd_pre->second.first.second),make_pair(ibeta_dssd_pre->second.second.first,ibeta_dssd_pre->second.second.second)));
                    if (mult_z<maxNposBeta){
                        AIDACluster *cluster=new AIDACluster;
                        double xx=ibeta_dssd_pre->second.first.first;
                        double yy=ibeta_dssd_pre->second.second.first;
                        if(ibeta_dssd_pre->second.first.second.second==1) xx=round(xx);
                        if(ibeta_dssd_pre->second.second.second.second==1) yy=round(yy);
                        cluster->SetHitPosition(xx,yy,z);
                        //cluster->SetHitPosition(ibeta_dssd_pre->second.first.first,ibeta_dssd_pre->second.second.first,z);
                        cluster->SetXEnergy(ibeta_dssd_pre->second.first.second.first);
                        cluster->SetYEnergy(ibeta_dssd_pre->second.second.second.first);
                        cluster->SetXMult(ibeta_dssd_pre->second.first.second.second);
                        cluster->SetYMult(ibeta_dssd_pre->second.second.second.second);
                        //if (ibeta_dssd_pre->second.first.second.second!=ibeta_dssd_pre->second.second.second.second) cout<<"eurica!"<<cluster->GetXMultiplicity()<<"-"<<cluster->GetYMultiplicity()<<endl;
                        //Timing
                        int hitx = (int) round(ibeta_dssd_pre->second.first.first);
                        int hity = (int) round(ibeta_dssd_pre->second.second.first);
                        if (dssdH_T_X[z][hitx][0]>0&&dssdH_T_Y[z][hity][0]>0){

                            //! Take the ealiest time stamp
                            if (dssdH_T_X[z][hitx][0]>dssdH_T_Y[z][hity][0]) {
                                cluster->SetTimestamp(dssdH_T_Y[z][hity][0]);
                                cluster->SetFastTimestamp(dssdF_T_Y[z][hity][0]);
                            }else{
                                cluster->SetTimestamp(dssdH_T_X[z][hitx][0]);
                                cluster->SetFastTimestamp(dssdF_T_X[z][hitx][0]);
                            }
                            //! if there is available fast time stamp
                            if (dssdF_T_X[z][hitx][0]==0||dssdF_T_Y[z][hity][0]==0){
                                cluster->SetFastTimestamp(dssdF_T_X[z][hitx][0] + dssdF_T_Y[z][hity][0]);
                            }

                        }else{
                            cout<<__PRETTY_FUNCTION__<<" something wrong!"<<endl;
                        }
                        this->AddCluster(cluster);
                        maxZ=z;
                    }
                    xindex.push_back(ibeta_dssd_pre->second.first.first);
                    yindex.push_back(ibeta_dssd_pre->second.second.first);
                    mult_z++;
                }
            }
        }//end of additional step...

    }//end of loop on all dssd
    this->SetMaxZ(maxZ);
    if(this->GetNClusters()>0) return true;
    return false;
}


bool AIDA::IonGetPosAllNew()//Clustering algorithm applied for Ion events
{
    int maxmult=1;
    Double_t dssdH_E_X[NumDSSD][NumStrX][maxmult];
    Double_t dssdH_E_Y[NumDSSD][NumStrY][maxmult];
    Long64_t dssdH_T_X[NumDSSD][NumStrX][maxmult];
    Long64_t dssdH_T_Y[NumDSSD][NumStrY][maxmult];

    Long64_t dssdF_T_X[NumDSSD][NumStrX][maxmult];
    Long64_t dssdF_T_Y[NumDSSD][NumStrY][maxmult];

    int mult_strip_x[NumDSSD][NumStrX];
    int mult_strip_y[NumDSSD][NumStrX];
    for (int i=0;i<NumDSSD;i++){
        for (int j=0;j<NumStrX;j++){
            mult_strip_x[i][j]=0;
            mult_strip_y[i][j]=0;
            for (int k=0;k<maxmult;k++){
                dssdH_E_X[i][j][k]=0;
                dssdH_E_Y[i][j][k]=0;
                dssdH_T_X[i][j][k]=0;
                dssdF_T_Y[i][j][k]=0;
                dssdF_T_X[i][j][k]=0;
                dssdF_T_Y[i][j][k]=0;
            }
        }
    }
    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        int xy = hit->GetXY();
        double energy = hit->GetEnergy();
        unsigned long long time = hit->GetTimestamp();
        unsigned long long fastime = hit->GetFastTimestamp();
        if (xy<128){
            if (mult_strip_x[z][xy]<maxmult) {
                dssdH_E_X[z][xy][mult_strip_x[z][xy]] = energy;
                dssdH_T_X[z][xy][mult_strip_x[z][xy]] = time;
                dssdF_T_X[z][xy][mult_strip_x[z][xy]] = fastime;
            }
            mult_strip_x[z][xy]++;
        }else{
            if (mult_strip_y[z][xy-128]<maxmult) {
                dssdH_E_Y[z][xy-128][mult_strip_y[z][xy-128]] = energy;
                dssdH_T_Y[z][xy-128][mult_strip_y[z][xy-128]] = time;
                dssdF_T_Y[z][xy-128][mult_strip_y[z][xy-128]] = fastime;
            }
            mult_strip_y[z][xy-128]++;
        }
    }

    //! Old code

    Int_t maxStrInCluster = 1000;
    Int_t maxNposBeta = 1000;
    Int_t maxNCluster = 1000;
    Int_t maxZ=-1;

    for(Int_t z=0; z<NumDSSD; z++){
        Int_t mult_z = 0;

        vector<pair<Double_t,pair<pair<Double_t,pair<Double_t,Int_t> >,pair<Double_t,pair<Double_t,Int_t> > >  > > beta_dssd_pre;
        vector<pair<Double_t,pair<pair<Double_t,pair<Double_t,Int_t> >,pair<Double_t,pair<Double_t,Int_t> > >  > >::iterator ibeta_dssd_pre;

        Double_t posX=-1;
        Double_t posY=-1;

        Int_t nClusterX=0;
        Int_t nClusterY=0;
        Int_t nStripsX=0;
        Double_t E_X=0;
        Double_t E_X_ch=0;

        for(Int_t x=0; x<NumStrX; x++){
            //cout << "x-strip found" << endl;
            if (dssdH_E_X[z][x][0]>0) {
                E_X+=dssdH_E_X[z][x][0];
                E_X_ch+=x*dssdH_E_X[z][x][0];
                nStripsX++;
            }
            if ((dssdH_E_X[z][x][0]<=0||x==NumStrX-1)&&nStripsX>0) {//end of X cluster
                if (nStripsX<=maxStrInCluster){ //if cluster of less than 3 strips
                    posX=E_X_ch/E_X;
                    //another loop for Y strips
                    nClusterY=0;
                    Int_t nStripsY=0;
                    Double_t E_Y=0;
                    Double_t E_Y_ch=0;
                    for(Int_t y=0; y<NumStrY; y++){
                        //cout << "x-strip found" << endl;
                        if (dssdH_E_Y[z][y][0]>0) {
                            E_Y+=(Double_t) dssdH_E_Y[z][y][0];
                            E_Y_ch+=(Double_t) y*dssdH_E_Y[z][y][0];
                            nStripsY++;
                        }
                        if ((dssdH_E_Y[z][y][0]<=0||y==NumStrY-1)&&nStripsY>0) {//end of X cluster
                            if (nStripsY<=maxStrInCluster){ //if cluster of less than 3 strips
                                posY=E_Y_ch/E_Y;
                                //cout<<"x"<<posX <<"y"<<posY<<"e" <<(E_Y/E_X-1)*(E_Y/E_X-1)<<endl;
                                if (E_X>E_Y)
                                beta_dssd_pre.push_back(make_pair(1-E_Y/E_X,make_pair(make_pair(posX,make_pair(E_X,nStripsX)),make_pair(posY,make_pair(E_Y,nStripsY)))));
                                else
                                beta_dssd_pre.push_back(make_pair(1-E_X/E_Y,make_pair(make_pair(posX,make_pair(E_X,nStripsX)),make_pair(posY,make_pair(E_Y,nStripsY)))));
                            }
                            E_Y=0;
                            E_Y_ch=0;
                            nStripsY=0;
                            nClusterY++;
                        }
                    }
                    //ok!
                }
                E_X=0;
                E_X_ch=0;
                nStripsX=0;
                nClusterX++;
            }
        }

        //Additional step to filter all possible combination:
        if (nClusterX>0&&nClusterX<maxNCluster&&nClusterY>0&&nClusterY<maxNCluster){
            //Record number of cluster here
            //sort(beta_dssd_pre.begin(),beta_dssd_pre.end());
            for (ibeta_dssd_pre=beta_dssd_pre.begin();ibeta_dssd_pre<beta_dssd_pre.end();++ibeta_dssd_pre){
                if (mult_z<maxNposBeta){
                    AIDACluster *cluster=new AIDACluster;
                    double xx=ibeta_dssd_pre->second.first.first;
                    double yy=ibeta_dssd_pre->second.second.first;
                    if(ibeta_dssd_pre->second.first.second.second==1) xx=round(xx);
                    if(ibeta_dssd_pre->second.second.second.second==1) yy=round(yy);
                    cluster->SetHitPosition(xx,yy,z);
                    //cluster->SetHitPosition(ibeta_dssd_pre->second.first.first,ibeta_dssd_pre->second.second.first,z);
                    cluster->SetXEnergy(ibeta_dssd_pre->second.first.second.first);
                    cluster->SetYEnergy(ibeta_dssd_pre->second.second.second.first);
                    cluster->SetXMult(ibeta_dssd_pre->second.first.second.second);
                    cluster->SetYMult(ibeta_dssd_pre->second.second.second.second);
                    //Timing
                    int hitx = (int) round(ibeta_dssd_pre->second.first.first);
                    int hity = (int) round(ibeta_dssd_pre->second.second.first);
                    if (dssdH_T_X[z][hitx][0]>0&&dssdH_T_Y[z][hity][0]>0){
                        //! Take the ealiest time stamp
                        if (dssdH_T_X[z][hitx][0]>dssdH_T_Y[z][hity][0]) {
                            cluster->SetTimestamp(dssdH_T_Y[z][hity][0]);
                            cluster->SetFastTimestamp(dssdF_T_Y[z][hity][0]);
                        }else {
                            cluster->SetTimestamp(dssdH_T_X[z][hitx][0]);
                            cluster->SetFastTimestamp(dssdF_T_X[z][hitx][0]);
                        }
                        /*
                        //! if there is available fast time stamp
                        if (dssdF_T_X[z][hitx][0]==0||dssdF_T_Y[z][hity][0]==0){
                            cluster->SetFastTimestamp(dssdF_T_X[z][hitx][0] + dssdF_T_Y[z][hity][0]);
                        }
                        */

                    }else{
                        cout<<__PRETTY_FUNCTION__<<" something wrong!"<<endl;
                    }
                    this->AddCluster(cluster);
                    maxZ=z;
                }
                mult_z++;
            }
        }//end of additional step...
    }//end of loop on all dssd
    this->SetMaxZ(maxZ);
    if(this->GetNClusters()>0) return true;
    return false;
}

