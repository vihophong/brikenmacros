#ifndef DATASTRUCT_H
#define DATASTRUCT_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>

#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"

using namespace std;

//! INPUT STRUCTs
struct TreeData {
    ULong64_t ts;
    ULong64_t sts;
    Double_t tof;
    Double_t zet;
    Double_t aoq;
    Double_t f5x;
    Double_t f11x;
    Double_t f11y;
    Double_t f11dt;
    Double_t beta;
};

class AIDASimpleStruct : public TObject {
public:
  //! default constructor
  AIDASimpleStruct(){
    Clear();
  }
  virtual ~AIDASimpleStruct(){}

  //! Clear the music information
  virtual void Clear(){
      fid=0.;
      fts=0;
      fevt=0;
      fx=-9999.;
      fy=-9999.;
      fz=-1;
      fex=-9999.;
      fey=-9999.;
      fncx=0;
      fncy=0;
      fnx=0;
      fny=0;
      fnz=0;
      fniz=0;
      frflag=0;
      fdrflag=0;
      ftw=-9999.;
      fdtion=-9999.;
      fdz=0;
      fmult=0;
      fmindy=9999;
      fmindx=9999;
      fsumexyrank=0;
      fmaxelastdssd=-9999;
      fxts=0;
      fyts=0;
      fxtsmax=0;
      fytsmax=0;
      fminx = -11;
      fmaxx = -11;
      fminy = -11;
      fmaxy = -11;
      ftdiffmin= -99999;
      ftdiffmax= -99999;
  }

  //! copy cluster
  virtual void Copy(AIDASimpleStruct& obj){
      obj.SetID(fid);
      obj.SetEventNumber(fevt);
      obj.SetTimestamp(fts);
      obj.SetMult(fmult);
      obj.SetHitPosition(fx,fy,fz);
      obj.SetXEnergy(fex);
      obj.SetYEnergy(fey);
      obj.SetXClusterMult(fncx);
      obj.SetYClusterMult(fncy);
      obj.SetXMult(fnx);
      obj.SetYMult(fny);
      obj.SetZMult(fnz);
      obj.SetStripMultFlag(fniz);
      obj.SetRankingFlag(frflag);
      obj.SetEDiffRankingFlag(fdrflag);
      obj.SetTimeWidth(ftw);
      obj.SetDtIon(fdtion);
      obj.SetDZ(fdz);
      obj.SetMinimumDistanceX(fmindx);
      obj.SetMinimumDistanceY(fmindy);
      obj.SetSumEXYRank(fsumexyrank);
      obj.SetMaxELastDSSD(fmaxelastdssd);
      obj.SetXTimestamp(fxts);
      obj.SetYTimestamp(fyts);
      obj.SetXMaxTimestamp(fxtsmax);
      obj.SetYMaxTimestamp(fytsmax);
      obj.SetXMinPos(fminx);
      obj.SetXMaxPos(fmaxx);
      obj.SetYMinPos(fminy);
      obj.SetYMaxPos(fmaxy);
      obj.SetTimeDifferenceMax(ftdiffmax);
      obj.SetTimeDifferenceMin(ftdiffmin);
  }

  //! Set XY
  void SetHitPosition(double x,double y,short z){
      fx=x;
      fy=y;
      fz=z;
  }

  //! xy area of cluster
  void SetXMinPos(short minx) {fminx=minx;}
  void SetXMaxPos(short maxx) {fmaxx=maxx;}
  void SetYMinPos(short miny) {fminy=miny;}
  void SetYMaxPos(short maxy) {fmaxy=maxy;}

  //! Set ID of evnt
  void SetID(unsigned char id){fid = id;}

  //! Set Evt number
  void SetEventNumber(unsigned int evt){fevt = evt;}

  //! Set the X strips sum energy
  void SetXEnergy(double xenergy){fex = xenergy;}
  //! Set the Y strips sum energy
  void SetYEnergy(double yenergy){fey = yenergy;}

  //! Set the X strips ealiest timestamp
  void SetXTimestamp(unsigned long long xts){fxts = xts;}
  //! Set the Y strips ealiest timestamp
  void SetYTimestamp(unsigned long long yts){fyts = yts;}

  //! Set the X strips latest time
  void SetXMaxTimestamp(double xts){fxtsmax = xts;}


  //! Set time diffrernce between X-Y
  void SetTimeDifferenceMin(int tdiffmin){ftdiffmin = tdiffmin;}
  void SetTimeDifferenceMax(int tdiffmax){ftdiffmax = tdiffmax;}


  //! Set the Y strips latest time
  void SetYMaxTimestamp(double yts){fytsmax = yts;}

  //! Set the event multiplicity
  void SetMult(unsigned short mult){fmult = mult;}

  //! Set the X strips multiplicity
  void SetXMult(unsigned short xmult){fnx = xmult;}
  //! Set the Y strips multiplicity
  void SetYMult(unsigned short ymult){fny = ymult;}

  //! Set the X strips multiplicity
  void SetXClusterMult(unsigned short xcmult){fncx = xcmult;}
  //! Set the Y strips multiplicity
  void SetYClusterMult(unsigned short ycmult){fncy = ycmult;}


  //! Set the DSSD  multiplicity
  void SetZMult(unsigned short zmult){fnz = zmult;}

  //! Set the timestamp
  void SetTimestamp(unsigned long long ts){fts = ts;}

  //! Set Energy ratio ranking flag
  void SetRankingFlag(unsigned char rankingflag){frflag = rankingflag;}

  //! Set Energy difference ranking flag
  void SetEDiffRankingFlag(unsigned char rankingflag){fdrflag = rankingflag;}

  //! Set Minimum distance X
  void SetMinimumDistanceX(double mindx){fmindx=mindx;}
  //! Set Minimum distance Y
  void SetMinimumDistanceY(double mindy){fmindy=mindy;}

  //! Set Multiple hit in strip flag
  void SetStripMultFlag(unsigned short niz){fniz = niz;}

  //! Set Event Time Width
  void SetTimeWidth(double tw){ftw=tw;}

  //! Set time distance to ion
  void SetDtIon(double dtion){fdtion=dtion;}

  //! Set Z correction
  void SetDZ(unsigned short dz){fdz=dz;}

  //! Set EX+EY rank
  void SetSumEXYRank(unsigned short sumexyrank){fsumexyrank = sumexyrank;}

  //! Set Max Energy of last layer
  void SetMaxELastDSSD(double maxelastdssd){fmaxelastdssd = maxelastdssd;}

  //! Get ID of evnt
  unsigned char GetID(){return fid;}

  unsigned int GetEventNumber(){return fevt;}

  double GetHitPositionX(){return fx;}
  double GetHitPositionY(){return fy;}
  short GetHitPositionZ(){return fz;}

  //! Get xy area of cluster
  short GetMinHitPositionX(){return fminx;}
  short GetMaxHitPositionX(){return fmaxx;}
  short GetMinHitPositionY(){return fminy;}
  short GetMaxHitPositionY(){return fmaxy;}

  //! Get the X strips sum energy
  double GetXEnergy(){return fex;}
  //! Get the Y strips sum energy
  double GetYEnergy(){return fey;}

  //! Get the X strips ealiest timestamp
  unsigned long long GetXMinTimestamp(){return fxts;}
  //! Get the Y strips ealiest timestamp
  unsigned long long GetYMinTimestamp(){return fyts;}

  //! Get the X strips latest timestamp
  unsigned long long GetXMaxTimestamp(){return fxtsmax;}
  //! Get the Y strips latest timestamp
  unsigned long long GetYMaxTimestamp(){return fytsmax;}

  //! Get the time diffrerence of Y and X strips (min)
  int GetTimeDifferenceMin(){return ftdiffmin;}
  //! Get the time diffrerence of Y and X strips (max)
  int GetTimeDifferenceMax(){return ftdiffmax;}


  //! Get the event multiplicity
  unsigned short GetMultiplicity(){return fmult;}


  //! Get the X clusters multiplicity
  unsigned short GetXClusterMultiplicity(){return fncx;}
  //! Get the Y clusters multiplicity
  unsigned short GetYClusterMultiplicity(){return fncy;}

  //! Get the X strips multiplicity
  unsigned short GetXMultiplicity(){return fnx;}
  //! Get the Y strips multiplicity
  unsigned short GetYMultiplicity(){return fny;}

  //! Get the Z strips multiplicity
  unsigned short GetZMultiplicity(){return fnz;}

  //! Get the timestamp
  unsigned long long GetTimestamp(){return fts;}

  //! Get Energy ratio ranking flag
  unsigned char GetRankingFlag(){return frflag;}

  //! Get Energy difference ranking flag
  unsigned char GetEDiffRankingFlag(){return fdrflag;}

  //! Get minimum cluster distance X
  double GetMinimumDistanceX(){return fmindx;}
  //! Get minimum cluster distance Y
  double GetMinimumDistanceY(){return fmindy;}

  //! Get Multiple hit in strip flag
  unsigned short GetStripMultFlag(){return fniz;}

  //! Get time distance to ion
  double GetDtIon(){return fdtion;}

  //! Get even time width
  double GetTimeWidth(){return ftw;}

  //! Get Multiple hit in strip flag
  unsigned short GetDZ(){return fdz;}

  //! Get ex+ey rank
  unsigned short GetSumEXYRank(){return fsumexyrank;}

  //! Set Max Energy of last layer
  double GetMaxELastDSSD(){return fmaxelastdssd;}


  //! Printing information
  void Print(Option_t *option = "") const {
    cout <<"\tts "<< fts;
    cout << "\tX " << fx;
    cout << "\tY " << fy;
    cout << "\tZ " << fz;
    cout << "\tX energy " << fex;
    cout << "\tY energy " << fey;
    cout << "\tX multiplicity " << fnx;
    cout << "\tY multiplicity " << fny;
    cout << "\ttimestamp " << fts;
    cout << "\tranking flag " << frflag;
    return;
  }

protected:
  unsigned long long fts;

  unsigned int fevt;

  //! id :ion or beta
  unsigned char fid;

  //! event multiplicity
  unsigned short fmult;

  //! DSSDs coordinator
  double fx;
  double fy;
  short fz;

  //! area of cluster
  short fminx;
  short fmaxx;
  short fminy;
  short fmaxy;

  //! store the hit number
  //vector<short*> fhitsno;

  //! the energy lab system
  double fex;
  double fey;
  //! xy hit multiplicity
  unsigned short fnx;
  unsigned short fny;
  //! xy cluster multiplicity
  unsigned short fncx;
  unsigned short fncy;
  //! number of clusters
  unsigned short fnz;
  //! number of event with mult>1 for 1 dssd
  unsigned short fniz;
  //! ranking flag
  unsigned char frflag;
  //! ranking flag EY-EX
  unsigned char fdrflag;

  //!ealiest time stamp in x
  unsigned long long fxts;
  //!ealiest time stamp in y
  unsigned long long fyts;
  //!lastest time stamp in x
  unsigned long long fxtsmax;
  //!lastest time stamp in y
  unsigned long long fytsmax;

  //! time diffrence between y and x cluster
  int ftdiffmin;
  int ftdiffmax;

  //! ex+ey rank
  unsigned short fsumexyrank=0;

  //! minimum cluster distance x
  double fmindx;
  //! minimum cluster distance y
  double fmindy;

  //! event time width in us
  double ftw;
  //! time distance with prev ion in us
  double fdtion;

  //! z correction (if applicable)
  unsigned short fdz;

  //! max energy of the last layer (for light particle veto method)
  double fmaxelastdssd;

  /// \cond CLASSIMP
  ClassDef(AIDASimpleStruct,1);
  /// \endcond
};

class BELENHit : public TObject{
public:
    BELENHit(){
        Clear();
    }
    BELENHit(Double_t posx, Double_t posy, Double_t posz, short daqid, short id, short ring,short type, unsigned long long ts, int adc, int en, unsigned short hitsadded)
    {
        fdaqid = daqid;
        fid = id;
        fring = ring;
        ftype = type;
        fts = ts;
        fadc = adc;
        fen = en;
        fpos.SetXYZ(posx,posy,posz);
        fhitsadded = hitsadded;
    }
    virtual void Clear(){
        fid = -1;
        fring = -1;
        ftype = -1;
        fdaqid = -1;
        fpos.SetXYZ(-9999,-9999,-9999);
        frndpos.SetXYZ(-9999,-9999,-9999);
        fts = 0;
        fadc = -1;
        fen = -1;
        fhitsadded = 0;

        //for veto
        fdvetotime = -999999.;
        ff11time = -999999.;
        fvetotime = -999999.;
    }
    virtual void Copy(BELENHit& obj){
        obj.SetID(fid);
        obj.SetRing(fring);
        obj.SetType(ftype);
        obj.SetDaqID(fdaqid);
        obj.SetTimestamp(fts);
        obj.SetPos(fpos.X(),fpos.Y(),fpos.Z());
        obj.SetRndPos(frndpos.X(),frndpos.Y(),frndpos.Z());
        obj.SetADC(fadc);
        obj.SetEnergy(fen);
        obj.SetHitsAdded(fhitsadded);
        //for veto
        obj.SetF11Time(ff11time);
        obj.SetDownstreamVetoTime(fdvetotime);
        obj.SetFinalVetoTime(fvetotime);
    }

    //! Set the energy
    void SetEnergy(double energy){fen = energy;}

    //! Set the raw ADC value
    void SetADC(int adc){fadc = adc;}

    //! Set the counter ID
    void SetID(short id){fid = id;}
    //! Set the counter daq ID
    void SetDaqID(short daqid){fdaqid = daqid;}
    //! Set the counter ring
    void SetRing(short ring){fring = ring;}
    //! Set the counter type
    void SetType(short type){ftype = type;}

    //! Set the timestamp
    void SetTimestamp(unsigned long long int ts){fts = ts;}

    //! Set the He3 position
    void SetPos(Double_t x, Double_t y, Double_t z){fpos.SetXYZ(x,y,z);}

    //! Set the He3 position
    void SetRndPos(Double_t x, Double_t y, Double_t z){frndpos.SetXYZ(x,y,z);}

    //! Set current hits
    void SetHitsAdded(unsigned short hitsadded){fhitsadded = hitsadded;}


    //! for veto
    void SetF11Time(double f11time){ff11time = f11time;}
    void SetDownstreamVetoTime(double dvetotime){fdvetotime = dvetotime;}
    void SetFinalVetoTime(double vetotime){fvetotime = vetotime;}


    //! Get the ID
    short GetID(){return fid;}
    //! Get the ID
    short GetDaqID(){return fdaqid;}
    //! Get the ring (my precious!)
    short GetMyPrecious(){return fring;}
    //! Get the type
    short GetType(){return ftype;}

    //! Get the energy
    double GetEnergy(){return fen;}
    //! Get the timestamp
    unsigned long long int GetTimestamp(){return fts;}
    //! Get the raw ADC value
    int GetADC(){return fadc;}
    //! Get 3He position
    TVector3 GetPosition(){return fpos;}
    TVector3 GetRndPosition(){return frndpos;}

    //! Get current hits
    unsigned short GetHitsAdded(){return fhitsadded;}


    //! for veto
    double GetF11Time(){return ff11time;}
    double GetDownstreamVetoTime(){return fdvetotime;}
    double GetFinalVetoTime(){return fvetotime;}

    //! Printing information
    void Print(Option_t *option = "") const {
      cout << "ID " << fid;
      cout << "daq ID " << fdaqid;
      cout << "ring " << fring;
      cout << "type" <<ftype;
      cout << "\tX pos " << fpos.X();
      cout << "\tY pos " << fpos.Y();
      cout << "\tZ pos " << fpos.Z();
      cout << "\tX pos " << frndpos.X();
      cout << "\tY pos " << frndpos.Y();
      cout << "\tZ pos " << frndpos.Z();
      cout << "\tadc " << fadc;
      cout << "\tenergy " << fen;
      cout << "\ttimestamp " << fts;
      cout << "\thits added " << fhitsadded << endl;
      return;
    }

protected:

    //! current hits
    unsigned short fhitsadded;
    //! Position of 3He counter
    TVector3 fpos;
    //! Position of 3He counter with random generator
    TVector3 frndpos;
    //! daq ID number of 3He counter
    short fdaqid;
    //! physical ID number of 3He counter
    short fid;
    //! ring number of the 3He counter
    short fring;
    //! type of tube: 0: riken, 1: upc 1 inch 2: ornl 1 inch 3: ornl 2 inch
    short ftype;

    //! ADC value
    int fadc;
    //! Energy calibrated value
    double fen;
    //! timestamp value
    unsigned long long fts;

    //! special added
    double ff11time;//in us
    double fdvetotime;//in us
    double fvetotime; //in us


    /// \cond CLASSIMP
    ClassDef(BELENHit,1);
    /// \endcond
    ///
};

class CloverHit : public TObject
{
public:
    CloverHit(){
        Clear();
    }
    virtual void Clear(){
        fid = -1;
        fdaqid = -1;
        fpos.SetXYZ(-1,-1,-1);
        fts = 0;
        fadc = -1;
        fen = -1;
        fhitsadded = 0;
        fclover = -1;
        fcloverleaf = -1;
    }

    virtual void Copy(CloverHit& obj){
        obj.SetID(fid);
        obj.SetClover(fclover);
        obj.SetCloverLeaf(fcloverleaf);
        obj.SetDaqID(fdaqid);
        obj.SetTimestamp(fts);
        obj.SetPos(fpos.X(),fpos.Y(),fpos.Z());
        obj.SetADC(fadc);
        obj.SetEnergy(fen);
        obj.SetHitsAdded(fhitsadded);
    }

    //! Set the energy
    void SetEnergy(double energy){fen = energy;}

    //! Set the raw ADC value
    void SetADC(int adc){fadc = adc;}

    //! Set the ID
    void SetID(short id){fid = id;}
    //! Set the Daq ID
    void SetDaqID(short daqid){fdaqid = daqid;}
    //! Set the clover
    void SetClover(short clover){fclover = clover;}
    //! Set the clover leaf
    void SetCloverLeaf(short cloverleaf){fcloverleaf = cloverleaf;}

    //! Set the timestamp
    void SetTimestamp(unsigned long long int ts){fts = ts;}

    //! Set the Clover position
    void SetPos(Double_t x, Double_t y, Double_t z){fpos.SetXYZ(x,y,z);}
    //! Set current hits
    void SetHitsAdded(unsigned short hitsadded){fhitsadded = hitsadded;}


    //! Get the ID
    short GetID(){return fid;}
    //! Get the Daq ID
    short GetDaqID(){return fdaqid;}
    //! Get the clover
    short GetClover(){return fclover;}
    //! Get the clover leaf
    short GetCloverLeaf(){return fcloverleaf;}

    //! Get the energy
    double GetEnergy(){return fen;}
    //! Get the timestamp
    unsigned long long GetTimestamp(){return fts;}
    //! Get the raw ADC value
    int GetADC(){return fadc;}

    //! Get clover position
    TVector3 GetPosition(){return fpos;}

    //! Get current hits
    unsigned short GetHitsAdded(){return fhitsadded;}

    //! Printing information
    void Print(Option_t *option = "") const {
      cout << "ID " << fid;
      cout << "daq ID " << fdaqid;
      cout << "clover leaf"<<fcloverleaf;
      cout << "clover no."<<fclover;
      cout << "\tX pos " << fpos.X();
      cout << "\tY pos " << fpos.Y();
      cout << "\tZ pos " << fpos.Z();
      cout << "\tadc " << fadc;
      cout << "\tenergy " << fen;
      cout << "\ttimestamp " << fts;
      cout << "\thits added " << fhitsadded << endl;
      return;
    }

protected:
    //! current hits (for addback)
    unsigned short fhitsadded;
    //! Position of the clover which respect to the center hole (in cm)
    TVector3 fpos;
    //! ID number of the clover counter
    short fid;
    short fdaqid;
    //! clover number 1: left side clover, 2: right side clover
    short fclover;
    //! leaf in the clover 1: black 2: red 3: green 4: blue
    short fcloverleaf;

    //! ADC value
    int fadc;
    //! Energy calibrated value
    double fen;
    //! timestamp value
    unsigned long long fts;

    /// \cond CLASSIMP
    ClassDef(CloverHit,1);
    /// \endcond
    ///
};

//! OUTPUT STRUCTs

#endif
