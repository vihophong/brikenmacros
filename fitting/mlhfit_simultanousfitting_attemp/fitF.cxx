/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "fitF.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(fitF) 

 fitF::fitF(const char *name, const char *title, 
                        RooAbsReal& _x,
//                        RooAbsCategory& _y,
                        RooAbsReal& _neueff,
                        RooAbsReal& _p0,
                        RooAbsReal& _p1,
                        RooAbsReal& _p2,
                        RooAbsReal& _p3,
                        RooAbsReal& _p4,
                        RooAbsReal& _p5,
                        RooAbsReal& _p6,
                        RooAbsReal& _p7,
                        RooAbsReal& _p8,
                        RooAbsReal& _p9,
                        RooAbsReal& _p10,
                        RooAbsReal& _p11,
                        RooAbsReal& _p12,
                        RooAbsReal& _p13,
                        RooAbsReal& _p14,
                        RooAbsReal& _p15,
                        RooAbsReal& _p16,
                        RooAbsReal& _p17,
                        RooAbsReal& _p18,
                        RooAbsReal& _p19,
                        RooAbsReal& _p20,
                        RooAbsReal& _p21,
                        RooAbsReal& _p22,
                        RooAbsReal& _p23) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
// y("y","y",this,_y),
   neueff("neueff","neueff",this,_neueff),
   p0("p0","p0",this,_p0),
   p1("p1","p1",this,_p1),
   p2("p2","p2",this,_p2),
   p3("p3","p3",this,_p3),
   p4("p4","p4",this,_p4),
   p5("p5","p5",this,_p5),
   p6("p6","p6",this,_p6),
   p7("p7","p7",this,_p7),
   p8("p8","p8",this,_p8),
   p9("p9","p9",this,_p9),
   p10("p10","p10",this,_p10),
   p11("p11","p11",this,_p11),
   p12("p12","p12",this,_p12),
   p13("p13","p13",this,_p13),
   p14("p14","p14",this,_p14),
   p15("p15","p15",this,_p15),
   p16("p16","p16",this,_p16),
   p17("p17","p17",this,_p17),
   p18("p18","p18",this,_p18),
   p19("p19","p19",this,_p19),
   p20("p20","p20",this,_p20),
   p21("p21","p21",this,_p21),
   p22("p22","p22",this,_p22),
   p23("p23","p23",this,_p23)
 { 
 } 


 fitF::fitF(const fitF& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
//   y("y",this,other.y),
   neueff("neueff",this,other.neueff),
   p0("p0",this,other.p0),
   p1("p1",this,other.p1),
   p2("p2",this,other.p2),
   p3("p3",this,other.p3),
   p4("p4",this,other.p4),
   p5("p5",this,other.p5),
   p6("p6",this,other.p6),
   p7("p7",this,other.p7),
   p8("p8",this,other.p8),
   p9("p9",this,other.p9),
   p10("p10",this,other.p10),
   p11("p11",this,other.p11),
   p12("p12",this,other.p12),
   p13("p13",this,other.p13),
   p14("p14",this,other.p14),
   p15("p15",this,other.p15),
   p16("p16",this,other.p16),
   p17("p17",this,other.p17),
   p18("p18",this,other.p18),
   p19("p19",this,other.p19),
   p20("p20",this,other.p20),
   p21("p21",this,other.p21),
   p22("p22",this,other.p22),
   p23("p23",this,other.p23)
 { 
 } 



 Double_t fitF::evaluate() const 
 { 
     Double_t par0[22];
     par0[0]=p0;
     par0[1]=p1;
     par0[2]=p2;
     par0[3]=p3;
     par0[4]=p4;
     par0[5]=p5;
     par0[6]=p6;
     par0[7]=p7;
     par0[8]=p8;
     par0[9]=p9;
     par0[10]=p10;
     par0[11]=p11;
     par0[12]=p12;
     par0[13]=p13;
     par0[14]=p14;
     par0[15]=p15;
     par0[16]=p16;
     par0[17]=p17;
     par0[18]=p18;
     par0[19]=p19;
     par0[20]=p20;
     par0[21]=p23;

     Double_t par1[22];
     par1[0]=p0;
     par1[1]=p1;
     par1[2]=p2;
     par1[3]=p3;
     par1[4]=p4;
     par1[5]=p5;
     par1[6]=p6;
     par1[7]=p7;
     par1[8]=p8;
     par1[9]=p9;
     par1[10]=p10;
     par1[11]=p11;
     par1[12]=p12;
     par1[13]=p13;
     par1[14]=p14;
     par1[15]=p15;
     par1[16]=p16;
     par1[17]=p17;
     par1[18]=p18;
     par1[19]=p19;
     par1[20]=p21;
     par1[21]=p23;

     Double_t par2[22];
     par2[0]=p0;
     par2[1]=p1;
     par2[2]=p2;
     par2[3]=p3;
     par2[4]=p4;
     par2[5]=p5;
     par2[6]=p6;
     par2[7]=p7;
     par2[8]=p8;
     par2[9]=p9;
     par2[10]=p10;
     par2[11]=p11;
     par2[12]=p12;
     par2[13]=p13;
     par2[14]=p14;
     par2[15]=p15;
     par2[16]=p16;
     par2[17]=p17;
     par2[18]=p18;
     par2[19]=p19;
     par2[20]=p22;
     par2[21]=p23;

     Double_t ret=0;

     ret = fcn_gen(x,par0);
     return ret ;
 } 

 //! Global function
 Double_t fitF::fcn_gen(Double_t t, Double_t *par) const{
     Int_t npaths;
     Int_t ndecay[kmaxpaths];
     Int_t decaymap[kmaxpaths][kmaxndecay];
     Int_t nneu[kmaxpaths][kmaxndecay];
     //! define decay map
     npaths=12;
     //! path 1(go for ri2)
     ndecay[0]=2;
     decaymap[0][0]=1;decaymap[0][1]=2;
     nneu[0][0]=0;
     //! path 2(go for ri4)
     ndecay[1]=3;
     decaymap[1][0]=1;decaymap[1][1]=2;decaymap[1][2]=4;
     nneu[1][0]=0;nneu[1][1]=0;
     //! path 3(go for ri7)
     ndecay[2]=4;
     decaymap[2][0]=1;decaymap[2][1]=2;decaymap[2][2]=4;decaymap[2][3]=7;
     nneu[2][0]=0;nneu[2][1]=0;nneu[2][2]=0;
     //! path4(go for ri3)
     ndecay[3]=2;
     decaymap[3][0]=1;decaymap[3][1]=3;
     nneu[3][0]=1;
     //!path5(go for ri6)
     ndecay[4]=2;
     decaymap[4][0]=1;decaymap[4][1]=6;
     nneu[4][0]=2;

     //! path6(go for ri5-route 1)
     ndecay[5]=3;
     decaymap[5][0]=1;decaymap[5][1]=2;decaymap[5][2]=5;
     nneu[5][0]=0;nneu[5][1]=1;
     //! path7(go for ri5-route 2)
     ndecay[6]=3;
     decaymap[6][0]=1;decaymap[6][1]=3;decaymap[6][2]=5;
     nneu[6][0]=1;nneu[6][1]=0;

     //! path8 (go for ri9-route 1)
     ndecay[7]=3;
     decaymap[7][0]=1;decaymap[7][1]=3;decaymap[7][2]=9;
     nneu[7][0]=1;nneu[7][1]=1;
     //! path9 (go for ri9-route 2)
     ndecay[8]=3;
     decaymap[8][0]=1;decaymap[8][1]=6;decaymap[8][2]=9;
     nneu[8][0]=2;nneu[8][1]=0;

     //! path10 (go for ri8-route 1)
     ndecay[9]=4;
     decaymap[9][0]=1;decaymap[9][1]=2;decaymap[9][2]=4;decaymap[9][3]=8;
     nneu[9][0]=0;nneu[9][1]=0;nneu[9][2]=1;

     //! path11 (go for ri8-route 2)
     ndecay[10]=4;
     decaymap[10][0]=1;decaymap[10][1]=2;decaymap[10][2]=5;decaymap[10][3]=8;
     nneu[10][0]=0;nneu[10][1]=1;nneu[10][2]=0;

     //! path12 (go for ri8-route 3)
     ndecay[11]=4;
     decaymap[11][0]=1;decaymap[11][1]=3;decaymap[11][2]=5;decaymap[11][3]=8;
     nneu[11][0]=1;nneu[11][1]=0;nneu[11][2]=0;

     //Double_t bkg=par[knri*2+2];
     //if (t<0) return bkg;
     //Double_t returnval=bkg;

     //if (t<0) return 0;
     Double_t returnval=0;

     Double_t* pn=&par[knri];
     Double_t* lamda=par;
     Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
     p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
     Double_t N0=par[knri*2+1]/par[0];

     //! Parent nuclei
     returnval+=lamda[0]*N0*exp(-lamda[0]*t);

     for (Int_t i=0;i<npaths;i++){
         returnval+=lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,t);
     }
     return returnval;
 }

 Double_t fitF::fcn_gen_w1neutron(Double_t t, Double_t *par) const{
     Int_t npaths;
     Int_t ndecay[kmaxpaths];
     Int_t decaymap[kmaxpaths][kmaxndecay];
     Int_t nneu[kmaxpaths][kmaxndecay];
     //! define decay map
     npaths=12;
     //! path 1(go for ri2)
     ndecay[0]=2;
     decaymap[0][0]=1;decaymap[0][1]=2;
     nneu[0][0]=0;
     //! path 2(go for ri4)
     ndecay[1]=3;
     decaymap[1][0]=1;decaymap[1][1]=2;decaymap[1][2]=4;
     nneu[1][0]=0;nneu[1][1]=0;
     //! path 3(go for ri7)
     ndecay[2]=4;
     decaymap[2][0]=1;decaymap[2][1]=2;decaymap[2][2]=4;decaymap[2][3]=7;
     nneu[2][0]=0;nneu[2][1]=0;nneu[2][2]=0;
     //! path4(go for ri3)
     ndecay[3]=2;
     decaymap[3][0]=1;decaymap[3][1]=3;
     nneu[3][0]=1;
     //!path5(go for ri6)
     ndecay[4]=2;
     decaymap[4][0]=1;decaymap[4][1]=6;
     nneu[4][0]=2;

     //! path6(go for ri5-route 1)
     ndecay[5]=3;
     decaymap[5][0]=1;decaymap[5][1]=2;decaymap[5][2]=5;
     nneu[5][0]=0;nneu[5][1]=1;
     //! path7(go for ri5-route 2)
     ndecay[6]=3;
     decaymap[6][0]=1;decaymap[6][1]=3;decaymap[6][2]=5;
     nneu[6][0]=1;nneu[6][1]=0;

     //! path8 (go for ri9-route 1)
     ndecay[7]=3;
     decaymap[7][0]=1;decaymap[7][1]=3;decaymap[7][2]=9;
     nneu[7][0]=1;nneu[7][1]=1;
     //! path9 (go for ri9-route 2)
     ndecay[8]=3;
     decaymap[8][0]=1;decaymap[8][1]=6;decaymap[8][2]=9;
     nneu[8][0]=2;nneu[8][1]=0;

     //! path10 (go for ri8-route 1)
     ndecay[9]=4;
     decaymap[9][0]=1;decaymap[9][1]=2;decaymap[9][2]=4;decaymap[9][3]=8;
     nneu[9][0]=0;nneu[9][1]=0;nneu[9][2]=1;

     //! path11 (go for ri8-route 2)
     ndecay[10]=4;
     decaymap[10][0]=1;decaymap[10][1]=2;decaymap[10][2]=5;decaymap[10][3]=8;
     nneu[10][0]=0;nneu[10][1]=1;nneu[10][2]=0;

     //! path12 (go for ri8-route 3)
     ndecay[11]=4;
     decaymap[11][0]=1;decaymap[11][1]=3;decaymap[11][2]=5;decaymap[11][3]=8;
     nneu[11][0]=1;nneu[11][1]=0;nneu[11][2]=0;


     Double_t randcoinf=par[knri*2+3];
     //Double_t bkg=par[knri*2+2];
     //if (t<0) return bkg;
     //Double_t returnval=bkg;
     //if (t<0) return 0;
     Double_t returnval=0;

     Double_t* pn=&par[knri];
     Double_t* lamda=par;
     Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
     p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
     Double_t N0=par[knri*2+1]/par[0];

     for (Int_t i=0;i<npaths;i++){
         returnval+=neueff*pn[decaymap[i][ndecay[i]-1]-1]*lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,t)*(1-randcoinf);
     }

     //! decay with neutron part
     returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*t)*(1-randcoinf);
     returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*t)*(1-randcoinf);
     return returnval;
 }


 Double_t fitF::fcn_gen_w2neutron(Double_t t, Double_t *par) const{
     Int_t npaths;
     Int_t ndecay[kmaxpaths];
     Int_t decaymap[kmaxpaths][kmaxndecay];
     Int_t nneu[kmaxpaths][kmaxndecay];
     //! define decay map
     npaths=12;
     //! path 1(go for ri2)
     ndecay[0]=2;
     decaymap[0][0]=1;decaymap[0][1]=2;
     nneu[0][0]=0;
     //! path 2(go for ri4)
     ndecay[1]=3;
     decaymap[1][0]=1;decaymap[1][1]=2;decaymap[1][2]=4;
     nneu[1][0]=0;nneu[1][1]=0;
     //! path 3(go for ri7)
     ndecay[2]=4;
     decaymap[2][0]=1;decaymap[2][1]=2;decaymap[2][2]=4;decaymap[2][3]=7;
     nneu[2][0]=0;nneu[2][1]=0;nneu[2][2]=0;
     //! path4(go for ri3)
     ndecay[3]=2;
     decaymap[3][0]=1;decaymap[3][1]=3;
     nneu[3][0]=1;
     //!path5(go for ri6)
     ndecay[4]=2;
     decaymap[4][0]=1;decaymap[4][1]=6;
     nneu[4][0]=2;

     //! path6(go for ri5-route 1)
     ndecay[5]=3;
     decaymap[5][0]=1;decaymap[5][1]=2;decaymap[5][2]=5;
     nneu[5][0]=0;nneu[5][1]=1;
     //! path7(go for ri5-route 2)
     ndecay[6]=3;
     decaymap[6][0]=1;decaymap[6][1]=3;decaymap[6][2]=5;
     nneu[6][0]=1;nneu[6][1]=0;

     //! path8 (go for ri9-route 1)
     ndecay[7]=3;
     decaymap[7][0]=1;decaymap[7][1]=3;decaymap[7][2]=9;
     nneu[7][0]=1;nneu[7][1]=1;
     //! path9 (go for ri9-route 2)
     ndecay[8]=3;
     decaymap[8][0]=1;decaymap[8][1]=6;decaymap[8][2]=9;
     nneu[8][0]=2;nneu[8][1]=0;

     //! path10 (go for ri8-route 1)
     ndecay[9]=4;
     decaymap[9][0]=1;decaymap[9][1]=2;decaymap[9][2]=4;decaymap[9][3]=8;
     nneu[9][0]=0;nneu[9][1]=0;nneu[9][2]=1;

     //! path11 (go for ri8-route 2)
     ndecay[10]=4;
     decaymap[10][0]=1;decaymap[10][1]=2;decaymap[10][2]=5;decaymap[10][3]=8;
     nneu[10][0]=0;nneu[10][1]=1;nneu[10][2]=0;

     //! path12 (go for ri8-route 3)
     ndecay[11]=4;
     decaymap[11][0]=1;decaymap[11][1]=3;decaymap[11][2]=5;decaymap[11][3]=8;
     nneu[11][0]=1;nneu[11][1]=0;nneu[11][2]=0;


     Double_t randcoinf=par[knri*2+3];
     //Double_t bkg=par[knri*2+2];
     //if (t<0) return bkg;
     //Double_t returnval=bkg;
     //if (t<0) return 0;
     Double_t returnval=0;

     Double_t* pn=&par[knri];
     Double_t* lamda=par;
     Double_t p2n[knri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
     p2n[0]=par[knri*2];//special for p2n of parrent, to be develope later
     Double_t N0=par[knri*2+1]/par[0];

     for (Int_t i=0;i<npaths;i++){
         returnval+=neueff*pn[decaymap[i][ndecay[i]-1]-1]*lamda[decaymap[i][ndecay[i]-1]-1]*corefcn(ndecay[i],decaymap[i],nneu[i],pn,p2n,lamda,N0,t)*randcoinf;
     }

     //! decay with neutron part
     returnval+=neueff*pn[0]*lamda[0]*N0*exp(-lamda[0]*t)*randcoinf;
     returnval+=2*(neueff*(1-neueff))*p2n[0]*lamda[0]*N0*exp(-lamda[0]*t)*randcoinf;
     //! decay with 2 neutron part
     returnval+=neueff*neueff*p2n[0]*lamda[0]*N0*exp(-lamda[0]*t);
     return returnval;
 }



 Double_t fitF::corefcn(Int_t ndecay,Int_t* decaymap,Int_t* nneu, Double_t* b1n,Double_t* b2n,Double_t* lamda,Double_t N0,Double_t t) const{
     Double_t fcnret=0;

     Double_t factor1=1.;
     //only parrent decay p2n
     for (int i=0;i<ndecay-1;i++){
         if (nneu[i]==0){
             factor1=factor1 * (1-b1n[decaymap[i]-1]-b2n[decaymap[i]-1])*lamda[decaymap[i]-1];
         }else if (nneu[i]==1){
             factor1=factor1 * b1n[decaymap[i]-1]*lamda[decaymap[i]-1];
         }else{
             factor1=factor1 * b2n[decaymap[i]-1]*lamda[decaymap[i]-1];
         }
     }
     Double_t factor2=0;
     for (int i=0;i<ndecay;i++){
         Double_t factor2i=exp(-lamda[decaymap[i]-1]*t);
         Double_t factor2ij=1;
         for (int j=0;j<ndecay;j++)
             if (j!=i) factor2ij=factor2ij*(lamda[decaymap[j]-1]-lamda[decaymap[i]-1]);
         factor2=factor2+factor2i/factor2ij;
     }

     fcnret=factor1*N0*factor2;
     return fcnret;
 }






