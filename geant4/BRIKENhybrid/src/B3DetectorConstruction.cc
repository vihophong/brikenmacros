//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B3DetectorConstruction.cc 71323 2013-06-13 16:54:23Z gcosmo $
//
/// \file B3DetectorConstruction.cc
/// \brief Implementation of the B3DetectorConstruction class

#include "B3DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4RunManager.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4UnionSolid.hh"

#include <fstream>
#include <stdio.h>

#include "Clover.hh"
using namespace std;	 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void B3DetectorConstruction::ring(G4String infilename,G4String outfilename,G4double f=1)
{
	ifstream infile(infilename);
	ofstream outfile(outfilename);
	G4int nRing,nTubes;
	infile>>nRing;
	infile>>nTubes;
	G4int ringID,tubeType,nTubesofRing,nTube,fdetGroup;
	G4double ringRadius,offDeg,posZ,xx,yy,zz,ringRadiusPre,ringRadiusPre2,afterReadRingRadius;
	G4int k=0;
	outfile<<nTubes<<"\n";
	for(G4int i=0; i<nRing; i++){
	        ringRadiusPre2=afterReadRingRadius;
		ringRadiusPre=ringRadius;
		infile>>ringID>>tubeType>>ringRadius>>nTubesofRing>>nTube>>offDeg>>posZ>>fdetGroup;
		afterReadRingRadius=ringRadius;
		if (fdetGroup==1){
		  ringRadius=f*(ringRadius-firstRingRadiusPre)+firstRingRadius;
		}
                if (fdetGroup>1){
		  ringRadius=f*(ringRadius-ringRadiusPre2)+ringRadiusPre;
		}
		G4double dphi=twopi/nTubesofRing;
		for (G4int j=0;j<nTube;j++){
			G4double phi=offDeg*twopi/360.+j*dphi;
			xx=ringRadius*cos(phi);
			yy=ringRadius*sin(phi);
			zz=posZ;
			outfile<<k<<" "<<tubeType<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<fdetGroup<<"\n";
			k++;
		}
	}
	infile.close();
	outfile.close();
}

void B3DetectorConstruction::ring2(G4String infilename,G4String outfilename)
{
	ifstream infile(infilename);
	ofstream outfile(outfilename);
	G4int nRing,nTubes;
	infile>>nRing;
	infile>>nTubes;
	G4int ringID,tubeType,nTubesofRing,nTube,fdetGroup;
	G4double ringRadius,offDeg,posZ,xx,yy,zz,ringRadiusPre,ringRadiusPre2,afterReadRingRadius;
	G4int k=0;
	outfile<<nTubes<<"\n";
	for(G4int i=0; i<nRing; i++){
		infile>>ringID>>tubeType>>ringRadius>>nTubesofRing>>nTube>>offDeg>>posZ>>fdetGroup;
		G4double dphi=twopi/nTubesofRing;
		for (G4int j=0;j<nTube;j++){
			G4double phi=offDeg*twopi/360.+j*dphi;
			xx=ringRadius*cos(phi);
			yy=ringRadius*sin(phi);
			zz=posZ;
			outfile<<k<<" "<<tubeType<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<fdetGroup<<"\n";
			k++;
		}
	}
	infile.close();
	outfile.close();
}

void B3DetectorConstruction::ring3(G4String infilename,G4String outfilename)
{
	ifstream infile(infilename);
	ofstream outfile(outfilename);
	G4int nRings,totalTubes;
	infile>>nRings>>totalTubes;

	outfile<<totalTubes<<"\n";

	G4int ringID,nGroup,nTubes;
	G4double ringRadius,xx,yy,zz,offDeg,startX,startY;
	G4int gg=0;
	for (G4int i=0;i<nRings;i++){
	  infile>>ringID>>nGroup>>nTubes;	  
	  G4int tubeType[100];G4int startN[100];G4int nTubesonRing[100];G4double ringZ[100];G4int detGroup[100];G4int isOPP[100];
	  for (G4int k=0;k<nGroup;k++){
	    infile>>tubeType[k]>>nTubesonRing[k]>>startN[k]>>ringZ[k]>>detGroup[k]>>isOPP[k];	    
	  }

	  ringRadius=sqrt(secondX0*secondX0+secondY0*secondY0)+sFactor*i;
	  offDeg=twopi/4-acos(secondY0/ringRadius);
	  G4double dphi=((twopi/4-offDeg)*2)/(nTubes-1);
	  for (G4int j=0;j<nTubes;j++){
	    G4double phi=offDeg+j*dphi;
	    xx=ringRadius*cos(phi);
	    yy=ringRadius*sin(phi);

	    for (G4int ma=0;ma<nGroup;ma++){
	      if (j>=startN[ma]&&j<startN[ma]+nTubesonRing[ma]){		  
		if (isOPP[ma]==1){
		  outfile<<gg<<" "<<tubeType[ma]<<" "<<xx<<" "<<-yy<<" "<<ringZ[ma]<<" "<<detGroup[ma]<<"\n";
		}else{
		  outfile<<gg<<" "<<tubeType[ma]<<" "<<xx<<" "<<yy<<" "<<ringZ[ma]<<" "<<detGroup[ma]<<"\n";
		}
		gg++;
	      }
	    }
	  }
	}

	infile.close();
	outfile.close();
}

void B3DetectorConstruction::rect1(G4String infilename,G4String outfilename,G4double h0=1,G4double f=1)
{
	ifstream infile(infilename);
	ofstream outfile(outfilename);
	G4int nRings,totalTubes;
	infile>>nRings;
	infile>>totalTubes;
	G4int ringID,tubeType,nTube,nTubeonRing,startNum,fdetGroup, isEdge;
	G4double ringRadius,ringZ,xx,yy,zz,preRingRadius,preRingRadius2,afterRingRadius;
	outfile<<totalTubes<<"\n";
	G4double offset,interval;
	for (G4int i=0;i<nRings;i++){
	  preRingRadius=ringRadius;
	  preRingRadius2=afterRingRadius;
	  infile>>ringID>>tubeType>>nTube>>nTubeonRing>>startNum>>ringRadius>>ringZ>>fdetGroup>>isEdge;
	  firstRingRadiusPre=ringRadius;
	  afterRingRadius=ringRadius;

	  if (fdetGroup==0) {
	    ringRadius=ringRadius*h0;
	    firstRingRadius=ringRadius;
	  }else{
	    ringRadius=f*(ringRadius-preRingRadius2)+preRingRadius;
	    firstRingRadius=ringRadius;
	  }
	  startNoFix:
	  offset=ringRadius/(sqrt(2)-1+nTube/4)-ringRadius/sqrt(2);
	  interval=sqrt(2)*ringRadius/(sqrt(2)-1+nTube/4);
	  G4int gg=0;
	  if (isEdge==1&&nTube!=4){
	    offset=-ringRadius/sqrt(2);
	    interval=sqrt(2)*ringRadius/(nTube/4-1);
	  }else if (isEdge==1&&nTube==4){
	    offset=-ringRadius/sqrt(2);
	  }
	  for (G4int k=0;k<4;k++){
	    for (G4int j=0;j<nTube/4;j++){	      
	      if (k==0){
		xx=offset+j*interval;
		yy=ringRadius/sqrt(2);
	      }else if (k==2){
		xx=-offset-j*interval;
		yy=-ringRadius/sqrt(2);
	      }else if (k==1){
		yy=-offset-j*interval;
		xx=ringRadius/sqrt(2);
	      }else if (k==3){
		yy=offset+j*interval;
		xx=-ringRadius/sqrt(2);
	      }
	      zz=ringZ;
	      if (gg>=startNum&&gg<startNum+nTubeonRing) outfile<<k<<" "<<tubeType<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<fdetGroup<<"\n";
	      gg++;
	    }
	  }

	}
	infile.close();
	outfile.close();
}


void B3DetectorConstruction::rect5(G4String infilename,G4String outfilename,G4double h0=1,G4double f=1)
{
	ifstream infile(infilename);
	ofstream outfile(outfilename);
	G4int nRings,totalTubes,nTubesFirstRing;
	G4double h0absolute,f0absolute;
	infile>>nRings>>totalTubes>>nTubesFirstRing>>h0absolute>>f0absolute;
	G4double h0f=h0absolute*h0;
	G4double f0f=f0absolute*f;
	sFactor=f0f+5.08;
	h0f=h0f+2.54/2;
	f0f=f0f+2.54;
	outfile<<totalTubes<<"\n";

	G4int ringID,nGroup;
	G4double ringRadius,xx,yy,zz,temp,temp2;
	G4int gg=0;
	for (G4int i=0;i<nRings;i++){
	  infile>>ringID>>nGroup;	  
	  G4int tubeType[100];G4int startN[100];G4int nTubesonRing[100];G4double ringZ[100];G4int detGroup[100];G4int isOPP[100];
	  for (G4int k=0;k<nGroup;k++){
	    infile>>tubeType[k]>>nTubesonRing[k]>>startN[k]>>ringZ[k]>>detGroup[k]>>isOPP[k];
	    
	  }
	  if (nTubesFirstRing%2!=0){//1 tube at center
	    for (G4int j=0;j<nTubesFirstRing;j++){
	      xx=-(nTubesFirstRing-1)/2*f0f+f0f*j+f0f*i/2;
	      yy=h0f+f0f*sin(twopi*60/360)*i;

	      if (i==0&&j==0) {
		secondX0=-xx+2.54/2+5.08/2+h0*f*0.05+0.4;//truc x
		secondY0=yy-2.54/2+5.08/2;//truc y
	      }

	      for (G4int ma=0;ma<nGroup;ma++){
		if (j>=startN[ma]&&j<startN[ma]+nTubesonRing[ma]){
		  
		  if (isOPP[ma]==1){
		    outfile<<gg<<" "<<tubeType[ma]<<" "<<xx<<" "<<-yy<<" "<<ringZ[ma]<<" "<<detGroup[ma]<<"\n";
		  }else{
		    outfile<<gg<<" "<<tubeType[ma]<<" "<<xx<<" "<<yy<<" "<<ringZ[ma]<<" "<<detGroup[ma]<<"\n";
		  }
		  gg++;
		}
	      }
	    }
	  }else{//no tube at center
	    for (G4int j=0;j<nTubesFirstRing;j++){
	      xx=-(nTubesFirstRing/2-0.5)*f0f+f0f*j+f0f*i/2;
	      yy=h0f+f0f*sin(twopi*60/360)*i;

	      if (i==0&&j==0) {
		secondX0=-xx+2.54/2+5.08/2+h0*f*0.05+0.4;//truc x
		secondY0=yy-2.54/2+5.08/2;//truc y
	      }

	      for (G4int ma=0;ma<nGroup;ma++){
		if (j>=startN[ma]&&j<startN[ma]+nTubesonRing[ma]){		  
		  if (isOPP[ma]==1){
		    outfile<<gg<<" "<<tubeType[ma]<<" "<<xx<<" "<<-yy<<" "<<ringZ[ma]<<" "<<detGroup[ma]<<"\n";
		  }else{
		    outfile<<gg<<" "<<tubeType[ma]<<" "<<xx<<" "<<yy<<" "<<ringZ[ma]<<" "<<detGroup[ma]<<"\n";
		  }
		  gg++;
		}
	      }
	    }
	  }//end if


	}
	infile.close();
	outfile.close();
}


void B3DetectorConstruction::rect4(G4String infilename,G4String outfilename,G4double h0=1,G4double f=1)
{
	ifstream infile(infilename);
	ofstream outfile(outfilename);
	G4int nRings,totalTubes;
	infile>>nRings;
	infile>>totalTubes;
	G4int ringID,tubeType,nTube,nTubeonRing,startNum,fdetGroup, isEdge,isOPP;
	G4double ringRadius,ringZ,xx,yy,zz,preRingRadius,preRingRadius2,afterRingRadius,fAngle;
	outfile<<totalTubes<<"\n";
	G4double offset,interval;
	for (G4int i=0;i<nRings;i++){
	  preRingRadius=ringRadius;
	  preRingRadius2=afterRingRadius;
	  infile>>ringID>>tubeType>>nTube>>nTubeonRing>>startNum>>ringRadius>>ringZ>>fdetGroup>>isEdge;
	  fAngle=fAngle*twopi/360;
	  firstRingRadiusPre=ringRadius;
	  afterRingRadius=ringRadius;

	  if (fdetGroup==0) {
	    ringRadius=ringRadius*h0;
	    firstRingRadius=ringRadius;
	  }else{
	    ringRadius=f*(ringRadius-preRingRadius2)+preRingRadius;
	    firstRingRadius=ringRadius;
	  }
	  startNoFix:

	  offset=ringRadius*sin(fAngle)/(nTube-1)-ringRadius*sin(fAngle);
	  interval=ringRadius*sin(fAngle)*2/nTube;
	  G4cout<<"aaa="<<interval<<G4endl;
	  G4int gg=0;
	  if (isEdge==1){
	    offset=-ringRadius*sin(fAngle);
	    interval=ringRadius*sin(fAngle)*2/(nTube-1);
	  }

	  for (G4int j=0;j<nTube;j++){
	    xx=offset+j*interval;
	    if (isOPP==1) yy=-ringRadius*cos(fAngle); else yy=ringRadius*cos(fAngle);
	    zz=ringZ;
	    if (gg>=startNum&&gg<startNum+nTubeonRing) outfile<<gg<<" "<<tubeType<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<fdetGroup<<"\n";
	    gg++;
	  }	 
	}
	infile.close();
	outfile.close();
}


void B3DetectorConstruction::rect2(G4String infilename,G4String outfilename,G4double h0=1.)
{
	ifstream infile(infilename);
	ofstream outfile(outfilename);
	G4int nRings,nTubes;
	infile>>nRings;
	infile>>nTubes;
	G4int ringID,tubeType,nTubeonRing,fdetGroup,rotate;
	G4double ringRadius,ringZ,startD,DeltaD,xx,yy,zz;
	G4int k=0;
	outfile<<nTubes<<"\n";
	for(G4int i=0; i<nRings; i++){
	  infile>>ringID>>nTubeonRing>>tubeType>>ringRadius>>ringZ>>startD>>DeltaD>>fdetGroup>>rotate;
	  if (fdetGroup==0&&h0!=1.) {
	    firstRingRadiusPre=abs(ringRadius/cos(twopi*45/360));
	    firstRingRadius=abs(h0*ringRadius/cos(twopi*45/360));
	    startD=-1.5*((firstRingRadius*sin(twopi*45/360)*2-5.*h0)/3);
	    DeltaD=(firstRingRadius*sin(twopi*45/360)*2-5.*h0)/3;
	    ringRadius=h0*ringRadius;
	  }else{
	    firstRingRadiusPre=abs(ringRadius/cos(twopi*45/360));
	    firstRingRadius=abs(h0*ringRadius/cos(twopi*45/360));
	  }
		for (G4int j=0;j<nTubeonRing;j++){
		  if (rotate==0){
			xx=startD+j*DeltaD;
			yy=ringRadius;
			zz=ringZ;
			outfile<<k<<" "<<tubeType<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<fdetGroup<<"\n";
		  }else{
			yy=startD+j*DeltaD;
			xx=ringRadius;
			zz=ringZ;
			outfile<<k<<" "<<tubeType<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<fdetGroup<<"\n";
		  }
			k++;
		}
	}
	infile.close();
	outfile.close();
}


void B3DetectorConstruction::readConfiguration(G4String infilename)
{
	ifstream infile(infilename);
	infile>>nofTubes;
	G4int tempd;
	for(G4int i=0; i<nofTubes; i++){
	  infile>>tempd>>tubeGroup[i]>>x[i]>>y[i]>>z[i]>>detGroup[i];
        }
        infile.close();
}

void B3DetectorConstruction::readConfiguration2(G4String infilename1,G4String infilename2)
{
  G4int nofTubes1,nofTubes2;
	ifstream infile1(infilename1);
	infile1>>nofTubes1;
	ifstream infile2(infilename2);
	infile2>>nofTubes2;
	nofTubes=nofTubes1+nofTubes2;

	G4int tempd;
	for(G4int i=0; i<nofTubes1; i++){
	  infile1>>tempd>>tubeGroup[i]>>x[i]>>y[i]>>z[i]>>detGroup[i];
        }
	for(G4int i=nofTubes1; i<nofTubes; i++){
	  infile2>>tempd>>tubeGroup[i]>>x[i]>>y[i]>>z[i]>>detGroup[i];
        }
        infile1.close();
        infile2.close();
}

void B3DetectorConstruction::readConfiguration3(G4String infilename1,G4String infilename2,G4String infilename3)
{
  G4int nofTubes1,nofTubes2,nofTubes3;
	ifstream infile1(infilename1);
	infile1>>nofTubes1;
	ifstream infile2(infilename2);
	infile2>>nofTubes2;
	ifstream infile3(infilename3);
	infile3>>nofTubes3;
	nofTubes=nofTubes1+nofTubes2+nofTubes3;

	G4int tempd;
	for(G4int i=0; i<nofTubes1; i++){
	  infile1>>tempd>>tubeGroup[i]>>x[i]>>y[i]>>z[i]>>detGroup[i];
        }
	for(G4int i=nofTubes1; i<nofTubes1+nofTubes2; i++){
	  infile2>>tempd>>tubeGroup[i]>>x[i]>>y[i]>>z[i]>>detGroup[i];
        }
	for(G4int i=nofTubes1+nofTubes2; i<nofTubes; i++){
	  infile3>>tempd>>tubeGroup[i]>>x[i]>>y[i]>>z[i]>>detGroup[i];
        }
        infile1.close();
        infile2.close();
        infile3.close();
}



void B3DetectorConstruction::setConfig(G4double h0,G4double f)
{
  fh0=h0;
  ff=f;
  //rect3("test2.brr","test2.br",h0,1);
  //ring("test1.brr","test1.br",f);
  //readConfiguration2("test1.br","test2.br");
  //rect5("test3.brr","test3.br",fh0,ff);
  //ring3("test4.brr","test4.br");
  readConfiguration3("test3.br","test4.br","test5.br");
  if(!((h0==1)&&(f==1))) G4RunManager::GetRunManager()->ReinitializeGeometry();
}

B3DetectorConstruction::B3DetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(true),
  nofTubes(0),
  fDetectorMessenger(0)
{

  fDetectorMessenger = new DetectorMessenger(this);
	detGroup=new G4int[1000];
	tubeGroup=new G4int[1000];
	x=new G4double[1000];
	y=new G4double[1000];
	z=new G4double[1000];

	secondX0=0;
	secondY0=0;
	sFactor=0;
//Read input Geo--
	/*fh0=1;
	ff=1;        
	G4String infilename;
	G4cout<<"Enter input Geometry:"<<G4endl;
	G4cin>>infilename;
	G4String outfilename=infilename+"r";
	ring2(infilename,outfilename);
	readConfiguration(outfilename);
	*/
	//readConfiguration2("test1.br","test2.br");
	setConfig(1,1);

//--
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::~B3DetectorConstruction()
{
  delete fDetectorMessenger;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B3DetectorConstruction::Construct()
{  
  // Cleanup old geometry                                                                                                                      
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  

  //
	//Temperature for all
	G4double temperature=298.15*kelvin;
	//-------------------Material used---------------
	G4NistManager* nist = G4NistManager::Instance();
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");

  //
  //  define Vacuum
  //
  G4Material* Vacuum = new G4Material("Galactic", 1., 1.01*g/mole, 
				      1.e-25*g/cm3,kStateGas,2.73*kelvin,3.e-18*pascal);


	G4Isotope* He3 = new G4Isotope("He3",2,3, 3.016*g/mole);
	G4Element* elHe3  = new G4Element("He3 enriched", "3He",1);

	elHe3->AddIsotope(He3, 100.*perCent);	
	G4Element* elC  = nist->FindOrBuildElement("C");
	//G4Element* elC = new G4Element("C_ele","C_ele",6.0,12.0*g/mole);
	G4Element* elTSH = new G4Element( "TS_H_of_Polyethylene","H_POLYETHYLENE" , 1.0 , 1.0079*g/mole ); 
	G4Material* CH2mat= new G4Material( "Polyethylene_TS" ,0.95*g/cm3 ,2);//, kStateUndefined, temperature); 
	CH2mat -> AddElement(elTSH,2);
	CH2mat -> AddElement(elC,1);

	G4Element* elS  = nist->FindOrBuildElement("S");
	G4Element* elMn  = nist->FindOrBuildElement("Mn");
	G4Element* elCr  = nist->FindOrBuildElement("Cr");
	G4Element* elP  = nist->FindOrBuildElement("P");
	G4Element* elSi  = nist->FindOrBuildElement("Si");
	G4Element* elNi  = nist->FindOrBuildElement("Ni");
	G4Element* elN  = nist->FindOrBuildElement("N");
	G4Element* elFe  = nist->FindOrBuildElement("Fe");
	G4Material* STL= new G4Material( "StainlessSteels" ,8.03*g/cm3 ,9);

	STL -> AddElement(elC,0.08*perCent);
	STL -> AddElement(elS,0.03*perCent);
	STL -> AddElement(elMn,2.*perCent);
	STL -> AddElement(elCr,20.*perCent);
	STL -> AddElement(elP,0.045*perCent);
	STL -> AddElement(elSi,0.75*perCent);
	STL -> AddElement(elNi,10.5*perCent);
	STL -> AddElement(elN,0.1*perCent);
	STL -> AddElement(elFe,66.495*perCent);

	//-------------------Define Geometry---------------

	G4double moderator_sizeX=90*cm;
	G4double moderator_sizeY=90*cm;
	G4double moderator_sizeZ=90*cm;
  //     
  // World
  //
  G4double world_sizeX = 1.2*moderator_sizeX;
  G4double world_sizeY = 1.2*moderator_sizeY;
  G4double world_sizeZ = 1.2*moderator_sizeZ;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ); //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        Vacuum,         //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);       // checking overlaps 
        //
	//Moderator Polyethylene
	//
	G4Box* solidModerator =    
    new G4Box("Moderator",                       //its name
       0.5*moderator_sizeX, 0.5*moderator_sizeY, 0.5*moderator_sizeZ); //its size
  G4LogicalVolume* logicModerator =                         
    new G4LogicalVolume(solidModerator,          //its solid
                        CH2mat,         //its material
                        "ModeratorLogic");            //its name
	// Place Moderator in the world
  G4VPhysicalVolume* physModerator = 
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    "ModeratorLogic",               //its name
                    logicModerator,            //its logical volume
                    physWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);         // checking overlaps
               

  //
  // AIDA
  //
  G4double AIDA_XY = 11.6*cm;
  G4double AIDA_Z = 90*cm;  

  G4double AIDA_XY2 = 14*cm;
  G4double AIDA_Z2 = 15.2*cm;     

	G4Box* solidAIDA =    
    new G4Box("AIDA1",                       //its name
       0.5*AIDA_XY, 0.5*AIDA_XY, 0.5*AIDA_Z); //its size
	G4Box* solidAIDA2 =    
    new G4Box("AIDA2",                       //its name
       0.5*AIDA_Z, 0.5*AIDA_XY, 0.5*AIDA_XY); //its size   
	G4Box* solidAIDA3 =    
    new G4Box("AIDA3",                       //its name
       0.5*AIDA_Z2, 0.5*AIDA_XY2, 0.5*AIDA_XY2); //its size
	G4Box* solidAIDA4 =    
    new G4Box("AIDA4",                       //its name
       0.5*AIDA_Z2, 0.5*AIDA_XY2, 0.5*AIDA_XY2); //its size
G4UnionSolid* solidAIDA_Sum=new G4UnionSolid("solidAIDA",solidAIDA,solidAIDA2,0,G4ThreeVector(0,0,0));
 G4UnionSolid* solidAIDA_Sum2=new G4UnionSolid("solidAIDAsum2",solidAIDA_Sum,solidAIDA3,0,G4ThreeVector((24+15.2/2+11.6/2)*cm,0,0));
 G4UnionSolid* solidAIDA_Sum3=new G4UnionSolid("solidAIDAsum3",solidAIDA_Sum2,solidAIDA4,0,G4ThreeVector(-(24+15.2/2+11.6/2)*cm,0,0));
  G4LogicalVolume* logicAIDA =                         
    new G4LogicalVolume(solidAIDA_Sum3,        //its solid
                        Vacuum,         //its material
                        "AIDALV");        //its name
               
  //
  // place AIDA in the Moderator
  //                    
  G4VPhysicalVolume* physAIDA=new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
									  "AIDA",               //its name
                    logicAIDA,            //its logical volume                    
                    physModerator,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);         // checking overlaps 

      
  //
  //Place Clover in Moderator
  //
        
  Clover* clv=new Clover();
  G4ThreeVector clvPos=G4ThreeVector((-11.6+5.8/2)*cm,0,-5.8/2*cm);
  G4RotationMatrix clvRot;
  clvRot.set(0,0,0);
  clvRot.rotateY(90.*deg);
  //clvRot.invert();
  clv->SetPosition(clvPos);
  clv->SetRotation(clvRot);
  //clv->CreateSolids();
  clv->Placement(0,physAIDA);
  
          
  //2nd one
  Clover* clv2=new Clover();
  clvPos.setX((11.6-5.8/2)*cm);
  //clvRot.invert();
  clvRot.rotateY(180.*deg);
  clv2->SetPosition(clvPos);
  clv2->SetRotation(clvRot);
  clv2->Placement(0,physAIDA);





  /*
      
  Clover* clv=new Clover();
  G4ThreeVector clvPos=G4ThreeVector(-fakeAIDA_X/2+2.8*cm,0.*cm,-5.5/2*cm);
  G4RotationMatrix clvRot;
  clvRot.set(0,0,0);
  //clvRot.invert();
  clvRot.rotateY(270.*deg);
  clv->SetPosition(clvPos);
  clv->SetRotation(clvRot);
  //clv->CreateSolids();
  clv->Placement(0,fakephysAIDA);
  */
  

  //
  // Place He3 rings in the Moderator
  //	
	G4double temp=3.01603/(temperature/kelvin)/82.05746;
	G4Material* He3matp10 = new G4Material("He3matp10", (10.0*temp)*g/cm3, 1, kStateGas, temperature, 10.0*atmosphere);
	He3matp10->AddElement(elHe3,100*perCent);
	G4Material* He3matp8 = new G4Material("He3matp8", (8.0*temp)*g/cm3, 1, kStateGas, temperature, 8.0*atmosphere);
	He3matp8->AddElement(elHe3,100*perCent);
	G4Material* He3matp5 = new G4Material("He3matp5", (5.0*temp)*g/cm3, 1, kStateGas, temperature, 5.13*atmosphere);
	He3matp5->AddElement(elHe3,100*perCent);
	G4Material* He3matp4 = new G4Material("He3matp4", (4.0*temp)*g/cm3, 1, kStateGas, temperature, 4.0*atmosphere);
	He3matp4->AddElement(elHe3,100*perCent);
	
  //
  // visualization attributes
  //

  G4VisAttributes* Red= 
    new G4VisAttributes(G4Colour(1.0,0.0,0.0,0.3));  // Red
    Red->SetVisibility(true);
    Red->SetForceSolid(true);

  G4VisAttributes* Green= 
    new G4VisAttributes(G4Colour(0.0,1.0,0.0,0.05));  // Green  
    Green->SetVisibility(true);
    Green->SetForceSolid(true);
    
  G4VisAttributes* Blue= 
    new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.3));  // Blue
    Blue->SetVisibility(true);
    Blue->SetForceSolid(true);
        
  G4VisAttributes* Yellow= 
    new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.3));  // Yellow
    Yellow->SetVisibility(true);
    Yellow->SetForceSolid(true);
    
  G4VisAttributes* Turquoise= 
    new G4VisAttributes(G4Colour(0.0,0.8,0.8,0.3));  // Turquoise
    Turquoise->SetVisibility(true);
    Turquoise->SetForceSolid(true);
    
  G4VisAttributes* Magenta= 
    new G4VisAttributes(G4Colour(1.0,0.0,1.0,0.3));  // Magenta
    Magenta->SetVisibility(true);
    Magenta->SetForceSolid(true);
    
  G4VisAttributes* Grey= 
    new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.08));  // Gray
    Grey->SetVisibility(true);
    Grey->SetForceSolid(true);


//---- Place tube in moderator----

for (G4int i=0;i<nofTubes;i++)
{
//--- Define Detector Group --
//--------------------------------------------------------------
	G4LogicalVolume* He3Tube;
	G4ThreeVector position;                    
	switch (tubeGroup[i]){
	case 0: //UPC
	{
	  G4LogicalVolume* He3Tube3=B3DetectorConstruction::DefineTubeGroup(i, 60.*cm,63.72*cm, 2.438/2*cm, 2.54/2*cm, 0.36*cm, He3matp8, STL);
		//DefineTubeGroup(G4int nTube, G4double He3Z,G4double Outer_Z, G4double He3R, G4double Outer_R, G4double He3OffsetPos, G4Material* He3mat, G4Material* Outermat)
		He3Tube=He3Tube3;
		He3Tube->SetVisAttributes(Red);
		position=G4ThreeVector(x[i]*cm,y[i]*cm,z[i]-0.36*cm);
		break;
	}
	case 1: //GSI
	{
	  G4LogicalVolume* He3Tube3=B3DetectorConstruction::DefineTubeGroup(i, 60.*cm,63.72*cm, 2.438/2*cm, 2.54/2*cm, 0.36*cm, He3matp10, STL);
		He3Tube=He3Tube3;
		He3Tube->SetVisAttributes(Blue);
		position=G4ThreeVector(x[i]*cm,y[i]*cm,z[i]-0.36*cm);
		break;
	}
	case 2: //ORNL 1 inch
	{
	  G4LogicalVolume* He3Tube3=B3DetectorConstruction::DefineTubeGroup(i, 60.96*cm,66.09*cm, 2.438/2*cm, 2.54/2*cm, 0.*cm, He3matp10,  STL);
		He3Tube=He3Tube3;
		He3Tube->SetVisAttributes(Yellow);
		position=G4ThreeVector(x[i]*cm,y[i]*cm,z[i]);
		break;
	}
	case 3: //ORNL 2 inch
	{
	  G4LogicalVolume* He3Tube3=B3DetectorConstruction::DefineTubeGroup(i, 60.96*cm,66.09*cm, 4.978/2*cm, 5.08/2*cm, 0.61*cm, He3matp10,  STL);
		He3Tube=He3Tube3;
		He3Tube->SetVisAttributes(Turquoise);
		position=G4ThreeVector(x[i]*cm,y[i]*cm,z[i]-0.61*cm);
		break;
	}
	case 4: //RIKEN (plus Z value)
	{
	  G4LogicalVolume* He3Tube3=B3DetectorConstruction::DefineTubeGroup(i, 30.0*cm,35.212*cm, 2.438/2*cm, 2.54/2*cm, -0.54*cm, He3matp5,  STL);
		He3Tube=He3Tube3;
		He3Tube->SetVisAttributes(Magenta);
		if (z[i]==0){
		  position=G4ThreeVector(x[i]*cm,y[i]*cm,z[i]+0.54*cm);
		}else{
		  position=G4ThreeVector(x[i]*cm,y[i]*cm,z[i]*cm);
		}
		break;
	}
	case 5: //RIKEN opposite (minus Z value)
	{
	  G4LogicalVolume* He3Tube3=B3DetectorConstruction::DefineTubeGroup(i, 30.0*cm,35.212*cm, 2.438/2*cm, 2.54/2*cm, 0.54*cm, He3matp5,  STL);
		He3Tube=He3Tube3;
		He3Tube->SetVisAttributes(Magenta);
		if (z[i]==0){
		  position=G4ThreeVector(x[i]*cm,y[i]*cm,z[i]-0.54*cm);
		}else{
		  position=G4ThreeVector(x[i]*cm,y[i]*cm,z[i]*cm);
		}
		break;
	}
	case 6: //JINR
	{
	  G4LogicalVolume* He3Tube3=B3DetectorConstruction::DefineTubeGroup(i, 50.0*cm,55.0*cm, 2.9/2*cm, 3.0/2*cm, 0.*cm, He3matp4,  STL);	 

		He3Tube=He3Tube3;
		He3Tube->SetVisAttributes(Grey);
		position=G4ThreeVector(x[i]*cm,y[i]*cm,z[i]*cm);
		break;
	}

	}


  new G4PVPlacement(0,
											position,             //position
                      "PhysicsTube",             //its name
                      He3Tube,            //its logical volume
                      physModerator,             //its mother  volume
                      false,                 //no boolean operation
                      i,                 //copy number
                      fCheckOverlaps);       // checking overlaps
}

 logicModerator->SetVisAttributes(Grey);                                        
 logicAIDA->SetVisAttributes(Green);
 //fakelogicAIDA->SetVisAttributes(Green);
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

  // Print materials
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl; 

  //always return the physical World
  //

  return physWorld;
}


G4LogicalVolume* B3DetectorConstruction::DefineTubeGroup(G4int nTube, G4double He3Z,G4double Outer_Z, G4double He3R, G4double Outer_R,G4double He3OffsetPos, G4Material* He3mat, G4Material* Outermat)
{

	char tempname[50];
	sprintf(tempname,"_%d",nTube);
	G4String identifier=(G4String) tempname;


  //     
  // Define Outer(TUBE)
  //
    G4Tubs* solidHe3Outer =
    new G4Tubs("He3out", 0, Outer_R, 0.5*Outer_Z, 0., twopi);
                     
    G4LogicalVolume* tubeLogic = 
    new G4LogicalVolume(solidHe3Outer,          //its solid
                        Outermat,         //its material
                        "Outertube"+identifier);        //its name

  //     
  // Define He3s
  //
  G4Tubs* solidHe3 =
    new G4Tubs("He3", 0, He3R, 0.5*He3Z, 0., twopi);
  G4LogicalVolume* logicHe3 = 
    new G4LogicalVolume(solidHe3,          //its solid
                        He3mat,         //its material
                        "Tube"+identifier);        //its name

	//
	// Define Inactive region
	//
  G4Tubs* wholeHe3= new G4Tubs("wholeHe3",0., He3R, 0.5*Outer_Z-(Outer_R-He3R),0., twopi);
    G4LogicalVolume* logicInactive = 
    new G4LogicalVolume(wholeHe3,          //its solid
                        He3mat,         //its material
                        "Inactive"+identifier);        //its name


    
    new G4PVPlacement(0,
		      G4ThreeVector(0.*mm,0.*mm,He3OffsetPos),             //position
                      logicHe3,            //its logical volume
		      "He3Physical"+identifier,             //its name			
                      logicInactive,             //its mother  volume
                      false,                 //no boolean operation
                      0,                 //copy number
                      fCheckOverlaps);       // checking overlaps
    
    new G4PVPlacement(0,
		      G4ThreeVector(0.*mm,0.*mm,0*mm),             //position
                      logicInactive,            //its logical volume
                      "InactivePV"+identifier,             //its name
                      tubeLogic,             //its mother  volume
                      false,                 //no boolean operation
                      0,                 //copy number
                      fCheckOverlaps);       // checking overlaps
    

    return tubeLogic;
}

