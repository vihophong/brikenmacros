//---------------------------------------------------------------------
// Create the solids defining an Eurogam Phase-II detector
//---------------------------------------------------------------------
#include "DSSD.hh"
//#include "MyMaterials.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Trap.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "globals.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
//#include "PhysicalConstants.h"

#include <stdio.h>
//using namespace std;
DSSD::DSSD()
{
  //some default DSSD detector parameters
// fTotalGeL     = 69.00 * mm;  //was 70
  fTotalGeL     = 69.00 * mm;  //was 70
  fCrystalR     = 24.25 * mm;  //was 25
  fEndCap2Ge    = 20.00 * mm;  //Distance from Al outer face to Ge
  //added to fudge PTG's efficiency measurements for STUK
  fFudge = 5.0*mm;
  fEndCap2Ge += fFudge;
  
  fHoleR            =  5.5 * mm; //was 5.0
  fPassiveThick     =  0.5 * mm;  
  fContactThick     =  0.5 * mm;  
  fGapBetweenLeaves =  0.6 * mm;

  // dssd dimension
  fDSSDdx=75.6 * mm;
  fDSSDdy=75.6 * mm;
  fDSSDdz=1. * mm;
  fDSSDspacing=1. * cm;


  // PCB dimension
  fPCBdx=95. * mm;
  fPCBdy=95. * mm;
  fPCBdz=1. * mm;

  // Rod dimension
  fRodPosX=42 * mm;
  fRodPosY=42 * mm;

  fRodR=2* mm;
  fRodL=43* cm;
  fRodoffL=17* cm;

  // Snout dimension
    fSnoutdx=10* cm;
    fSnoutdy=10* cm;
    fSnoutdz=45* cm;
    fSnoutThickness=1* mm;
    fSnoutZ=17* cm;




  //---------------------------------
  //make the required materials and assign them
  //fMat = MyMaterials::GetInstance();
G4NistManager* nist = G4NistManager::Instance();
  //assign default materials.....

  airMaterial=nist->FindOrBuildMaterial("G4_AIR");
  vacuumMaterial = new G4Material("glc", 1., 1.01*g/mole, 
				      1.e-25*g/cm3,kStateGas,2.73*kelvin,3.e-18*pascal);


  SnoutMaterial=nist->FindOrBuildMaterial("G4_Al");
  DSSDMaterial=nist->FindOrBuildMaterial("G4_Si");
  RodMaterial=nist->FindOrBuildMaterial("G4_Ti");

  //
  // materials for rad-source setup
  //
  //from http://www.physi.uni-heidelberg.de/~adler/TRD/TRDunterlagen/RadiatonLength/tgc2.htm
  //Epoxy (for FR4 )
  G4Element* elH  = nist->FindOrBuildElement("H");
  G4Element* elC  = nist->FindOrBuildElement("C");
  G4Material* SiO2=nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

  G4Material* Epoxy = new G4Material("Epoxy" , 1.2*g/cm3, 2);
  Epoxy->AddElement(elH, 2);
  Epoxy->AddElement(elC, 2);

  //FR4 (Glass + Epoxy)
  PCBMaterial = new G4Material("FR4"  , 1.86*g/cm3, 2);
  PCBMaterial->AddMaterial(SiO2, 0.528);
  PCBMaterial->AddMaterial(Epoxy, 0.472);
  //create the solids.....
  CreateSolids();

}

//Destructor
DSSD::~DSSD() {
  //
}

void DSSD::SetPosition(G4ThreeVector thisPos) {
  
  position = thisPos*mm;
  G4cout << " ----> A Phase-II will be placed at " << position/mm << " mm" << G4endl;
 
}

void DSSD::SetRotation(G4RotationMatrix thisRot) {
  
  rotation = thisRot;
  
  //G4cout << " ----> Single Ge will be placed at " << position/mm << " mm" << G4endl;
 
}


//---------------------------------------------------------------------
// Create the solids defining Phase-II DSSDs
//---------------------------------------------------------------------
void  DSSD::CreateSolids()
{

  G4Box*        dssdBox[numDSSD];
  G4Box*        pcbBox[numDSSD];
  G4SubtractionSolid*        pcbBoxSub[numDSSD];
  G4SubtractionSolid*        pcbBoxSub1[numDSSD];
  G4SubtractionSolid*        pcbBoxSub2[numDSSD];
  G4SubtractionSolid*        pcbBoxSub3[numDSSD];
  G4SubtractionSolid*        pcbBoxSub4[numDSSD];
  G4Tubs*  pcbRodHole1[numDSSD];
  G4Tubs*  pcbRodHole2[numDSSD];
  G4Tubs*  pcbRodHole3[numDSSD];
  G4Tubs*  pcbRodHole4[numDSSD];

  G4Tubs*  pcbRod[4];


  char tempname[50];
  for (G4int i=0;i<numDSSD;i++){
      sprintf(tempname,"_%d",i);
      G4String identifier=(G4String) tempname;
      dssdBox[i] = new G4Box("dssdBox"+identifier,fDSSDdx/2.,fDSSDdy/2.,fDSSDdz/2.);
      //G4Tubs*  endCapBox=new G4Tubs("endCapBox", 0,  encapRodR, 0.5*encapRodL, 0.*degree, 360.*degree);
      logicDSSD[i] = new G4LogicalVolume(dssdBox[i],          //its solid
                          DSSDMaterial,         //its material
                          "logicDSSD"+identifier);        //its name
      pcbBox[i] = new G4Box("pcbBox"+identifier,fPCBdx/2.,fPCBdy/2.,fPCBdz/2.);

      pcbRodHole1[i]=new G4Tubs("pcbRodHole1"+identifier, 0,  2* mm, 0.5*fPCBdz, 0.*degree, 360.*degree);
      pcbRodHole2[i]=new G4Tubs("pcbRodHole2"+identifier, 0,  2* mm, 0.5*fPCBdz, 0.*degree, 360.*degree);
      pcbRodHole3[i]=new G4Tubs("pcbRodHole3"+identifier, 0,  2* mm, 0.5*fPCBdz, 0.*degree, 360.*degree);
      pcbRodHole4[i]=new G4Tubs("pcbRodHole4"+identifier, 0,  2* mm, 0.5*fPCBdz, 0.*degree, 360.*degree);

      pcbBoxSub1[i]=new G4SubtractionSolid("pcbBoxsub1"+identifier,pcbBox[i],pcbRodHole1[i],0,G4ThreeVector(fRodPosX,fRodPosY,0));
      pcbBoxSub2[i]=new G4SubtractionSolid("pcbBoxsub2"+identifier,pcbBoxSub1[i],pcbRodHole2[i],0,G4ThreeVector(-fRodPosX,-fRodPosY,0));
      pcbBoxSub3[i]=new G4SubtractionSolid("pcbBoxsub2"+identifier,pcbBoxSub2[i],pcbRodHole3[i],0,G4ThreeVector(fRodPosX,-fRodPosY,0));
      pcbBoxSub4[i]=new G4SubtractionSolid("pcbBoxsub2"+identifier,pcbBoxSub3[i],pcbRodHole4[i],0,G4ThreeVector(-fRodPosX,fRodPosY,0));

      pcbBoxSub[i]=new G4SubtractionSolid("pcbBoxsub"+identifier,pcbBoxSub4[i],dssdBox[i]);

      logicPCB[i] = new G4LogicalVolume(pcbBoxSub[i],          //its solid
                          PCBMaterial,         //its material
                          "logicPCB"+identifier);        //its name
  }
       pcbRod[0]=new G4Tubs("pcbRod0", 0,  fRodR, 0.5*fRodL, 0.*degree, 360.*degree);
       logicRod[0]=new G4LogicalVolume(pcbRod[0],
                                       RodMaterial,         //its material
                                       "logicRod0"
                                       );
       pcbRod[1]=new G4Tubs("pcbRod1", 0,  fRodR, 0.5*fRodL, 0.*degree, 360.*degree);
       logicRod[1]=new G4LogicalVolume(pcbRod[1],
                                       RodMaterial,         //its material
                                       "logicRod1"
                                       );
       pcbRod[2]=new G4Tubs("pcbRod2", 0,  fRodR, 0.5*fRodL, 0.*degree, 360.*degree);
       logicRod[2]=new G4LogicalVolume(pcbRod[2],
                                       RodMaterial,         //its material
                                       "logicRod2"
                                       );
       pcbRod[3]=new G4Tubs("pcbRod3", 0,  fRodR, 0.5*fRodL, 0.*degree, 360.*degree);
       logicRod[3]=new G4LogicalVolume(pcbRod[3],
                                       RodMaterial,         //its material
                                       "logicRod3"
                                       );

       G4Box*        SnoutBoxOut=new G4Box("SnoutBoxOut",fSnoutdx/2.,fSnoutdy/2.,fSnoutdz/2.);
       G4Box*        SnoutBoxIn=new G4Box("SnoutBoxIn",fSnoutdx/2.-fSnoutThickness,fSnoutdy/2.-fSnoutThickness,fSnoutdz/2.);

       G4SubtractionSolid* SnoutBoxSub=new G4SubtractionSolid("SnoutBoxSub",SnoutBoxOut,SnoutBoxIn,0,G4ThreeVector(0,0,0));
       logicSnout=new G4LogicalVolume(SnoutBoxSub,
                                       SnoutMaterial,         //its material
                                       "logicSnout"
                                       );


}


//------------------------------------------------------------------
void DSSD::Placement(G4int copyNo, G4VPhysicalVolume* physiMother)
{

    char tempname[50];
    for (G4int i=0;i<numDSSD;i++){
        sprintf(tempname,"_%d",i);
        G4String identifier=(G4String) tempname;
        physiDSSD[i] = new G4PVPlacement(&rotation,
                           G4ThreeVector(position.x(),position.y(),position.z()+fDSSDspacing*i-fDSSDspacing*numDSSD/2),
                           "physiDSSD"+identifier,       //its name
                           logicDSSD[i],//its logical volume
                           physiMother,         //its mother
                           true,               //no boolean operat
                           copyNo+i,             //copy number
                           true);              //overlap check
        physiPCB[i] = new G4PVPlacement(&rotation,
                           G4ThreeVector(position.x(),position.y(),position.z()+fDSSDspacing*i-fDSSDspacing*numDSSD/2),
                           "physiDSSD"+identifier,       //its name
                           logicPCB[i],//its logical volume
                           physiMother,         //its mother
                           true,               //no boolean operat
                           copyNo+i,             //copy number
                           true);              //overlap check
    }

    physiRod[0] = new G4PVPlacement(&rotation,
                       G4ThreeVector(position.x()+fRodPosX,position.y()+fRodPosY,position.z()+fRodoffL),
                       "physiRod0",       //its name
                       logicRod[0],//its logical volume
                       physiMother,         //its mother
                       true,               //no boolean operat
                       copyNo,             //copy number
                       true);              //overlap check
    physiRod[1] = new G4PVPlacement(&rotation,
                       G4ThreeVector(position.x()-fRodPosX,position.y()-fRodPosY,position.z()+fRodoffL),
                       "physiRod1",       //its name
                       logicRod[1],//its logical volume
                       physiMother,         //its mother
                       true,               //no boolean operat
                       copyNo+1,             //copy number
                       true);              //overlap check
    physiRod[2] = new G4PVPlacement(&rotation,
                       G4ThreeVector(position.x()+fRodPosX,position.y()-fRodPosY,position.z()+fRodoffL),
                       "physiRod2",       //its name
                       logicRod[2],//its logical volume
                       physiMother,         //its mother
                       true,               //no boolean operat
                       copyNo+2,             //copy number
                       true);              //overlap check
    physiRod[3] = new G4PVPlacement(&rotation,
                       G4ThreeVector(position.x()-fRodPosX,position.y()+fRodPosY,position.z()+fRodoffL),
                       "physiRod3",       //its name
                       logicRod[3],//its logical volume
                       physiMother,         //its mother
                       true,               //no boolean operat
                       copyNo+3,             //copy number
                       true);              //overlap check


    physiSnout = new G4PVPlacement(&rotation,
                       G4ThreeVector(position.x(),position.y(),position.z()+fSnoutZ),
                       "physiSnout",       //its name
                       logicSnout,//its logical volume
                       physiMother,         //its mother
                       true,               //no boolean operat
                       copyNo,             //copy number
                       true);              //overlap check

  G4VisAttributes* visAttDSSD = new G4VisAttributes( G4Colour(0.9,1.,0.3) );
  visAttDSSD->SetVisibility(true);
 visAttDSSD->SetForceSolid(true);

 G4VisAttributes* visAttPCB = new G4VisAttributes( G4Colour(0.,0.7,0.3) );
  visAttPCB->SetVisibility(true);
  visAttPCB->SetForceWireframe(true);
  
  G4VisAttributes* visAttRod = new G4VisAttributes( G4Colour(0.9,0.9,0.9) );
  visAttRod->SetForceSolid(true);
  visAttRod->SetVisibility(true);

  G4VisAttributes* visAttSnout = new G4VisAttributes( G4Colour(0.2,1.,1.) );
  visAttSnout->SetForceWireframe(true);
  visAttSnout->SetVisibility(true);

  for (G4int i=0;i<numDSSD;i++){
      logicDSSD[i]->SetVisAttributes(visAttDSSD);
      logicPCB[i]->SetVisAttributes(visAttPCB);
      if (i<4) logicRod[i]->SetVisAttributes(visAttRod);
  }
  logicSnout->SetVisAttributes(visAttSnout);


}
