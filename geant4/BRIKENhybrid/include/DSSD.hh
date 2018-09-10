#ifndef DSSD_H
#define DSSD_H 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
//class G4Box;
//class G4Tubs;
//class G4Polycone;
class G4UnionSolid;
class G4SubtractionSolid;
//class G4IntersectionSolid;
//class G4Polyhedra;
class G4LogicalVolume;
class G4VPhysicalVolume;
class MyDetectorMessenger;
//class MyMaterials;

#define numDSSD 6

class DSSD {
  
public:
  DSSD();
  DSSD(G4String);
  ~DSSD();

public:
  void SetPosition( G4ThreeVector );
  void SetRotation( G4RotationMatrix );
  void Placement(G4int, G4VPhysicalVolume*);

  inline G4double GetTaperedEndCapL() {return fEndCapTaperL/mm;};

private:
  //General materials....
//  MyMaterials*   fMat;
  G4Material* endCapMaterial;
  G4Material* contactMaterial;
  G4Material* geMaterial;

  G4Material* airMaterial;
  G4Material* vacuumMaterial;

  G4Material* DSSDMaterial;
  G4Material* PCBMaterial;
  G4Material* RodMaterial;
  G4Material* SnoutMaterial;

  G4ThreeVector        position;
  G4RotationMatrix     rotation;

  G4double             fDSSDdx;
  G4double             fDSSDdy;
  G4double             fDSSDdz;
  G4double             fPCBdx;
  G4double             fPCBdy;
  G4double             fPCBdz;

  G4double             fRodPosX;
  G4double             fRodPosY;


  G4double             fRodL;
  G4double             fRodoffL;
  G4double             fRodR;


  G4double             fSnoutdx;
  G4double             fSnoutdy;
  G4double             fSnoutdz;
  G4double             fSnoutThickness;
  G4double             fSnoutZ;




  G4double             fDSSDspacing;



  G4double             fTotalGeL;
  G4double             fHoleR;


  G4double             fCrystalR;
  G4double             fContactThick;
  G4double             fPassiveThick;
  G4double             fEndCapTaperL;
  G4double             fEndCapBoxL; //Add
  G4double             fEndCapThickness;
  G4double             fEndCap2Ge;
  G4double             fFudge;
  G4double             fVacuumPosZ;
  G4double             fContact_dZ;
  G4double             fGeLeafPosZ;
  G4double             fGapBetweenLeaves;
  G4double             fGeLeaf_dX;
  G4double             fHole_dX;
  G4double             fHole_dY;
  //Add
  //G4double rodR;
  //G4double rodL;

  G4UnionSolid*        solidEndCap;
  G4UnionSolid*        solidVacuum;
  G4UnionSolid*        solidGeLeaf;
  G4UnionSolid*        solidPassivated;
  G4UnionSolid*        solidContact;
  G4UnionSolid*        solidBoreHole;

  G4LogicalVolume*     logicEndCap;
  G4VPhysicalVolume*   physiEndCap;
  G4LogicalVolume*     logicVacuum;
  G4VPhysicalVolume*   physiVacuum;
  G4LogicalVolume*     logicGeLeaf[4];
  G4VPhysicalVolume*   physiGeLeaf[4];
  G4LogicalVolume*     logicPassivated[4];
  G4VPhysicalVolume*   physiPassivated[4];
  G4LogicalVolume*     logicContact[4];
  G4VPhysicalVolume*   physiContact[4];
  G4LogicalVolume*     logicBoreHole[4];
  G4VPhysicalVolume*   physiBoreHole[4];


  G4LogicalVolume*   logicDSSD[numDSSD];
  G4VPhysicalVolume*   physiDSSD[numDSSD];

  G4LogicalVolume*   logicRod[4];
  G4VPhysicalVolume*   physiRod[4];

  G4LogicalVolume*   logicPCB[numDSSD];
  G4VPhysicalVolume*   physiPCB[numDSSD];


  G4LogicalVolume*   logicSnout;
  G4VPhysicalVolume*   physiSnout;

  //add more
  //G4LogicalVolume*     logicSupport;
  //G4VPhysicalVolume*   physiSupport;


private:
  void CreateSolids();
  void MakeMaterials();

};

#endif
