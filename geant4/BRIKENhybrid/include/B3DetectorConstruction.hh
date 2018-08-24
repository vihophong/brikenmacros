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
// $Id: B3DetectorConstruction.hh 71323 2013-06-13 16:54:23Z gcosmo $
//
/// \file B3DetectorConstruction.hh
/// \brief Definition of the B3DetectorConstruction class

#ifndef B3DetectorConstruction_h
#define B3DetectorConstruction_h 1

#include "B3DetectorConstruction.hh"

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

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class DetectorMessenger;
/// Detector construction class to define materials and geometry.
///
/// Crystals are positioned in Ring, with an appropriate rotation matrix. 
/// Several copies of Ring are placed in the full detector.

class B3DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B3DetectorConstruction();
    virtual ~B3DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    void setConfig(G4double ,G4double);
		G4int GetnofTube()const{return nofTubes;}
		G4double Geth0()const{return fh0;}
		G4double Getf()const{return ff;}
		G4int* GetDetGroup()const{return detGroup;}
		G4int* GetTubeGroup()const{return tubeGroup;}
		G4double* Getx()const{return x;}
		G4double* Gety()const{return y;}
		G4double* Getz()const{return z;}

        G4double Getxi(G4int i)const{return x[i];}
        G4double Getyi(G4int i)const{return y[i];}
        G4double Getzi(G4int i)const{return z[i];}
  private:
  void readConfiguration(G4String infilename);
  void readConfiguration2(G4String infilename1,G4String infilebname2);
  void readConfiguration3(G4String infilename1,G4String infilebname2,G4String infilename3);
  void rect1(G4String infilename,G4String outfilename, G4double, G4double);
  void rect2(G4String infilename,G4String outfilename,G4double h0);
  void rect3(G4String infilename,G4String outfilename, G4double, G4double);
  void rect4(G4String infilename,G4String outfilename, G4double, G4double);
  void rect5(G4String infilename,G4String outfilename, G4double, G4double);
  void ring(G4String infilename,G4String outfilename,G4double f);
  void ring2(G4String infilename,G4String outfilename);
  void ring3(G4String infilename,G4String outfilename);

  G4double firstRingRadius,firstRingRadiusPre,fh0,ff,secondX0,secondY0,sFactor;

  G4bool fCheckOverlaps;
		G4int nofTubes;
		G4int* detGroup;
		G4int* tubeGroup;
		G4double* x;
		G4double* y;
		G4double* z;
  G4LogicalVolume* DefineTubeGroup(G4int, G4double, G4double, G4double, G4double, G4double, G4Material*, G4Material*);

  DetectorMessenger* fDetectorMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

