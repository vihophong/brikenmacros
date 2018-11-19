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
// $Id: B1SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class
#include <string>
#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B3DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "Analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  nofTube(0)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (nofTube==0) { 
    const B3DetectorConstruction* detectorConstruction
      = static_cast<const B3DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    nofTube = detectorConstruction->GetnofTube();   
  }

  //G4AnalysisManager* mann = G4AnalysisManager::Instance();


  G4String volumeName = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
  G4ThreeVector firePos=step->GetPreStepPoint()->GetPosition();
  if (volumeName.substr(0,4) != "Tube") return;
	G4int volumenum=atoi((volumeName.substr(5,volumeName.length()-4)).c_str());
	G4double edepStep = step->GetTotalEnergyDeposit();
  // collect Global time in this pre step
	if (edepStep>0.){
        G4double timeinStep = step->GetPreStepPoint()->GetGlobalTime();
        fEventAction->AddTime(volumenum,timeinStep/ns);
		fEventAction->AddEdep(volumenum,edepStep/keV);
        fEventAction->AddX(volumenum,firePos.x()/cm);
        fEventAction->AddY(volumenum,firePos.y()/cm);
        fEventAction->AddZ(volumenum,firePos.z()/cm);


    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

