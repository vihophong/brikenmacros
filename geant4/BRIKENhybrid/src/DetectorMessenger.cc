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
/// \file hadronic/Hadr00/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
// $Id: DetectorMessenger.cc 77210 2013-11-22 01:58:38Z adotti $
//
//
/////////////////////////////////////////////////////////////////////////
//
// DetectorMessenger
//
// Created: 20.06.08 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

#include "DetectorMessenger.hh"

#include "B3DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(B3DetectorConstruction * Det)
:G4UImessenger(), fDetector(Det)
{
  ftestDir = new G4UIdirectory("/geo/");
  ftestDir->SetGuidance("Geometry parameters");

  fUIcmd = new G4UIcommand("/geo/parms",this);
  fUIcmd->SetGuidance("Geometry parameter h0 f");

  G4UIparameter* symbPrm = new G4UIparameter("h0",'d',false);
  symbPrm->SetGuidance("h0");
  fUIcmd->SetParameter(symbPrm); 

  G4UIparameter* symbPrm1 = new G4UIparameter("f",'d',false);
  symbPrm1->SetGuidance("f");
  fUIcmd->SetParameter(symbPrm1); 

  fUIcmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fUIcmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fUIcmd;
  delete ftestDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fUIcmd ) {
    G4double h0,f;
    std::istringstream is(newValue);
    is >> h0 >> f;
    fDetector->setConfig(h0,f);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

