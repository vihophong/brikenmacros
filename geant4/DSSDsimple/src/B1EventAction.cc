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
// $Id: B1EventAction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
#include "B3Run.hh"
#include "B3DetectorConstruction.hh"
#include "B3PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"

#include "Analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction()
: G4UserEventAction()
{
    const B3DetectorConstruction* detectorConstruction
      = static_cast<const B3DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    nofTube = detectorConstruction->GetnofTube();
		detGroup= detectorConstruction->GetDetGroup();
		tubeGroup= detectorConstruction->GetTubeGroup();
		x=new G4double[nofTube];
		y=new G4double[nofTube];
		z=new G4double[nofTube];

		totalHit=0;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{
  delete[] x;
  delete[] y;
  delete[] z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event* evt)
{     

  if (evt->GetEventID()==0) totalHit=0;
  for (G4int i=0;i<nofTube;i++)
	{
		fTime[i]=0;
		fEdep[i]=0;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void B1EventAction::AddTime(G4int nTube,G4double time)
{
	if (time>fTime[nTube]){
		fTime[nTube]=time;
	}
}
void B1EventAction::AddEdep(G4int nTube,G4double eDeph)
{
	fEdep[nTube]=eDeph+fEdep[nTube];
}
void B1EventAction::AddX(G4int nTube,G4double fX)
{
	x[nTube]=fX;
}
void B1EventAction::AddY(G4int nTube,G4double fY)
{
	y[nTube]=fY;
}
void B1EventAction::AddZ(G4int nTube,G4double fZ)
{
	z[nTube]=fZ;
}

void B1EventAction::EndOfEventAction(const G4Event* event)
{
G4AnalysisManager* mann = G4AnalysisManager::Instance();

  G4PrimaryParticle* primaryParticle = event->GetPrimaryVertex()->GetPrimary();
  G4double ke = primaryParticle->GetKineticEnergy();

  G4int mult=0;
  G4double dssdedep=0;
  for (G4int i=0;i<nofTube;i++)
	{
			if (fEdep[i]>0.)
			{
              //G4cout<<"FIRE!!"<<fEdep[i]<<G4endl;
              /*if(fEdep[i]>180&&fEdep[i]<800) {
			    totalHit++;
                mann->FillH1(1,fEdep[i]);
			  }
              */
              if (i==144){
                  mann->FillH1(1,fEdep[i]);
              }
              if (i==145){
                  mann->FillH1(2,fEdep[i]);
              }
              if (i==146){
                  mann->FillH1(3,fEdep[i]);
              }
              if (i==147){
                  //mann->FillH1(4,fEdep[i]);
              }
              if (i<289) {
                  dssdedep+=fEdep[i];
                  mult++;
              }


			  /*
				mann->FillNtupleDColumn(0,i);
	  		mann->FillNtupleDColumn(1,detGroup[i]);
				mann->FillNtupleDColumn(2,tubeGroup[i]);
				mann->FillNtupleDColumn(3,x[i]);
				mann->FillNtupleDColumn(4,y[i]);
				mann->FillNtupleDColumn(5,z[i]);
				mann->FillNtupleDColumn(6,fEdep[i]);
	  		mann->FillNtupleDColumn(7,fTime[i]);
	  		mann->FillNtupleDColumn(8,ke/MeV);
				mann->AddNtupleRow();
				*/
			}
	}
   mann->FillH1(4,dssdedep);
   mann->FillH1(5,mult);


}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
