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
// $Id: B3RunAction.cc 71323 2013-06-13 16:54:23Z gcosmo $
//
/// \file B3RunAction.cc
/// \brief Implementation of the B3RunAction class
#include <fstream>
#include <stdio.h>
#include "B3DetectorConstruction.hh"
#include "B3RunAction.hh"
#include "B1EventAction.hh"
#include "exGPSPrimaryGeneratorAction.hh"
#include "B3Run.hh"
#include "G4Threading.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Analysis.hh"

using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3RunAction::B3RunAction()
 : G4UserRunAction()
{

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3RunAction::~B3RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* B3RunAction::GenerateRun()
{ return new B3Run; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3RunAction::BeginOfRunAction(const G4Run* run)
{ 
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

  G4AnalysisManager* mann = G4AnalysisManager::Instance();
  mann->SetVerboseLevel(1);
  mann->SetFirstHistoId(1);

  char tempname[500];
  sprintf(tempname,"tempHist%d",run->GetRunID());


  mann->OpenFile(tempname);
  mann->CreateH1("1","NeutronE", 100,0,1000);
  mann->CreateH1("2","NeutronTime", 100,0,1000);

  mann->CreateNtuple("BRIKEN","BRIKEN");
  //! First neutron hits
  mann->CreateNtupleDColumn("id");
  mann->CreateNtupleDColumn("detGroup");
  mann->CreateNtupleDColumn("tubeGroup");
  mann->CreateNtupleDColumn("x");
  mann->CreateNtupleDColumn("y");
  mann->CreateNtupleDColumn("z");
  mann->CreateNtupleDColumn("E");
  mann->CreateNtupleDColumn("T");
  mann->CreateNtupleDColumn("partEnergy");
  //! Second neutron hits
  mann->CreateNtupleDColumn("id2");
  mann->CreateNtupleDColumn("detGroup2");
  mann->CreateNtupleDColumn("tubeGroup2");
  mann->CreateNtupleDColumn("x2");
  mann->CreateNtupleDColumn("y2");
  mann->CreateNtupleDColumn("z2");
  mann->CreateNtupleDColumn("E2");
  mann->CreateNtupleDColumn("T2");
  mann->CreateNtupleDColumn("partEnergy2");

  mann->FinishNtuple();

  // Open an output file
//char temp[50];
//sprintf(temp,"Hist%d",run->GetRunID());
//mann->OpenFile(temp);

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.

  //Get efficiency


    G4int runID=run->GetRunID();

    G4String partName;
    G4double partEnergy;

    //if (G4Threading::G4GetThreadId()==-1){
        ofstream outfile("outGrid.txt",std::ofstream::out | std::ofstream::app);
        const B3DetectorConstruction* detectorConstruction
                = static_cast<const B3DetectorConstruction*>
                (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        G4AnalysisManager* aM = G4AnalysisManager::Instance();

        const exGPSPrimaryGeneratorAction* generatorAction
                = static_cast<const exGPSPrimaryGeneratorAction*>(
                    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
        if (generatorAction)
        {
            partEnergy = generatorAction->GetEbeam();
            partName = generatorAction->GetNameofParticle();
            outfile<<generatorAction->GetPartPos().getX()/cm<<"\t"<<generatorAction->GetPartPos().getY()/cm<<"\t"<<
                     generatorAction->GetPartPos().getZ()/cm<<"\t"<<partEnergy/MeV<<"\t"<<
                     aM->GetH1(1)->entries()<<"\t"<<nofEvents<<endl;
        }
    //}
    //if (runID%10==0) outfile<<"\n"<<detectorConstruction->Geth0()<<" "<<detectorConstruction->Getf()<<" ";
    //  outfile<<(G4double)aM->GetH1(1)->entries()/(G4double)nofEvents*100<<" ";
    //  outfile.close();
    //}
    //G4cout<<"Thread ID is: "<<G4Threading::G4GetThreadId() << G4endl;
    //G4cout<<"FIRE on tubes!!!  --- "<<aM->GetH1(1)->entries()<<"/"<<nofEvents<<G4endl;





  //print
  //
  if (IsMaster())
  {
    G4cout
     << "\n--------------------End of Global Run-----------------------"
     << " \n The run was " << nofEvents << " events "<<"Energy is:"<<partEnergy/MeV<<"MeV";
  }
  else
  {
    G4cout
     << "\n--------------------End of Local Run------------------------"
     << " \n The run was " << nofEvents << " "<< partName<<G4endl;
  }

  G4AnalysisManager* mann = G4AnalysisManager::Instance();
  mann->Write();
  mann->CloseFile();

  //Complete clean up;
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

