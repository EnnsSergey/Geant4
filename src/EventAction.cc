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
//
/// \file B1/src/EventAction.cc
/// \brief Implementation of the B1::EventAction class

#include "EventAction.hh"
#include "RunAction.hh"
//#include "G4SystemOfEvents.hh"
#include "Hit.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include <array>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <fmt/core.h>
#include <fmt/format.h>
#include "RunAction.hh"
#include "G4Types.hh"
#include "ReadOut.hh"
float threshold = 0;
//ReadOut readOut;


namespace B1
{


	std::ofstream outFile;

	GEMSensitiveDetector::GEMSensitiveDetector (G4String SDname) : G4VSensitiveDetector(SDname)		
	{
		//outFile.open("hit.txt",std::ios_base::trunc);
		outFile.open("hit.txt",std::ios_base::trunc);

		outFile << fmt::format("{:15} {:>30} {:^32} {:>10}\n", "номер события", "Выделение энергии, keV", "Координаты", "частица");
	};
	GEMSensitiveDetector::~GEMSensitiveDetector() {
		outFile.close();
	}




	G4bool	GEMSensitiveDetector::ProcessHits(G4Step *step, G4TouchableHistory *R0hist){

		::Hit hit;
                auto eventManager = G4EventManager::GetEventManager();
                auto eventAction = static_cast<EventAction*>(eventManager->GetUserEventAction());
		hit.pos = (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition())/2;
		hit.energyDeposit = step->GetTotalEnergyDeposit();
		hit.touchable = step->GetPreStepPoint()->GetTouchableHandle();
		////auto copyNo = touchable -> GetVolume(0) -> GetCopyNo();
		hit.particle = step->GetTrack()->GetParticleDefinition()->GetParticleName();
		auto ionizationEnergy = 30*eV;
		float q = hit.energyDeposit/ionizationEnergy;
		if (hit.particle == "e-"){q = -1.0 * q;}
		G4int evtNo = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
		eventAction->readOut.AddCharge(hit.pos.x(), hit.pos.y(), q);
		/*		for (int i=0;i<ReadOut::YSIZE;i++)
				{
				for(int j=0;j<ReadOut::XSIZE;j++)

		{
			G4cout<<eventAction->readOut.charge[i][j]<< "charge"<<G4endl;
		}
	}
*/
		if (true){  
			outFile << fmt::format("{:>15} {:>30.6f} {:>10.3f} {:>10.3f} {:>10.3f} {:>10} {:>10.3f}\n", evtNo, hit.energyDeposit/CLHEP::keV, hit.pos.x(), hit.pos.y(), hit.pos.z(), hit.particle, step->GetTrack()->GetTotalEnergy()/CLHEP::keV);
			//outFile<<"Выделение энергии: "<<" Координаты хита:  "<<hitPos<<" Частица: "<<std::endl;

		};
		outFile.flush();

		return true;
	}
	//void Initialize(G4HCofThisEvent* HCE);
	//void EndOfEvent(G4HCofThisEvent* HCE);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
	: fRunAction(runAction)
{
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::BeginOfEventAction(const G4Event*)
{
	fEdep = 0.;
	for (int i=0;i<ReadOut::YSIZE;i++)
	{
		for(int j=0;j<ReadOut::XSIZE;j++)

		{
			readOut.charge[i][j]=0;
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
	auto runManager = G4RunManager::GetRunManager();
	auto runAction = static_cast<const RunAction*>(runManager->GetUserRunAction());

	// accumulate statistics in run action
	for (int i=0;i<ReadOut::YSIZE;i++)
	{
		for(int j=0;j<ReadOut::XSIZE;j++)

		{
			if (abs(readOut.charge[i][j])>0){runAction->histogram[i][j]+=1;}


		}
	}
	fRunAction->AddEdep(fEdep);
	if (fEdep>0){
		fRunAction->N++;
		std::cout << "fRunAction->N.GetValue() = " << fRunAction->N.GetValue() << std::endl;
	}
	/*fRunAction->outFile<<fEdep/CLHEP::keV<<" "<< fRunAction->N << "  " << event->GetEventID() <<std::endl;
	  fRunAction->outFile.flush();*/
	G4cout
		<< G4endl
		<< " N=" << fRunAction->N.GetValue() << " "
		<< event->GetEventID() 
		<< G4endl;
}


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

