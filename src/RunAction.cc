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
/// \file B1/src/RunAction.cc
/// \brief Implementation of the B1::RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
// #include "Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <atomic>
#include "ReadOut.hh"
#include "EventAction.hh"

namespace B1
{

	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	RunAction::RunAction()
	{
		// add new units for dose
		//
		const G4double milligray = 1.e-3*gray;
		const G4double microgray = 1.e-6*gray;
		const G4double nanogray  = 1.e-9*gray;
		const G4double picogray  = 1.e-12*gray;

		new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
		new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
		new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
		new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);
		CLHEP::HepRandom::setTheSeed((unsigned)clock());
		// Register accumulable to the accumulable manager
		G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
		accumulableManager->RegisterAccumulable(fEdep);
		accumulableManager->RegisterAccumulable(fEdep2);
		accumulableManager->RegisterAccumulable(N);

		for (int i=0;i<ReadOut::YSIZE; i++)
		{
			for (int j=0;j<ReadOut::XSIZE; j++)
			{
				accumulableManager->RegisterAccumulable(histogram[i][j]);
			}
		}
	}

	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	void RunAction::BeginOfRunAction(const G4Run*)
	{
		/*outFile.open ("edep.txt", std::ios::trunc);
		  if (!outFile){
		  std::cout<<"Ошибка"<<std::endl;
		  }*/
		//N=0;
		// inform the runManager to save random number seed
		G4RunManager::GetRunManager()->SetRandomNumberStore(false);

		// reset accumulables to their initial values
		G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
		accumulableManager->Reset();
		av_X = 0;
		av_Y = 0;
		av_sqrX = 0;
		av_sqrY = 0;
		math_exp_X = 0;
		math_exp_Y = 0;
		varX = 0;
		varY = 0;

	}

	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	void RunAction::EndOfRunAction(const G4Run* run)
	{
		G4int nofEvents = run->GetNumberOfEvent();
		if (nofEvents == 0) return;

		// Merge accumulables
		G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
		accumulableManager->Merge();

		// Compute dose = total energy deposit in a run and its variance
		//
		G4double edep  = fEdep.GetValue();
		G4double edep2 = fEdep2.GetValue();
		G4double nn    =     N.GetValue();

		G4double rms = edep2 - edep*edep/nofEvents;
		if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

		const auto detConstruction = static_cast<const DetectorConstruction*>(
				G4RunManager::GetRunManager()->GetUserDetectorConstruction());
		G4double mass = detConstruction->GetScoringVolume()->GetMass();
		G4double dose = edep/mass;
		G4double rmsDose = rms/mass;

		// Run conditions
		//  note: There is no primary generator action object for "master"
		//        run manager for multi-threaded mode.
		const auto generatorAction = static_cast<const PrimaryGeneratorAction*>(
				G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
		G4String runCondition;
		if (generatorAction)
		{
			const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
			runCondition += particleGun->GetParticleDefinition()->GetParticleName();
			runCondition += " of ";
			G4double particleEnergy = particleGun->GetParticleEnergy();
			runCondition += G4BestUnit(particleEnergy,"Energy");
		}
		// Print
		//
		if (IsMaster()) {
			G4cout << "Total number of events: " << nofEvents << " Total number of interacted: " << N.GetValue() << "   N = " << nn << G4endl
				<< "--------------------End of Global Run-----------------------" << G4endl;
			for (int i=0;i<ReadOut::YSIZE; i++)
			{
				for (int j=0;j<ReadOut::XSIZE; j++)
				{  double hist = static_cast<double>(histogram[i][j].GetValue());
					SumOfEvents += hist;
					av_X += pad_x(j)*hist;
					av_Y += pad_y(i)*hist;
					av_sqrX += pad_x(j)*pad_x(j)*hist;
					av_sqrY += pad_y(i)*pad_y(i)*hist;
				}
				math_exp_X = av_X/SumOfEvents;
				math_exp_Y = av_Y/SumOfEvents;
				varX = av_sqrX/SumOfEvents - math_exp_X*math_exp_X;
				varY = av_sqrY/SumOfEvents - math_exp_Y*math_exp_Y;
				G4cout << "x variance: " <<varX<<G4endl;
				G4cout<<"y variance: "<<varY<<G4endl;

			}
			G4cout<<"BEGIN"<<G4endl;
			for (int i=0;i<ReadOut::YSIZE; i++)
			{
				for (int j=0;j<ReadOut::XSIZE; j++)
				{  
					G4cout<<histogram[i][j].GetValue()<<" ";
				}
				G4cout<<G4endl;

			}
			G4cout<<"END"<<G4endl;

		}
		else {
			G4cout << "Energy deposition " << fEdep.GetValue()/CLHEP::MeV << " MeV "<< G4endl;
			G4cout << "Number of interacted particles " << N.GetValue() << G4endl;
			G4cout << G4endl << "--------------------End of Local Run------------------------" << G4endl;
		}

		

	}

	double RunAction::pad_x(int j)
	{
		return (static_cast<double>(j)+0.5 - 0.5*ReadOut::XSIZE)*(ReadOut::X_SIZE/ReadOut::XSIZE);
	}

	double RunAction::pad_y(int i)
	{
		return (static_cast<double>(i)+0.5 - 0.5*ReadOut::YSIZE)*(ReadOut::Y_SIZE/ReadOut::YSIZE);
	}

	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	void RunAction::AddEdep(G4double edep)
	{
		fEdep  += edep;
		fEdep2 += edep*edep;
	}

	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
