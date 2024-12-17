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
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class
#include <iostream>
#include "DetectorConstruction.hh"
#include "Hit.hh"
#include "ReadOut.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include <fmt/core.h>
#include <fmt/format.h>
#include <vector>

#include "EventAction.hh"
//struct Hit
//{
//	int a;
//
//  //G4ThreeVector pos; 
//  //G4double energyDeposit;
//  //G4TouchableHandle touchable;
//  // //copyNo = touchable -> GetVolume(0) -> GetCopyNo();
//  //G4String particle; 
//};

namespace B1
{
/*
class GEMSensitiveDetector: public G4VSensitiveDetector
{
		std::ofstream outFile;
	public:
		GEMSensitiveDetector (G4String SDname) : G4VSensitiveDetector(SDname)
	        {
		 outFile.open("hit.txt",std::ios_base::trunc);
		 
		 outFile << fmt::format("{:15} {:>30} {:^32} {:>10}\n", "номер события", "Выделение энергии, keV", "Координаты", "частица");
		};
		~GEMSensitiveDetector() override {
		outFile.close();
		}
                
            
		
	public:         
			G4bool ProcessHits(G4Step *step, G4TouchableHistory *R0hist) override {
			::Hit hit;
			hit.pos = (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition())/2;
			hit.energyDeposit = step->GetTotalEnergyDeposit();
			hit.touchable = step->GetPreStepPoint()->GetTouchableHandle();
			////auto copyNo = touchable -> GetVolume(0) -> GetCopyNo();
		        hit.particle = step->GetTrack()->GetParticleDefinition()->GetParticleName();
                        auto ionizationEnergy = 30*eV;
                        float q = hit.energyDeposit/ionizationEnergy;
			G4int evtNo = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
			readOut.AddCharge(hit.pos.x(), hit.pos.y(), q);
			if (true){  
		        	outFile << fmt::format("{:>15} {:>30.6f} {:>10.3f} {:>10.3f} {:>10.3f} {:>10} {:>10.3f}\n", evtNo, hit.energyDeposit/CLHEP::keV, hit.pos.x(), hit.pos.y(), hit.pos.z(), hit.particle, step->GetTrack()->GetTotalEnergy()/CLHEP::keV);
			//outFile<<"Выделение энергии: "<<" Координаты хита:  "<<hitPos<<" Частица: "<<std::endl;
			};
	                outFile.flush();

			return true;
		}
	//void Initialize(G4HCofThisEvent* HCE);
	//void EndOfEvent(G4HCofThisEvent* HCE);
      
};*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructSDandField()
{
 auto Sens = new  GEMSensitiveDetector("Detector");
 auto sdman = G4SDManager::GetSDMpointer();
 sdman->AddNewDetector(Sens);
 logicDet->SetSensitiveDetector(Sens);

}
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //

  G4Material* conv_mat = nist->FindOrBuildMaterial("G4_W");
  G4Material* det_mat = nist->FindOrBuildMaterial("G4_Ar");
  

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  struct GeometryConfig {
	  //converter
	    G4double conv_sizeX = 161.4*mm;
	    G4double conv_sizeY = 50.35*mm;
	    G4double conv_sizeZ = 8*mm;
	   //detector
	    G4double det_sizeX = 161.4*mm;
	    G4double det_sizeY = 50.35*mm;
	    G4double det_sizeZ = 3*mm;
	   //расстояние между детектором и конвертером
	    G4double conv_det_dist = 0*mm;
	  G4double GetZsize (void)const{
		return det_sizeZ + conv_sizeZ + conv_det_dist; 

		  }
	  };
  GeometryConfig geom_config;
  G4double world_sizeX = 1.5*geom_config.conv_sizeX;
  G4double world_sizeY = 1.5*geom_config.conv_sizeY;
  G4double world_sizeZ  = 10*geom_config.GetZsize();
  
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  

  auto solidWorld = new G4Box("World",                           // its name
    0.5 * world_sizeX, 0.5 * world_sizeY, 0.5 * world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
    world_mat,                                       // its material
    "World");                                        // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                           // at (0,0,0)
    logicWorld,                                // its logical volume
    "World",                                   // its name
    nullptr,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    checkOverlaps);                            // overlaps checking

  //
  // Envelope
  //
  auto solidConv = new G4Box("Converter",                    // its name
    0.5 * geom_config.conv_sizeX, 0.5 * geom_config.conv_sizeY, 0.5 * geom_config.conv_sizeZ);  // размеры конвертера
  auto solidDet = new G4Box("Detector", 0.5 * geom_config.det_sizeX, 0.5 * geom_config.det_sizeY, 0.5 * geom_config.det_sizeZ); //размеры детектора 
    

  auto logicConv = new G4LogicalVolume(solidConv,  // its solid
    conv_mat,                                     // its material
    "Converter");                                 // its name
    
  logicDet = new G4LogicalVolume(solidDet, det_mat, "Detector" );
    

  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(0, 0, 0.5*geom_config.conv_sizeZ),     
    logicConv,                 // its logical volume
    "Converter",               // its name
    logicWorld,               // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking
    
new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(0, 0, -0.5*geom_config.det_sizeZ),     
    logicDet,                 // its logical volume
    "Detector",               // its name
    logicWorld,               // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking
	
 

  //
  // Shape 1
  //
  //G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  /*G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Pb");
  G4ThreeVector pos1 = G4ThreeVector(0, 2*cm, -7*cm);

  // Conical section shape
  G4double shape1_rmina =  1.*cm, shape1_rmaxa = 2.*cm;
  G4double shape1_rminb =  2.*cm, shape1_rmaxb = 4.*cm;
  G4double shape1_hz = 3.*cm;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  auto solidShape1 = new G4Cons("Shape1", shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb,
    shape1_hz, shape1_phimin, shape1_phimax);

  auto logicShape1 = new G4LogicalVolume(solidShape1,  // its solid
    shape1_mat,                                        // its material
    "Shape1");                                         // its name

  new G4PVPlacement(nullptr,  // no rotation
    pos1,                     // at position
    logicShape1,              // its logical volume
    "Shape1",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

  //
  // Shape 2
  //
  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  G4ThreeVector pos2 = G4ThreeVector(0, -1*cm, 7*cm);

  // Trapezoid shape
  G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  G4double shape2_dz  = 6*cm;
  auto solidShape2 = new G4Trd("Shape2",  // its name
    0.5 * shape2_dxa, 0.5 * shape2_dxb, 0.5 * shape2_dya, 0.5 * shape2_dyb,
    0.5 * shape2_dz);  // its size

  auto logicShape2 = new G4LogicalVolume(solidShape2,  // its solid
    shape2_mat,                                        // its material
    "Shape2");                                         // its name

  new G4PVPlacement(nullptr,  // no rotation
    pos2,                     // at position
    logicShape2,              // its logical volume
    "Shape2",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

  // Set Shape2 as scoring volume
  //*/
  fScoringVolume = logicConv ;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
