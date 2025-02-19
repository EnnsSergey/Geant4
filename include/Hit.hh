#ifndef HIT_H
#define HIT_H
#include "G4Step.hh"
struct Hit
{
  G4ThreeVector pos; 
  G4double energyDeposit;
  G4TouchableHandle touchable;
   //copyNo = touchable -> GetVolume(0) -> GetCopyNo();
  G4String particle; 
};
#endif
