
#ifndef qsimDetectorConstruction_h
#define qsimDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4ios.hh"

#include "G4RotationMatrix.hh"

class G4Box;
class G4Trd;
class G4VSolid;
class qsimMaterials;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class qsimDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    qsimDetectorConstruction();
   ~qsimDetectorConstruction();
		//void StandModeSet();
		void DetModeSet(G4int );
		void StandModeSet(G4int );
   public:
    G4VPhysicalVolume* Construct();
    void UpdateGeometry();
    void PrintParameters();
    
    // Set Material Commands for World and Scintillator
    void SetWorldMaterial         (G4String);
    void SetCoupleMaterial        (G4String);
    void SetScintMaterial         (G4String);
    void SetPMTMaterial           (G4String);
    
    void SetPhotonDetPolish     (G4double);
    void SetPhotonDetReflectivity (G4double);
    void SetScintThickness      (G4double);

    
    G4Material* FindMaterial(G4String);
  

  private:
    G4double quartz_x;
    G4double quartz_y;
    G4double quartz_z;
	//G4int fStandMode;
	G4int fDetMode;
	G4int fStandMode;

	G4double quartz_zPos;

    G4double cone_rmin1;
    G4double cone_rmax1;
    G4double cone_rmin2;
    G4double cone_rmax2;
    G4double cone_z;
    G4double cone_sphi;
    G4double cone_fphi;
	
	G4double rin;
	G4double rout;
	G4double lngth;
    //defining parameters for standalone scintillator design
    qsimMaterials* materials;
    G4Material* detMaterial;
    G4Material* pmtMaterial;
    G4Material* wrapMaterial;
    
    G4Box*              solidWorld;
    G4LogicalVolume*    logicWorld;
    G4VPhysicalVolume*  physiWorld;

    G4Trd*              solidLightGuide1;
    G4LogicalVolume*    logicLightGuide1;
    G4VPhysicalVolume*  physiLightGuide1;  

    G4Trd*              solidLightGuide2;
    G4LogicalVolume*    logicLightGuide2;
    G4VPhysicalVolume*  physiLightGuide2;  

    G4Trd*              solidLightGuide1_sub;
    G4VSolid*           solidLightGuide1_subtract;

    G4Trd*              solidLightGuide2_sub;
    G4VSolid*           solidLightGuide2_subtract;
    
    G4Box*              solidScintillator;
    G4LogicalVolume*    logicScintillator;
    G4VPhysicalVolume*  physiScintillator;
    
    G4Box*              solidPhotonDet;
    G4LogicalVolume*    logicPhotonDet;
    G4VPhysicalVolume*  physiPhotonDet;

    G4Box*              solidPhotonDet1;
    G4LogicalVolume*    logicPhotonDet1;
    G4VPhysicalVolume*  physiPhotonDet1;

    G4Box*              solidPhotonDet2;
    G4LogicalVolume*    logicPhotonDet2;
    G4VPhysicalVolume*  physiPhotonDet2;

    G4Box*              solidPhotonDet3;
    G4LogicalVolume*    logicPhotonDet3;
    G4VPhysicalVolume*  physiPhotonDet3;

    G4Box*              solidPhotonDet4;
    G4LogicalVolume*    logicPhotonDet4;
    G4VPhysicalVolume*  physiPhotonDet4;

    G4Box*              solidPhotonDet5;
    G4LogicalVolume*    logicPhotonDet5;
    G4VPhysicalVolume*  physiPhotonDet5;
    
    G4double           worldSizeX;
    G4double           worldSizeY;
    G4double           worldSizeZ;
    
    G4double           scintX;
    G4double           scintY;
    G4double           scintZ;
    
    G4double            pmtLength;
    
    G4double pmtPolish;
    G4double pmtReflectivity;

  
  //add two private detector construction routines
  //1. qsimConstruct() and copy current detector construction code here
  //2. cosmicPiConstruct() and implement the cosmicPi scintillator here

  G4VPhysicalVolume* qsimConstruct();
  G4VPhysicalVolume* cosmicPiConstruct();
  void UpdateGeometryParameters();
  public:
	G4double fDetAngle, fQuartzPolish;
	// POSSCAN
	G4double fDetPosX, fDetPosY;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*qsimDetectorConstruction_h*/
