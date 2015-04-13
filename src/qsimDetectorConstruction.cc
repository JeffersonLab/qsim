
#include "qsimDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


#include "qsimDetector.hh"
#include "qsimScintDetector.hh"
#include "G4SDManager.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4GenericTrap.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"

#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

qsimDetectorConstruction::qsimDetectorConstruction()
{

  det_x = det_y = det_z = 275*cm;
  quartz_x = 1.75*cm; 
  quartz_y = 7.*cm; 
  //Change quartz thickness here. 
  quartz_z = 0.5*cm;

  quartz_zPos = -.0*cm;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

qsimDetectorConstruction::~qsimDetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* qsimDetectorConstruction::Construct()
{

//	------------- Materials -------------

  G4double a, z, density;
  G4int nelements;

// Air
// 
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nelements=1);
  Air->AddElement(N, 100.*perCent);
//  Air->AddElement(O, 30.*perCent);

// Quartz
// 
  G4Element* Si = new G4Element("Silicon", "Si", z=14 , a=28*g/mole);

  G4Material* Quartz = new G4Material("Quartz", density= 2.203*g/cm3, nelements=2);
  Quartz->AddElement(Si, 1);
  Quartz->AddElement(O, 2);

// Mirror
// 
  G4Element* Al = new G4Element("Aluminum", "Al", z=13 , a=27*g/mole);

  G4Element* Pb = new G4Element("Lead", "Pb", z=82 , a=207.2*g/mole);

  G4Material* Alu_Mat = new G4Material("Alu_Mat", 2.7*g/cm3, nelements=1);
  Alu_Mat->AddElement(Al, 1);

  G4Material* Pb_Mat = new G4Material("Pb_Mat", 11.34*g/cm3, nelements=1);
  Pb_Mat->AddElement(Pb, 1);

	//G4Material* Pb_Mat=Air; // To remove lead bricks, uncomment.
	
  G4Material* Mirror = new G4Material("Mirror", density= 2.7*g/cm3, nelements=1);
  Mirror->AddElement(Al, 1);


// Let us make cathode from a special metal (reflectivity 0, efficiency of photoelectrons 25%)
  G4Material* CATH = new G4Material("CATH", density= 2.7*g/cm3, nelements=1);
  CATH->AddElement(Al, 1);


//
// ------------ Generate & Add Material Properties Table ------------
//



const G4int nEntries = 190;

	G4double PhotonEnergy[nEntries] =
		{  2.4,2.42,2.44,2.46,2.48,2.5,2.52,2.54,2.56,2.58,
		2.6,2.62,2.64,2.66,2.68,2.7,2.72,2.74,2.76,2.78,
		2.8,2.82,2.84,2.86,2.88,2.9,2.92,2.94,2.96,2.98,
		3,3.02,3.04,3.06,3.08,3.1,3.12,3.14,3.16,3.18,
		3.2,3.22,3.24,3.26,3.28,3.3,3.32,3.34,3.36,3.38,
		3.4,3.42,3.44,3.46,3.48,3.5,3.52,3.54,3.56,3.58,
		3.6,3.62,3.64,3.66,3.68,3.7,3.72,3.74,3.76,3.78,
		3.8,3.82,3.84,3.86,3.88,3.9,3.92,3.94,3.96,3.98,
		4,4.02,4.04,4.06,4.08,4.1,4.12 ,4.14,4.16,4.18,  //Glass cuts off above 4.135eV, 87 entries
		4.2,4.22,4.24,4.26,4.28,4.3,4.32,4.34,4.36,4.38,
		4.4,4.42,4.44,4.46,4.48,4.5,4.52,4.54,4.56,4.58,
		4.6,4.62,4.64,4.66,4.68,4.7,4.72,4.74,4.76,4.78,
		4.8,4.82,4.84,4.86,4.88,4.9,4.92,4.94,4.96,4.98, //  Cut off -> 4.96eV ~ 250nm
		5,5.02,5.04   ,   5.06,5.08,5.1,5.12,5.14,5.16,5.18,   // 5.04eV = 246 nm is the 30% cutoff, 133 entries
		5.2,5.22,5.24,5.26,5.28,5.3,5.32,5.34,5.36,5.38,
		5.4,5.42,5.44,5.46,5.48,5.5,5.52,5.54,5.56,5.58,	
		5.6,5.62,5.64,5.66,5.68,5.7,5.72,5.74,5.76,5.78,
		5.8,5.82,5.84,5.86,5.88,5.9,5.92,5.94,5.96,5.98,
		6,6.02,6.04,6.06,6.08,6.1,6.12,6.14,6.16,6.18   };  // 200 nm

	G4double RefractiveIndex1[nEntries];
	G4double Absorption1[nEntries];
	G4double RefractiveIndex2[nEntries];
	G4double RefractiveIndex3[nEntries];
	G4double Reflectivity4[nEntries];
	G4double Efficiency4[nEntries];
	G4double Reflectivity3[nEntries];

 	for (int i = 0; i < nEntries; i++) {
		RefractiveIndex1[i]= 1.455 -(.005836*PhotonEnergy[i])+(.003374*PhotonEnergy[i]*PhotonEnergy[i]);
		PhotonEnergy[i] = PhotonEnergy[i]*eV;

//Aluminum
//		Reflectivity3[i] = 0; //.6;
	
//Aluminum Real
   if (PhotonEnergy[i] < 4.135*eV) Reflectivity3[i] = .75;  // regularly .75, .7 below  .56/.53/.46 tunes to 50 PEs
		else if (PhotonEnergy[i] >= 4.135*eV && PhotonEnergy[i] < 6.203*eV) Reflectivity3[i] = .7;
   else Reflectivity3[i] = .6;		// .6
		
//ALZAK		
//		if (PhotonEnergy[i] < 3.26*eV) {
//			Reflectivity3[i]=.93; }
//		else { Reflectivity3[i] = 0;}

// No Mirror
//		Reflectivity3[i] = 0;
		
//		Absorption1[i] = 50.*cm;  //Uniform
   
        Absorption1[i] = (exp(4.325)*exp(1.191*PhotonEnergy[i]/eV)*exp(-.213*PhotonEnergy[i]*PhotonEnergy[i]/(eV*eV))*exp(-.04086*PhotonEnergy[i]*PhotonEnergy[i]*PhotonEnergy[i]/(eV*eV*eV)))*m;

       if (Absorption1[i] > 25*m) {Absorption1[i] = 25*m;}

		RefractiveIndex2[i]=1.0+(298.*1e-6);
		RefractiveIndex3[i]=0;
		Reflectivity4[i]=0;

                if( PhotonEnergy[i] > 2.7*eV ){
                    Efficiency4[i]=.25;
                } else if ( 2.7*eV >= PhotonEnergy[i]  ) {
                    Efficiency4[i]=.2;
                } 

	}
	
//QUARTZ
	
  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();
  myMPT1->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1,nEntries);
  myMPT1->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption1,     nEntries);
  
  Quartz->SetMaterialPropertiesTable(myMPT1);

//
// Air
//

/*  G4double RefractiveIndex2[nEntries] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };  */

  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex2, nEntries);
  myMPT2->AddConstProperty("SCINTILLATIONYIELD", 140./MeV);


  const G4int NUMENTRIES = 2;
  G4double Scnt_PP[NUMENTRIES] = { 2.*eV, 6*eV };

  G4double Scnt_FAST[NUMENTRIES] = { 0.5, 0.5 };
  G4double Scnt_SLOW[NUMENTRIES] = { 0.5, 0.5 };

  myMPT2->AddProperty("FASTCOMPONENT", Scnt_PP, Scnt_FAST, NUMENTRIES);
  myMPT2->AddProperty("SLOWCOMPONENT", Scnt_PP, Scnt_SLOW, NUMENTRIES);

  myMPT2->AddConstProperty("RESOLUTIONSCALE", 2.0);
  myMPT2->AddConstProperty("FASTTIMECONSTANT",  1.*ns);
  myMPT2->AddConstProperty("SLOWTIMECONSTANT", 10.*ns);
  myMPT2->AddConstProperty("YIELDRATIO", 1.0);
  
  Air->SetMaterialPropertiesTable(myMPT2);





//
// Mirror (refractive index = 0) 
//


/*  G4double RefractiveIndex3[nEntries] =
            { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00 };   */


 

//  G4MaterialPropertiesTable* myMPT3 = new G4MaterialPropertiesTable();
//  myMPT3->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex3, nEntries);
  // myMPT3->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption3, nEntries);  
  //  myMPT3->AddProperty("REFLECTIVITY", PhotonEnergy, Reflectivity3, nEntries);
  //  myMPT3->AddProperty("EFFICIENCY",    PhotonEnergy, Efficiency3, nEntries);  

//  Mirror->SetMaterialPropertiesTable(myMPT3);

//
// CATH
//

/*
 G4double Reflectivity4[nEntries] =
            { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00 };

  G4double Efficiency4[nEntries] =
           {0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25, 0.25 };  */



 
  //G4MaterialPropertiesTable* myMPT4 = new G4MaterialPropertiesTable();
  //myMPT4->AddProperty("REFLECTIVITY",       PhotonEnergy, Reflectivity4,nEntries);
  //myMPT4->AddProperty("EFFICIENCY",    PhotonEnergy, Efficiency4, nEntries);  

  //CATH->SetMaterialPropertiesTable(myMPT4);

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

//
//	------------- Volumes --------------

// The detector
//
  G4Box* det_box = new G4Box("World",det_x,det_y,det_z);

  G4LogicalVolume* det_log
    = new G4LogicalVolume(det_box,Air,"World",0,0,0);

  det_log->SetVisAttributes(G4VisAttributes::Invisible);


  G4VPhysicalVolume* det_phys
    = new G4PVPlacement(0,G4ThreeVector(),det_log,"World",0,false,0);

  double PMT_rad = 2.54*cm;

  double tubelen = 12*2.54*cm;

  double mirthick = 2.54*cm*0.027;
  double wallthick = 0.5*cm;

  double cutboxsize = 10.0*cm;

  // Mirror ends

  G4Tubs *mirbase =  new G4Tubs("mirbase", 0*cm, PMT_rad+wallthick, 5*cm, 0*deg, 360*deg);
  G4Box  *mircut =  new G4Box("mirsub", cutboxsize, cutboxsize, cutboxsize);

  G4RotationMatrix* rotcut = new G4RotationMatrix;
  rotcut->rotateX(M_PI/4.*rad);

  G4RotationMatrix* rotcut2 = new G4RotationMatrix;
  rotcut2->rotateX(-M_PI/4.*rad);

  G4RotationMatrix* pmtcut = new G4RotationMatrix;
  pmtcut->rotateX(M_PI/2.*rad);
  G4Tubs *pmtcutbase =  new G4Tubs("pmtcutbase", 0*cm, PMT_rad, 5*cm, 0*deg, 360*deg);

  G4RotationMatrix* rotmir = new G4RotationMatrix;
  rotmir->rotateY(M_PI*rad);

  G4SubtractionSolid *cut1 = new G4SubtractionSolid("cut1", mirbase, mircut, rotcut, G4ThreeVector(0.0,(cutboxsize+mirthick/2)/sqrt(2.0), (cutboxsize+mirthick/2)/sqrt(2.0)));
  G4SubtractionSolid *mirror= new G4SubtractionSolid("mirror", cut1, mircut, rotcut, G4ThreeVector(0.0, -(cutboxsize+mirthick/2)/sqrt(2.0), -(cutboxsize+mirthick/2)/sqrt(2.0)));

  G4LogicalVolume* mirror_log
    = new G4LogicalVolume(mirror,Alu_Mat,"mirror_log",0,0,0);

  G4VPhysicalVolume* mir_phys_up
    = new G4PVPlacement(0 ,G4ThreeVector(0.0*mm,0, -tubelen-mirthick*sqrt(2.0)/2),mirror_log,"mirr_phys_up",
                        det_log,false,0);
  G4VPhysicalVolume* mir_phys_down
    = new G4PVPlacement(rotmir,G4ThreeVector(0.0*mm,0, tubelen+mirthick*sqrt(2.0)/2),mirror_log,"mirr_phys_down",
                        det_log,false,0);
  
  // Tube Wall

  G4Tubs *tubebase =  new G4Tubs("tubebase", PMT_rad, PMT_rad+wallthick, tubelen-PMT_rad-(wallthick+mirthick/2)/sqrt(2.0), 0*deg, 360*deg);

  G4SubtractionSolid *tcut1 = new G4SubtractionSolid("tcut1", tubebase, mircut, rotcut, G4ThreeVector(0.0,-(cutboxsize)/sqrt(2.0), (cutboxsize)/sqrt(2.0) + tubelen ));
  G4SubtractionSolid *tcut2= new G4SubtractionSolid("tcut2", tcut1, mircut, rotcut, G4ThreeVector(0.0, -(cutboxsize)/sqrt(2.0), -(cutboxsize)/sqrt(2.0) - tubelen));

  G4SubtractionSolid *tcut3= new G4SubtractionSolid("tcut3", tubebase, pmtcutbase, pmtcut, G4ThreeVector(0.0, 0.0, -tubelen));
  G4SubtractionSolid *tcut4= new G4SubtractionSolid("tcut4", tcut3, pmtcutbase, pmtcut, G4ThreeVector(0.0, 0.0, tubelen));


  G4LogicalVolume* tube_log
    = new G4LogicalVolume(tubebase,Alu_Mat,"tube_log",0,0,0);

   G4VisAttributes * tubeVisAtt
       = new G4VisAttributes(G4Colour(0.4,0.4,0.4));

   tube_log->SetVisAttributes(tubeVisAtt);

  G4VPhysicalVolume* tube_phys
    = new G4PVPlacement(0,G4ThreeVector(0.0*mm,0, 0),tube_log,"tube_phy",
                        det_log,false,0);

  //pmt
  
  G4double plngth = 1.5*mm;    
   G4double clngth = 0.1*mm;
  G4Tubs* pmt = new G4Tubs("PMTvol",0.0,PMT_rad,plngth,0,360*deg);

  G4LogicalVolume* pmt_log
      = new G4LogicalVolume(pmt,Air,"PMT",0,0,0);

  G4Tubs* cath = new G4Tubs("CATH",0,PMT_rad,clngth,0,360*deg);
  G4LogicalVolume* cath_log1
      = new G4LogicalVolume(cath,CATH,"CATHlog1",0,0,0);
  G4LogicalVolume* cath_log2
      = new G4LogicalVolume(cath,CATH,"CATHlog2",0,0,0);

  qsimDetector* cathSD1 = new qsimDetector("cathup", 1);
  qsimDetector* cathSD2 = new qsimDetector("cathdn", 2);

  SDman->AddNewDetector(cathSD1);
  cath_log1->SetSensitiveDetector(cathSD1);

  SDman->AddNewDetector(cathSD2);
  cath_log2->SetSensitiveDetector(cathSD2);

   G4VisAttributes * cathVisAtt
       = new G4VisAttributes(G4Colour(0.7,0.7,1.0));

   cath_log1->SetVisAttributes(cathVisAtt);
   cath_log2->SetVisAttributes(cathVisAtt);

  // Make PMT Sensitive
	
  G4String DetSDname = "pmt_up";
   /*
  
  G4String DetSDname = "pmt_up";

  qsimDetector* trackerSD = new qsimDetector(DetSDname, 1);
  
  SDman->AddNewDetector(trackerSD);
  pmt_log->SetSensitiveDetector(trackerSD);
  */

  G4VPhysicalVolume* pmt_phys;

    pmt_phys = new G4PVPlacement(pmtcut,G4ThreeVector(0.0,PMT_rad+wallthick,-tubelen),  
                        pmt_log,"PMT",
                        det_log,false,0);	

    pmt_phys = new G4PVPlacement(pmtcut,G4ThreeVector(0.0,PMT_rad+wallthick,tubelen),  
                        pmt_log,"PMT",
                        det_log,false,0);	

  G4VPhysicalVolume* cath_phys;

    cath_phys = new G4PVPlacement(pmtcut,G4ThreeVector(0.0,PMT_rad+wallthick+plngth*2+clngth,-tubelen),  
                        cath_log1,"CATH_phys_up",
                        det_log,false,0);	

    cath_phys = new G4PVPlacement(pmtcut,G4ThreeVector(0.0,PMT_rad+wallthick+plngth*2+clngth,tubelen),  
                        cath_log2,"CATH_phys_dn",
                        det_log,false,0);	


  G4OpticalSurface* MOpSurface = new G4OpticalSurface("MirrorOpSurface");
  G4OpticalSurface* CTHOpSurface = new G4OpticalSurface("CathodeOpSurface");
  G4OpticalSurface* TubeOpSurface = new G4OpticalSurface("TubeOpSurface");

  MOpSurface -> SetType(dielectric_metal);
  MOpSurface -> SetFinish(ground);
  MOpSurface -> SetModel(glisur);

  CTHOpSurface -> SetType(dielectric_metal);
  CTHOpSurface -> SetFinish(polishedlumirrorair);
  CTHOpSurface -> SetModel(glisur);

  TubeOpSurface -> SetType(dielectric_metal);
  TubeOpSurface -> SetFinish(ground);
  TubeOpSurface -> SetModel(glisur);

  //  G4double polish = 0.8;

/*  G4double Reflectivity3[nEntries] =
            { 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,
              0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,
              0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,
              0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,
              0.90, 0.90, 0.90, 0.90 };  */

  const G4int num = 2;
  G4double Ephoton[num] = {2.038*eV, 4.144*eV};
  G4double Reflectivity5[num] = {0,0};

  G4MaterialPropertiesTable* MOpSurfaceProperty = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable* COpSurfaceProperty = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable* TubeSurfaceProperty = new G4MaterialPropertiesTable();

  MOpSurfaceProperty -> AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity3,nEntries);
  MOpSurface -> SetMaterialPropertiesTable(MOpSurfaceProperty);

  COpSurfaceProperty -> AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity4,nEntries);
  COpSurfaceProperty -> AddProperty("EFFICIENCY",PhotonEnergy,Efficiency4,nEntries);

  CTHOpSurface -> SetMaterialPropertiesTable(COpSurfaceProperty);

  TubeSurfaceProperty -> AddProperty("REFLECTIVITY",Ephoton,Reflectivity5,2);

  TubeOpSurface -> SetMaterialPropertiesTable(TubeSurfaceProperty);
 
  G4LogicalSkinSurface* FrontSurface_1 = new
      G4LogicalSkinSurface("FrontMirrorOpS_1",mirror_log,MOpSurface);

  G4LogicalSkinSurface* CathSurface = new
      G4LogicalSkinSurface("CathOpS1", cath_log1,CTHOpSurface);
  new
      G4LogicalSkinSurface("CathOpS2", cath_log2,CTHOpSurface);

  G4LogicalSkinSurface* TubeSurface = new
      G4LogicalSkinSurface("TubeOpS", tube_log,TubeOpSurface);



  return det_phys;



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
