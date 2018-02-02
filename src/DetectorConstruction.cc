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
// ********************************************************************
// G4SiDetector.cc
//
// Description: Definition of the Canberra PD450-15-500AM Passivated
// Implanted Planar Silicon (PIPS) detector used by McMaster University
// to perform dosimetry measurements for the lens of the eye.
//
// NOTE1: McMaster is actually using the ORTEC CR-020-450-500 detector
// but the two models are essentially identical in terms of response.
//
// ********************************************************************

#include "DetectorConstruction.hh"
#include <cmath>

// Units and constants
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// Manager classes
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4GeometryManager.hh"
#include "G4SDManager.hh"

// Store classes
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

// Geometry classes
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"

// Primitive geometry types
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"

// Boolean operations on volumes
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

// Regions
#include "G4Region.hh"

// Messenger classes
#include "G4GenericMessenger.hh"

// Scoring Components
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4PSPassageTrackLength.hh"
#include "G4PSSphereSurfaceCurrent.hh"
#include "G4PSIncidentKineticEnergy.hh"
#include "G4SDParticleFilter.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(): G4VUserDetectorConstruction(), fCheckOverlaps(true),
WorldPhysical(0)
{	
	// Geometry Parameters (Default)
	
    // Housing
	HOUSING_OD = 32.0*mm;					// Housing outer diameter
    HOUSING_H = 12.3*mm;                    // Housing height
    HOUSING_FRONT_T = 0.5*mm;               // Housing front thickness
    HOUSING_FRONT_OD = 23.9*mm;             // Housing front opening outer diameter
    HOUSING_BACK_T = 4.0*mm;                // Housing back cover thickness
    HOUSING_SIDE_T = 1.0*mm;                // Housing side thickness

	// Elastomer Ring
	ELAST_RING_ID = HOUSING_FRONT_OD;		// Inner diameter
	ELAST_RING_OD = ELAST_RING_ID + 2*0.5*mm;
	ELAST_RING_H = 0.5*mm;					// Thickness

	// LBI
	LBI_T = 0.3*mm;							// Thickness
	LBI_H = 2.0*mm;							// Height

	// Silicon Chip
	SiChip_OD = 28.*mm;						// Silicon chip outer diameter
	SiChip_H = 500.*um;						// Silicon chip thickness in microns

	// Elastomer on the back of the Silicon chip
	ELAST_OD = ELAST_RING_OD;				
	ELAST_H = 0.5*mm;

	// Brass contact on the back of the elastomer
	RC_OD = 26.*mm;
	RC_H = 1.*mm;

	// RI 
	RI_OD = RC_OD;
	RI_ID = 4.*mm;
	RI_H = 4.*mm;

	// BNC Connector
	BNC_OD = 7.*mm;
	BNC_H = 7.*mm;

	// BNC Insulator (Ceramic)
	BNC_INS_OD = 6.25*mm;

	// Rotation Angle
	rotX = 0.0*deg;		

	// Source Radius
	sourceRadius = 20.*cm;
			 
	// Define Materials
	DefineMaterials();
	
	// Define commands to control the geometry
   	DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    
    // NIST Manager
	G4NistManager* nistManager = G4NistManager::Instance();
	nistManager->SetVerbose(0);

  	// NIST materials
  	//G4Material* galactic = nistManager->FindOrBuildMaterial("G4_Galactic");
	G4Material* Air = nistManager->FindOrBuildMaterial("G4_AIR");
	G4Material* Vacuum = nistManager->FindOrBuildMaterial("G4_Galactic");
	G4Material* Si = nistManager->FindOrBuildMaterial("G4_Si");
    G4Material* Stainless_Steel = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material* Brass = nistManager->FindOrBuildMaterial("G4_BRASS");
	G4Material* Polyethylene = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
	G4Material* Ceramic = nistManager->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
	  
  	// Set the materials for the Geometry
  	//fMatWorld = Air;
	fMatWorld = Vacuum;
    fMatHousing = Stainless_Steel;
	fMatInnerHousing = Air;
    fMatSiChip = Si;
    fMatElastomer = Polyethylene;
    fMatLBI = Polyethylene;
    fMatRC = Brass;
    fMatRI = Polyethylene;
	fMatBNCIns = Ceramic;
  	
  	// Print materials
	//G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 	
	// Cleanup old geometry
  	G4GeometryManager::GetInstance()->OpenGeometry();
  	G4PhysicalVolumeStore::GetInstance()->Clean();
  	G4LogicalVolumeStore::GetInstance()->Clean();
  	G4SolidStore::GetInstance()->Clean();	
  	
	////////////////////////////////////////////////////////////////////////
	// Construct The World Volume

	G4double world_X = 2*(sourceRadius + 1.*cm);
	G4double world_Y = world_X;
	G4double world_Z = world_X;
	
	G4Box* WorldSolid = new G4Box(	"World", world_X/2, world_Y/2, world_Z/2);
  
	WorldLogical = 
		new G4LogicalVolume(WorldSolid,						// The Solid
							fMatWorld,						// Material
							"World");						// Name
  
	WorldPhysical = 
		new G4PVPlacement(	0,								// Rotation
							G4ThreeVector(),				// Translation vector
							WorldLogical,					// Logical volume
							"World",						// Name
							0,								// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// Construct the Source sphere
	// Note: The actual radius of the Source solid will be slightly smaller (0.1 mm) than
	// specified in the macro files in order to allow tracking the incident kinetic energy
	// of particles.
	G4VSolid* SourceSolid = new G4Sphere("SourceSolid", 0., sourceRadius/2, 0., 360.0*degree, 0., 180.0*degree);

	SourceLogical = 
		new G4LogicalVolume(SourceSolid,						// The Solid
							fMatWorld,		    				// Material
							"SourceLogical");	     			// Name

	SourcePhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,0),
							SourceLogical,					// Logical volume
							"SourcePhysical",				// Name
							WorldLogical,					// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check
							
	////////////////////////////////////////////////////////////////////////
	// Construct the detector housing	
	
	G4VSolid* Housing_Cyl = new G4Tubs("Housing_Cyl", 0., HOUSING_OD/2, HOUSING_H/2, 0., 360.*deg);
	G4VSolid* BNC_Cyl = new G4Tubs("BNC_Cyl", 0., BNC_OD/2, BNC_H/2, 0., 360.*deg);

	// Make the inner housing solid through unions
	G4VSolid* Housing_Solid = new G4UnionSolid("HousingSolid", Housing_Cyl, BNC_Cyl, 0, G4ThreeVector(0,0,-HOUSING_H/2 - BNC_H/2));	
	
	HousingLogical = 
		new G4LogicalVolume(Housing_Solid,					// The Solid
							fMatHousing,    				// Material
							"Housing_Logical");     		// Name

    // Rotation matrix for the detector
    G4RotationMatrix Housing_Rotation = G4RotationMatrix();
    Housing_Rotation.rotateX(rotX);

	// Rotation, Translation, and Transformation of the detector			
	G4RotationMatrix Housing_Rot = G4RotationMatrix();
	Housing_Rot.rotateX(rotX);
	
	G4ThreeVector Housing_Trans = G4ThreeVector(0,0,0);

						
	HousingPhysical = 
		new G4PVPlacement(	G4Transform3D(Housing_Rot,Housing_Trans),	// Translation
							HousingLogical,					// Logical volume
							"Housing_Physical",		        // Name
							SourceLogical,					// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check
	
	////////////////////////////////////////////////////////////////////////
	// Construct the inner detector housing (i.e. the air gap)

	G4VSolid* IH_Cyl1 = new G4Tubs("IH_Cyl1", 0., HOUSING_OD/2 - HOUSING_SIDE_T, (HOUSING_H-HOUSING_FRONT_T-HOUSING_BACK_T)/2, 0., 360.*deg);
	G4VSolid* IH_Cyl2 = new G4Tubs("IH_Cyl2", 0., HOUSING_FRONT_OD/2, HOUSING_FRONT_T/2, 0., 360.*deg);
	G4VSolid* IH_Cyl3 = new G4Tubs("IH_Cyl3", 0., BNC_INS_OD/2, (BNC_H + HOUSING_BACK_T)/2 + 0.01*mm, 0., 360.*deg);
	
	// Make the inner housing solid through unions
	G4VSolid* InnerHousingSolid = new G4UnionSolid("PV_Solid", IH_Cyl1, IH_Cyl2, 0, G4ThreeVector(0,0,(HOUSING_H-HOUSING_FRONT_T-HOUSING_BACK_T)/2 + HOUSING_FRONT_T/2));	
	InnerHousingSolid = new G4UnionSolid("PV_Solid", InnerHousingSolid, IH_Cyl3, 0, G4ThreeVector(0,0,-(HOUSING_H-HOUSING_FRONT_T-HOUSING_BACK_T)/2 - (BNC_H + HOUSING_BACK_T)/2 + 0.01*mm));

	InnerHousingLogical = 
		new G4LogicalVolume(InnerHousingSolid,				// The Solid
							fMatInnerHousing,    			// Material
							"InnerHousingLogical");     	// Name

	InnerHousingPhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,(HOUSING_BACK_T-HOUSING_FRONT_T)/2),
							InnerHousingLogical,			// Logical volume
							"InnerHousingPhysical",			// Name
							HousingLogical,					// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check
	
	////////////////////////////////////////////////////////////////////////
	// Construct the elastomer ring (front)
	G4VSolid* ElastomerRingSolid = new G4Tubs("ElastomerRingSolid", ELAST_RING_ID/2, ELAST_RING_OD/2, ELAST_RING_H/2, 0., 360.*deg);

	ElastomerRingLogical = 
		new G4LogicalVolume(ElastomerRingSolid,				// The Solid
							fMatElastomer,    				// Material
							"ElastomerRingLogical");     	// Name

	ElastomerRingPhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,(HOUSING_H-HOUSING_FRONT_T-HOUSING_BACK_T)/2 - ELAST_RING_H/2),
							ElastomerRingLogical,			// Logical volume
							"ElastomerRingPhysical",		// Name
							InnerHousingLogical,			// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check

	
	////////////////////////////////////////////////////////////////////////
	// Construct the LBI polymer
	G4VSolid* LBI_Cyl = new G4Tubs("LBI_Cyl", ELAST_RING_OD/2, HOUSING_OD/2 - HOUSING_SIDE_T, LBI_H/2, 0., 360.*deg);
	G4VSolid* LBI_Cyl_Inner = new G4Tubs("LBI_Cyl_Inner", 0., HOUSING_OD/2 - HOUSING_SIDE_T - LBI_T, (LBI_H - LBI_T)/2, 0., 360.*deg);

	// Make the LBI solid through a subtraction operation
	G4VSolid* LBISolid = new G4SubtractionSolid("LBI_Solid", LBI_Cyl, LBI_Cyl_Inner, 0, G4ThreeVector(0,0,-(LBI_T+0.05*mm)/2));	

	LBILogical = 
		new G4LogicalVolume(LBISolid,						// The Solid
							fMatLBI,    					// Material
							"LBILogical");     				// Name

	LBIPhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,(HOUSING_H-HOUSING_FRONT_T-HOUSING_BACK_T)/2 - LBI_H/2),
							LBILogical,						// Logical volume
							"LBIPhysical",					// Name
							InnerHousingLogical,			// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check
	
	////////////////////////////////////////////////////////////////////////
	// Construct the silicon chip
	G4VSolid* SiChipSolid = new G4Tubs("SiChipSolid", 0., SiChip_OD/2, SiChip_H/2, 0., 360.*deg);

	SiChipLogical = 
		new G4LogicalVolume(SiChipSolid,					// The Solid
							fMatSiChip,    					// Material
							"SiChipLogical");     			// Name

	SiChipPhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,(HOUSING_H-HOUSING_FRONT_T-HOUSING_BACK_T)/2 - ELAST_RING_H - SiChip_H/2),
							SiChipLogical,					// Logical volume
							"SiChipPhysical",				// Name
							InnerHousingLogical,			// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check

	// Create a region for the Si chip so we can apply the PAI model to it
	G4Region* regSiChip = new G4Region("Region_Si_Chip");
  	SiChipLogical->SetRegion(regSiChip);
	regSiChip->AddRootLogicalVolume(SiChipLogical);

	////////////////////////////////////////////////////////////////////////
	// Construct the elastomer pad on the back of the Si chip
	G4VSolid* ElastSolid = new G4Tubs("ElastSolid", 0., ELAST_OD/2, ELAST_H/2, 0., 360.*deg);

	ElastLogical = 
		new G4LogicalVolume(ElastSolid,						// The Solid
							fMatElastomer,    				// Material
							"ElastLogical");     			// Name

	ElastPhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,(HOUSING_H-HOUSING_FRONT_T-HOUSING_BACK_T)/2 - ELAST_RING_H - SiChip_H - ELAST_H/2),
							ElastLogical,					// Logical volume
							"ElastPhysical",				// Name
							InnerHousingLogical,			// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// Construct the brass contact on the back of the elastomer
	G4VSolid* RC_Cyl1 = new G4Tubs("RC_Cyl1", 0., RC_OD/2, RC_H/2, 0., 360.*deg);
	G4VSolid* RC_Cyl2 = new G4Tubs("RC_Cyl2", 0., RC_H/2, (HOUSING_H + BNC_H - HOUSING_FRONT_T - ELAST_RING_H - SiChip_H - ELAST_H - RC_H)/2, 0., 360.*deg);

	G4VSolid* RCSolid = new G4UnionSolid("RCSolid", RC_Cyl1, RC_Cyl2, 0, G4ThreeVector(0,0,-RC_H/2 - (HOUSING_H + BNC_H - HOUSING_FRONT_T - ELAST_RING_H - SiChip_H - ELAST_H - RC_H)/2));

	RCLogical = 
		new G4LogicalVolume(RCSolid,						// The Solid
							fMatRC,		    				// Material
							"RCLogical");	     			// Name

	RCPhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,(HOUSING_H-HOUSING_FRONT_T-HOUSING_BACK_T)/2 - ELAST_RING_H - SiChip_H - ELAST_H - RC_H/2),
							RCLogical,						// Logical volume
							"RCPhysical",					// Name
							InnerHousingLogical,			// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// Construct the RI polymer
	G4VSolid* RISolid = new G4Tubs("RISolid", RI_ID/2, RI_OD/2, RI_H/2, 0., 360.*deg);

	RILogical = 
		new G4LogicalVolume(RISolid,						// The Solid
							fMatRI,		    				// Material
							"RILogical");	     			// Name

	RIPhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,(HOUSING_H-HOUSING_FRONT_T-HOUSING_BACK_T)/2 - ELAST_RING_H - SiChip_H - ELAST_H - RC_H - RI_H/2),
							RILogical,						// Logical volume
							"RIPhysical",					// Name
							InnerHousingLogical,			// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check
	
	////////////////////////////////////////////////////////////////////////
	// Construct the BNC insulator
	G4VSolid* BNCInsSolid = new G4Tubs("BNCInsSolid", RC_H/2, BNC_INS_OD/2, (BNC_H + HOUSING_BACK_T)/2 + 0.01*mm, 0., 360.*deg);

	BNCInsLogical = 
		new G4LogicalVolume(BNCInsSolid,						// The Solid
							fMatBNCIns,		    				// Material
							"BNCInsLogical");	     			// Name

	BNCInsPhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,-(HOUSING_H-HOUSING_FRONT_T-HOUSING_BACK_T)/2 - (BNC_H + HOUSING_BACK_T)/2 + 0.01*mm),
							BNCInsLogical,					// Logical volume
							"BNCInsPhysical",				// Name
							InnerHousingLogical,			// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
  	// Visualisation attributes
  	
  	// World Volume (White)
  	G4VisAttributes* Vis_World = new G4VisAttributes(G4Colour(1.,1.,1.,0.1));
  	Vis_World->SetForceWireframe(true);
  	WorldLogical->SetVisAttributes(Vis_World);

	// Source Volume (Light Yellow)
    G4VisAttributes* Vis_Source = new G4VisAttributes(G4Colour(1.,1.,1.,0.));
    Vis_Source->SetForceWireframe(true);
    SourceLogical->SetVisAttributes(Vis_Source);

    // Housing Volume (Gray)
    G4VisAttributes* Vis_Housing = new G4VisAttributes(G4Colour(0.5,0.5,0.5,.2));
    Vis_Housing->SetForceWireframe(false);
    HousingLogical->SetVisAttributes(Vis_Housing);

	// Inner Housing Volume (Cyan)
    G4VisAttributes* Vis_Inner_Housing = new G4VisAttributes(G4Colour(0.,1.0,1.0,0.2));
    Vis_Inner_Housing->SetForceWireframe(false);
    InnerHousingLogical->SetVisAttributes(Vis_Inner_Housing);

	// Elastomer Ring Volume (Magenta)
    G4VisAttributes* Vis_ElastomerRing = new G4VisAttributes(G4Colour(1.,0.,1.,1.));
    Vis_ElastomerRing->SetForceWireframe(false);
    ElastomerRingLogical->SetVisAttributes(Vis_ElastomerRing);

	// LBI Polymer Volume (Blue)
    G4VisAttributes* Vis_LBI = new G4VisAttributes(G4Colour(0.,0.,1.,0.3));
    Vis_LBI->SetForceWireframe(false);
    LBILogical->SetVisAttributes(Vis_LBI);

	// Silicon Chip Volume (Orange)
    G4VisAttributes* Vis_SiChip = new G4VisAttributes(G4Colour(1.,1.,0.,1.));
    Vis_SiChip->SetForceWireframe(false);
    SiChipLogical->SetVisAttributes(Vis_SiChip);

	// Elastomer Volume (Magenta)
    G4VisAttributes* Vis_Elast = new G4VisAttributes(G4Colour(1.,0.,1.,1.));
    Vis_Elast->SetForceWireframe(false);
    ElastLogical->SetVisAttributes(Vis_Elast);

	// RC Volume (Gold)
    G4VisAttributes* Vis_RC = new G4VisAttributes(G4Colour(1.,.2,.0,1.));
    Vis_RC->SetForceWireframe(false);
    RCLogical->SetVisAttributes(Vis_RC);

	// RI Volume (Dark Gray)
    G4VisAttributes* Vis_RI = new G4VisAttributes(G4Colour(0.1,0.1,0.1,1.));
    Vis_RI->SetForceWireframe(false);
    RILogical->SetVisAttributes(Vis_RI);

	// BNC Insulator Volume (Light Yellow)
    G4VisAttributes* Vis_BNCIns = new G4VisAttributes(G4Colour(0.3,0.3,0.,.5));
    Vis_BNCIns->SetForceWireframe(false);
    BNCInsLogical->SetVisAttributes(Vis_BNCIns);

	////////////////////////////////////////////////////////////////////////
	// Return world volume
	return WorldPhysical; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
	G4String filterName, particleName;
  
  	G4SDParticleFilter* gammaFilter = new G4SDParticleFilter(filterName="gammaFilter",particleName="gamma");
  	G4SDParticleFilter* electronFilter = new G4SDParticleFilter(filterName="electronFilter",particleName="e-");
	
	////////////////////////////////////////////////////////////////////////
	// Construct the Multi Functional Detector for the Si chip
	
	G4MultiFunctionalDetector* SiScorer = new G4MultiFunctionalDetector("Si");
	G4SDManager::GetSDMpointer()->AddNewDetector(SiScorer);	
	G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
	SiChipLogical->SetSensitiveDetector(SiScorer);
 	
 	G4PSEnergyDeposit* eDep = new G4PSEnergyDeposit("eDep");
    SiScorer->RegisterPrimitive(eDep);
	
	G4MultiFunctionalDetector* SourceScorer = new G4MultiFunctionalDetector("Source");
	G4SDManager::GetSDMpointer()->AddNewDetector(SourceScorer);	
	G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
	WorldLogical->SetSensitiveDetector(SourceScorer);

	G4VPrimitiveScorer* kinEGamma = new G4PSIncidentKineticEnergy("kinEGamma");
	kinEGamma->SetFilter(gammaFilter);
    SourceScorer->RegisterPrimitive(kinEGamma);

	G4VPrimitiveScorer* kinEElectron = new G4PSIncidentKineticEnergy("kinEElectron");
	kinEElectron->SetFilter(electronFilter);
    SourceScorer->RegisterPrimitive(kinEElectron);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorAngle(G4double val)
{
	if(WorldPhysical) {    
    	G4Exception ("DetectorConstruction::SetDetectorAngle()", "G4SiDetector", 
                 	JustWarning, 
                 	"Attempt to change already constructed geometry is ignored");
  	} else {
   		rotX = val;
  	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSourceRadius(G4double val)
{
	if(WorldPhysical) {    
    	G4Exception ("DetectorConstruction::SetSourceRadius()", "G4SiDetector", 
                 	JustWarning, 
                 	"Attempt to change already constructed geometry is ignored");
  	} else {
   		sourceRadius = val;
  	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorAngle()
{
	// Return the detector angle
	return rotX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetSourceRadius()
{
	// Return the detector angle
	return sourceRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineCommands()
{
    // Define command directory using generic messenger class
    fMessenger = new G4GenericMessenger(this, "/G4SiDetector/", "Geometry control");

    // Detector Angle Command
    G4GenericMessenger::Command& DetectorAngleCmd
      = fMessenger->DeclareMethodWithUnit("DetectorAngle","deg",
                                  &DetectorConstruction::SetDetectorAngle, 
                                  "Set the angle of the detector within the world volume.");
    DetectorAngleCmd.SetParameterName("angle", true);
    DetectorAngleCmd.SetDefaultValue("0.0");

	// Source Radius Command
	G4GenericMessenger::Command& SourceRadiusCmd
      = fMessenger->DeclareMethodWithUnit("SourceRadius","cm",
                                  &DetectorConstruction::SetSourceRadius, 
                                  "Set the radius of the source volume within the world volume.");
    SourceRadiusCmd.SetParameterName("radius", true);
    SourceRadiusCmd.SetDefaultValue("20.0");
}
