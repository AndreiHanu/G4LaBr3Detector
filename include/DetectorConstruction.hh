#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4GenericMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
    public:
  	    // Constructor
        DetectorConstruction();
        // Destructor
        virtual ~DetectorConstruction();

        // Defines the detector geometry and returns a pointer to the physical World Volume
        virtual G4VPhysicalVolume* Construct();
    
        // Sensitive Detector 
	    virtual void ConstructSDandField();
    
        // Set Methods
        void SetDetectorAngle(G4double val);
        void SetSourceRadius(G4double val);
    
        // Get Methods
        G4double GetDetectorAngle();
        G4double GetSourceRadius();
    
    private:
        // Defines all the detector materials
        void DefineMaterials();
    
        // Define commands to change the geometry
        void DefineCommands();
    
        G4GenericMessenger* fMessenger;
        G4bool fCheckOverlaps;
    
        // Standard Materials
        G4Material* fMatWorld;
        G4Material* fMatHousing;
        G4Material* fMatInnerHousing;
        G4Material* fMatSiChip;
        G4Material* fMatElastomer;
        G4Material* fMatLBI;
        G4Material* fMatRC;
        G4Material* fMatRI;
        G4Material* fMatBNCIns;
    
        // Logical Volumes
        G4LogicalVolume* WorldLogical;
        G4LogicalVolume* SourceLogical;
        G4LogicalVolume* HousingLogical;
        G4LogicalVolume* InnerHousingLogical;
        G4LogicalVolume* ElastomerRingLogical;
        G4LogicalVolume* LBILogical;
        G4LogicalVolume* SiChipLogical;
        G4LogicalVolume* ElastLogical;
        G4LogicalVolume* RCLogical;
        G4LogicalVolume* RILogical; 
        G4LogicalVolume* BNCInsLogical;
    
        // Physical Volumes
        G4VPhysicalVolume* WorldPhysical;
        G4VPhysicalVolume* SourcePhysical;
        G4VPhysicalVolume* HousingPhysical;
        G4VPhysicalVolume* InnerHousingPhysical;
        G4VPhysicalVolume* ElastomerRingPhysical;
        G4VPhysicalVolume* LBIPhysical;
        G4VPhysicalVolume* SiChipPhysical;
        G4VPhysicalVolume* ElastPhysical;
        G4VPhysicalVolume* RCPhysical;
        G4VPhysicalVolume* RIPhysical;
        G4VPhysicalVolume* BNCInsPhysical;
    
        // Geometry Parameters
        G4double HOUSING_OD;
        G4double HOUSING_H;
        G4double HOUSING_FRONT_T;
        G4double HOUSING_FRONT_OD;
        G4double HOUSING_BACK_T;
        G4double HOUSING_SIDE_T;
        G4double ELAST_RING_ID;
        G4double ELAST_RING_OD;
        G4double ELAST_RING_H;
        G4double LBI_T;
        G4double LBI_H;
        G4double SiChip_OD;
        G4double SiChip_H;
        G4double ELAST_OD;
        G4double ELAST_H;
        G4double RC_OD;
        G4double RC_H;
        G4double RI_OD;
        G4double RI_ID;
        G4double RI_H;
        G4double BNC_OD;
        G4double BNC_H;
        G4double BNC_INS_OD;

	    // Rotation Angles
	    G4double rotX;
        G4double sourceRadius;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

