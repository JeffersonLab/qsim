
#ifndef qsimPrimaryGeneratorAction_h
#define qsimPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"//abstract base class of the user's mandatory action class form primary vertex/particle generation. Has only one pure virtual method GeneratePrimaries() which is invoked from G4RunManager during theevent loop. Class NOT intended for generating primary vertex/particle by itself.
#include "G4String.hh"

class G4ParticleGun;//shoots a particle of a given type into a given direction with either a given momentum. Position and time of primary particle must be set by the corresponding set methods of G4VPrimaryGenerator base class.
class G4Event;//represents an event. Constructed and/or deleted by G4RunManager.Must have one or more primary vertexes and primary particles associated to thosevertexes
class qsimIO;
class qsimEvent;

class qsimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    qsimPrimaryGeneratorAction();
    ~qsimPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent); //Generates Primaries(Argument is the pointer to the event)
    G4ParticleGun* GetParticleGun();// gets the particle gun
    void SetIO( qsimIO *io ){ fIO = io; }
		bool Thetaspectrum(double);
		bool pspectrum(double );
		//void SourceModeSet();
		void SourceModeSet(G4int );
              

  
	private:
    G4ParticleGun* fParticleGun;
		G4int fSourceMode;
                
   
                


    qsimEvent *fDefaultEvent;
    qsimIO *fIO;

  public:
		G4double fXmin, fXmax, fYmin, fYmax;
    G4double fZ;
    G4double fEmin, fEmax;
    
    G4double fthetaMin, fthetaMax;
    G4double fTheta, fPhi;
 
};

#endif


