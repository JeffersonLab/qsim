#ifndef __QSIMDETECTORHIT_HH
#define __QSIMDETECTORHIT_HH

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class qsimDetectorHit : public G4VHit {
    public:
	qsimDetectorHit(G4int, G4int);
	~qsimDetectorHit();

	qsimDetectorHit(const qsimDetectorHit &right);
	const qsimDetectorHit& operator=(const qsimDetectorHit &right);
	G4int operator==(const qsimDetectorHit &right) const;
        
        //overload of the operator new operator(mimic malloc function in C). Size_t is the number of bytes to be allocated
	inline void *operator new(size_t);
        //overload of delete operator which returns the pointer of the aHit object
       	inline void operator delete(void *aHit);
	//Can be used to invoke object's constructor or arbritrary piece of memory
        void *operator new(size_t,void*p){return p;}

    private:

    public:
	G4int fDetID;
	G4int fCopyID;

	// Position and momentum in lab coordinates
	G4ThreeVector f3X;
	G4ThreeVector f3P;
	// Total momentum, energy, mass
	G4double fP, fE, fM;
	// Origin
	G4ThreeVector f3V;
	G4ThreeVector f3D;
	// Geant4 track ID, particle type, and mother ID
	G4int    fTrID, fPID, fmTrID;
	// Process generator type
	G4int    fGen;
};

//typedef creates alias name for another data type
//G4THitsCollection template class value type from 
typedef G4THitsCollection<qsimDetectorHit> qsimDetectorHitsCollection;

extern G4Allocator<qsimDetectorHit> qsimDetectorHitAllocator;

//Malloc and Free methods to be used when overloading
inline void* qsimDetectorHit::operator new(size_t){
    void *aHit;
    aHit = (void *) qsimDetectorHitAllocator.MallocSingle();
    return aHit;
}

inline void qsimDetectorHit::operator delete(void *aHit){
    qsimDetectorHitAllocator.FreeSingle( (qsimDetectorHit*) aHit);
}

#endif//__QSIMDETECTORHIT_HH
