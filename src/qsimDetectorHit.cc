#include "qsimDetectorHit.hh"

G4Allocator<qsimDetectorHit> qsimDetectorHitAllocator;
//defines hit information to be recorded on cathode
//constructor has two arguments, detector ID and the hit ID
qsimDetectorHit::qsimDetectorHit(G4int det, G4int copy){
    fDetID  = det;
    fCopyID = copy;

    f3X = G4ThreeVector(-1e9, -1e9, -1e9);//position in lab frame
    f3P = G4ThreeVector(-1e9, -1e9, -1e9);//momentum in lab frame
    f3V = G4ThreeVector(-1e9, -1e9, -1e9);//origin
    f3D = G4ThreeVector(-1e9, -1e9, -1e9);

    fP  = -1.0;//total momentum
    fE  = -1.0;//total energy
    fM  = -1.0;//total mass

    fTrID  = -1;//G4 track ID
    fPID   = (G4int) 1e9;//particle type ID 
    fmTrID = -1;//mother ID

    fGen   = 1;//process generator type
}
//Destructor
qsimDetectorHit::~qsimDetectorHit(){
}

qsimDetectorHit::qsimDetectorHit(const qsimDetectorHit &right) : G4VHit(){
    // copy constructor

    fDetID  = right.fDetID;
    fCopyID = right.fCopyID;
    f3X     = right.f3X;
    f3P     = right.f3P;
    f3D     = right.f3D;
    f3V     = right.f3V;

    fP      = right.fP;
    fE      = right.fE;
    fM      = right.fM;

    fTrID   = right.fTrID;
    fPID    = right.fPID;
    fmTrID  = right.fmTrID;
    fGen    = right.fGen;
}
//this pointer returns its own address
const qsimDetectorHit& qsimDetectorHit::operator =(const qsimDetectorHit &right){
    (*this) = right;
    return *this;
}

G4int qsimDetectorHit::operator==(const qsimDetectorHit &right ) const {
    return (this==&right) ? 1 : 0;
}
