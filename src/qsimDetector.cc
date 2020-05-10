#include "qsimDetector.hh"
#include "G4SDManager.hh"

//constructor with two arguments name and detector number
qsimDetector::qsimDetector( G4String name, G4int detnum ) : G4VSensitiveDetector(name){
    
    //array of chars for what??
    char colname[255];
    //detector number 
    fDetNo = detnum;
    assert( fDetNo > 0 );

    fHCID = -1;//hits collection ID??

//    fTrackSecondaries = false;
    fTrackSecondaries = true;
//sprintf is a C function which stores as a C string
    sprintf(colname, "genhit_%s_%d", name.data(), detnum);
    collectionName.insert(G4String(colname));

}

qsimDetector::~qsimDetector(){
}
//Hits collection created by the sensitive detector which is invoked at the beginning
void qsimDetector::Initialize(G4HCofThisEvent *){

    fHitColl = new qsimDetectorHitsCollection( SensitiveDetectorName, collectionName[0] );
}

///////////////////////////////////////////////////////////////////////
//Processes hits, must implement this for generating hits using info of G4Step Object
G4bool qsimDetector::ProcessHits( G4Step *step, G4TouchableHistory *){
    G4bool badedep = false;//bad energy deposition
    G4bool badhit  = false;//bad hit


    // Get touchable volume info (Touchable volume is a geometrical volume which has a unique placement in a detector description)
    G4TouchableHistory *hist = 
	(G4TouchableHistory*)(step->GetPostStepPoint()->GetTouchable());  // Pre->Post
    G4int  copyID = hist->GetReplicaNumber();

    G4StepPoint *poststep = step->GetPostStepPoint();   //Pre->Post
    G4Track     *track   = step->GetTrack();

//    G4Material* material = track->GetMaterial();

//    printf("Standard detector %d hit by %s!\n", fDetNo, track->GetParticleDefinition()->GetParticleName().data());

//    G4double edep = step->GetTotalEnergyDeposit();

    //  Make pointer to new hit if it's a valid track
    qsimDetectorHit *thishit;
    if( !badhit ){
	thishit = new qsimDetectorHit(fDetNo, copyID);
	fHitColl->insert( thishit );
    }

    if( !badhit ){
	// Hit
	thishit->f3X = poststep->GetPosition();
	thishit->f3V = track->GetVertexPosition();
	thishit->f3D = track->GetVertexMomentumDirection();
	thishit->f3P = track->GetMomentum();

	thishit->fP = track->GetMomentum().mag();
	thishit->fE = track->GetTotalEnergy();
	thishit->fM = track->GetDefinition()->GetPDGMass();

	thishit->fTrID  = track->GetTrackID();
	thishit->fmTrID = track->GetParentID();
	thishit->fPID   = track->GetDefinition()->GetPDGEncoding();

	// FIXME - Enumerate encodings
	thishit->fGen   = (long int) track->GetCreatorProcess();
    }

    return !badedep && !badhit;
}

///////////////////////////////////////////////////////////////////////

void qsimDetector::EndOfEvent(G4HCofThisEvent*HCE) {
    G4SDManager *sdman = G4SDManager::GetSDMpointer();

    if(fHCID<0){ fHCID = sdman->GetCollectionID(collectionName[0]); }

    HCE->AddHitsCollection( fHCID, fHitColl );
    return;
}


