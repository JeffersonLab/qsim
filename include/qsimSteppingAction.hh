
#ifndef __QSIMSTEPPINGACTION_HH
#define __QSIMSTEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class qsimSteppingAction : public G4UserSteppingAction
{//G4UserSteppingAction class represents actions taken place by the user at each end of stepping
  public:
    qsimSteppingAction();
    virtual ~qsimSteppingAction(){};// virtual destructors useful when you can delete an instance of a derived class through a pointer to a base class.

    virtual void UserSteppingAction(const G4Step*);

  private:
    G4bool drawFlag;

  public:
    inline void SetDrawFlag(G4bool val)
    { drawFlag = val; };
};

#endif//__QSIMSTEPPINGACTION_HH
