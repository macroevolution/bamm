#ifndef ____PreservationRateProposal__



#define ____PreservationRateProposal__
 

#include "Proposal.h"
   

class Random;
class Settings;
class Model;
class Prior;

class PreservationRateProposal : public Proposal
{
    
public:
    
    PreservationRateProposal(Random& random, Settings& settings,
                                Model& model, Prior& prior);
    
    virtual void propose();
    virtual void accept();
    virtual void reject();
    
    virtual double acceptanceRatio();
    
    
private:
 
    double getCurrentParameterValue();
    void setProposedParameterValue();
    void revertToOldParameterValue();
    
    double computeLogPriorRatio();
    double computeLogQRatio();
    
    Random& _random;
    Settings& _settings;
    Model& _model;
    Prior& _prior;
    
    
    double _currentParameterValue;
    double _proposedParameterValue;
    double _currentLogLikelihood;
    double _currentLogPrior;
    
    int _proposedEventCount;
    double _proposedLogLikelihood;
    double _proposedLogPrior;
    
    double _updatePreservationRateScale;
    double _cterm;  

};
 

#endif /* defined(____PreservationRateProposal__) */
