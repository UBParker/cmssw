import FWCore.ParameterSet.Config as cms

from DQMOffline.Trigger.HLTEGTnPMonitor_cfi import egmGsfElectronIDsForDQM,egHLTDQMOfflineTnPSource,egmPhotonIDSequenceForDQM,egHLTElePhoDQMOfflineTnPSource,photonIDValueMapProducer,egmPhotonIDsForDQM

egammaMonitorHLT = cms.Sequence(
    egmGsfElectronIDsForDQM*
    egHLTDQMOfflineTnPSource*
    egmPhotonIDSequenceForDQM*
    egHLTElePhoDQMOfflineTnPSource
)