import FWCore.ParameterSet.Config as cms

from HeavyIonsAnalysis.PhotonAnalysis.ElectronVID_cff import *

ggHiNtuplizer = cms.EDAnalyzer("ggHiNtuplizer",
    doGenParticles     = cms.bool(True),
    doSuperClusters    = cms.bool(False),
    doElectrons        = cms.bool(False),
    doPhotons          = cms.bool(True),
    doMuons            = cms.bool(False),
    runOnParticleGun   = cms.bool(False),
    useValMapIso       = cms.bool(True),
    doElectronVID      = cms.bool(False),
    doEleERegression   = cms.bool(False),
    doEffectiveAreas   = cms.bool(False),
    doPhoERegression   = cms.bool(False),
    doRecHitsEB        = cms.bool(False),
    doRecHitsEE        = cms.bool(False),
    recHitsEB          = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEB"),
    recHitsEE          = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEE"),
    pileupCollection   = cms.InputTag("addPileupInfo"),
    genParticleSrc     = cms.InputTag("genParticles"),
    superClustersEB    = cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel"),
    superClustersEE    = cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower"),
    gsfElectronLabel   = cms.InputTag("gedGsfElectrons"),
    recoPhotonSrc      = cms.InputTag("islandPhotons"),
    electronVetoID     = electronVetoID25nsV1,
    electronLooseID    = electronLooseID25nsV1,
    electronMediumID   = electronMediumID25nsV1,
    electronTightID    = electronTightID25nsV1,
    recoPhotonHiIsolationMap = cms.InputTag("photonIsolationHIProducerppIsland"),
    recoMuonSrc        = cms.InputTag("muons"),
    VtxLabel           = cms.InputTag("offlinePrimaryVerticesRecovery"),
    rho                = cms.InputTag("fixedGridRhoFastjetAll"),
    beamSpot           = cms.InputTag('offlineBeamSpot'),
    conversions        = cms.InputTag('allConversions'),
    effAreasConfigFile = cms.FileInPath('HeavyIonsAnalysis/PhotonAnalysis/data/EffectiveAreas_94X_v0'),
    doPfIso            = cms.bool(False),
    calcIDTrkIso       = cms.bool(False),
    trackSrc = cms.InputTag("generalTracks"),
    mvaSrc = cms.InputTag("generalTracks","MVAValues"),
    collSystemTag = cms.string("pbpb18"),
    particleFlowCollection = cms.InputTag("particleFlow"),
    removePhotonPfIsoFootprint = cms.bool(False),
    particleBasedIsolationPhoton = cms.InputTag("DUMMY"),
    saveAssociatedPFcands = cms.bool(False),
    doEvtPlane = cms.bool(False),
    evtPlane = cms.InputTag("hiEvtPlane"),
    indexEvtPlane = cms.int32(8),
    calcEvtPlanePF = cms.bool(False),
)

ggHiNtuplizerGED = ggHiNtuplizer.clone(
    doElectrons              = True,
    doMuons                  = True,
    recoPhotonSrc            = 'gedPhotons',
    recoPhotonHiIsolationMap = 'photonIsolationHIProducerppGED',
    doPfIso                  = True,
    removePhotonPfIsoFootprint = True,
    calcIDTrkIso             = True,
    particleBasedIsolationPhoton = cms.InputTag("particleBasedIsolation", "gedPhotons"),
    saveAssociatedPFcands = True,
    doEvtPlane = False,
)
