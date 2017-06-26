import FWCore.ParameterSet.Config as cms
import os
import sys
import commands

process = cms.Process("PracowniaJ5")

# MessageLogger & co.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

# input files (up to 255 files accepted)
process.source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring( 
    'file:/cms/cms/akalinow/CMS/OverlapTrackFinder/Crab/SingleMuFullEtaTestSample/750_FullEta_v2/data/SingleMu_16_p_1_1_eJz.root', 
                                  ),
skipEvents =  cms.untracked.uint32(0)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10))
 
# PostLS1 geometry used
process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

####Event Setup Producer, OMTF
process.load('L1Trigger.L1TMuonOverlap.fakeOmtfParams_cff')
process.load('L1Trigger.L1TMuonOverlap.simOmtfDigis_cfi')
process.simOmtfDigis.dumpResultToXML=cms.bool(True)

process.cwiczenie= cms.EDAnalyzer("Cwiczenie")

process.MyPath = cms.Path(process.simOmtfDigis + process.cwiczenie)
process.schedule = cms.Schedule(process.MyPath)
