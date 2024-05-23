import FWCore.ParameterSet.Config as cms
import os
import sys
import subprocess

process = cms.Process("MojaAnaliza")

# MessageLogger & co.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.MessageLogger.suppressWarning  = cms.untracked.vstring('Geometry','AfterSource','L1T')
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

'''
dataDir='/scratch_cmsse/konec/data/2023D_ParkingDoubleMuonLowMass/'
lsCommand='ls -1 '+dataDir+'|grep root'
print ('command: ',lsCommand)
dir=subprocess.Popen(lsCommand, stdout=subprocess.PIPE,shell=True,text=True)
lsOutput=dir.communicate()[0]
files=[]
for f in lsOutput.split():
  print ('file:'+dataDir+f)
  files.append('file:'+dataDir+f)

print (files)
print (len(files))
'''  

# input files (up to 255 files accepted)
process.source = cms.Source('PoolSource',
  fileNames = cms.untracked.vstring(
    'file:data.root',
#    'root://cms-xrd-global.cern.ch//store/data/Run2023D/ParkingDoubleMuonLowMass0/MINIAOD/22Sep2023_v1-v1/2550000/0419eec5-0ae4-4732-8f06-6d72dd25a149.root',
#   'root://xrootd-cms.infn.it//store/data/Run2023D/ParkingDoubleMuonLowMass0/MINIAOD/22Sep2023_v1-v1/2550000/0419eec5-0ae4-4732-8f06-6d72dd25a149.root',
  ),
)
#process.source.fileNames = files
process.source.skipEvents = cms.untracked.uint32(0)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10))

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data', '')


process.analiza= cms.EDAnalyzer("Analysis",
  muonSrc = cms.InputTag("slimmedMuons"),
  outHist = cms.string("histos.root"),
  debug = cms.bool(True)
)

process.MyPath = cms.Path(process.analiza)
process.schedule = cms.Schedule(process.MyPath)
