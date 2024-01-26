import FWCore.ParameterSet.Config as cms
process = cms.Process('FAKE')
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
#process.source = cms.Source("EmptySource", fileNames = cms.untracked.vstring())
process.source = cms.Source("EmptySource", firstTime = cms.untracked.uint64(1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.firstEvent = cms.untracked.PSet(input = cms.untracked.int32(0))
