import FWCore.ParameterSet.Config as cms
import sys

ifname='/afs/cern.ch/work/h/halil/STTrigger/crab_workspace/crab_skim_miniaod_ST/results/skim_007E248B-54F2-E311-A601-848F69FD4586_1.root'
ofname='nTuple_ST_multi.root'

if len(sys.argv) > 2: ifname=sys.argv[2]
if len(sys.argv) > 3: ofname=sys.argv[3]

print ifname, '->', ofname


if ifname[0:4] != 'file:': ifname = 'file:' + ifname
if ofname[0:4] != 'file:': ofname = 'file:' + ofname



skipEvents=0
maxEvents=-1

##########################################################

process = cms.Process("nTupler")

# Default values for module parameters
process.load("UserCode.nTupleMaker.nTupleMaker_cfi")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('nTupleMaker')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)

#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( maxEvents ) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/h/halil/STTrigger/skims_st_25ns/skim_1.root',
        'file:/afs/cern.ch/work/h/halil/STTrigger/skims_st_25ns/skim_2.root'
    ),
    skipEvents=cms.untracked.uint32(skipEvents)
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string(ofname)
)

#process.nTupler = cms.EDAnalyzer('nTupleMaker',
#   minTracks = cms.untracked.uint32(1000)
#)

# process.nTupler.minTracks=1000

##########################################################
process.nTupler.addMuons=False
process.nTupler.addElec=True
# process.nTupler.addPhot=True
process.nTupler.addJets=False
process.nTupler.addMET=False
#process.nTupler.labelMuons = "globalMuons"


process.p = cms.Path(process.nTupler)
