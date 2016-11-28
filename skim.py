import FWCore.ParameterSet.Config as cms
import sys

ifname="007E248B-54F2-E311-A601-848F69FD4586.root"
ofname="skim_pat.root"

print sys.argv

if len(sys.argv) > 2: ifname=sys.argv[2]
if len(sys.argv) > 3: ofname=sys.argv[3]

print ifname, '->', ofname

processName = "SKIM"
process = cms.Process(processName)

# this inputs the input files
process.source = cms.Source ("PoolSource",
                        fileNames=cms.untracked.vstring(
                ifname
        )
                )

# tell the process to only run over 100 events (-1 would mean run over
#  everything
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32 (1000)

)

# talk to output module
process.out = cms.OutputModule("PoolOutputModule",
                fileName = cms.untracked.string(ofname),
                outputCommands = cms.untracked.vstring(
                        'drop *',
                        'keep *_slimmedJets_*_*',
                        'keep *_slimmedElectrons_*_*',
                        'keep *_slimmedMuons_*_*',
                        'keep *_slimmedMETs_*_*',
                        'keep *_addPileupInfo_*_*',
                        'keep *_generator_*_*',
                        'keep *_slimmedGenJets_*_*',
                        'keep *_patTrigger_*_*',
                        'keep *_selectedPatTrigger_*_*',
                        'keep *_TriggerResults_*_*',
                )
)

process.outpath = cms.EndPath(process.out)
