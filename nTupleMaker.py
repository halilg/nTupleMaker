import FWCore.ParameterSet.Config as cms
import sys

ifname='file:skim_pat.root'
ofname='nTuple.root'

if len(sys.argv) > 2: ifname=sys.argv[2]
if len(sys.argv) > 3: ofname=sys.argv[3]

print ifname, '->', ofname


# if ifname[0:4] != 'file:': ifname = 'file:' + ifname
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
        ifname
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
process.nTupler.addMuons=True
process.nTupler.addElec=True
# process.nTupler.addPhot=True
process.nTupler.addJets=True
process.nTupler.addMET=True
#process.nTupler.labelMuons = "globalMuons"
process.nTupler.bTagD="pfCombinedInclusiveSecondaryVertexV2BJetTags"

process.nTupler.dumpHLT=cms.untracked.vstring(
        "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2",
        "HLT_DoubleMu4_3_Bs_v2",
        "HLT_DoublePhoton60_v1",
        "HLT_DoublePhoton85_v2",
        "HLT_Ele22_eta2p1_WPLoose_Gsf_v3",
        "HLT_Ele22_eta2p1_WPTight_Gsf_v3",
        "HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v3",
        "HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v2",
        "HLT_Ele23_WPLoose_Gsf_v3",
        "HLT_Ele23_WPLoose_Gsf_TriCentralPFJet50_40_30_v2",
        "HLT_Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV_p063_v1",
        "HLT_Ele27_WPLoose_Gsf_v1",
        "HLT_Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV_p063_v1",
        "HLT_Ele27_WPLoose_Gsf_TriCentralPFJet50_40_30_v1",
        "HLT_Ele27_eta2p1_WPLoose_Gsf_v2",
        "HLT_Ele27_eta2p1_WPTight_Gsf_v2",
        "HLT_Ele32_eta2p1_WPTight_Gsf_v2",
        "HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_v1",
        "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v3",
        "HLT_Ele105_CaloIdVT_GsfTrkIdT_v3",
        "HLT_HT200_v1",
        "HLT_HT275_v1",
        "HLT_HT325_v1",
        "HLT_HT425_v1",
        "HLT_HT575_v1",
        "HLT_HT450to470_v1",
        "HLT_HT470to500_v1",
        "HLT_HT500to550_v1",
        "HLT_HT550to650_v1",
        "HLT_HT650_v2",
        "HLT_Mu16_eta2p1_MET30_v1",
        "HLT_IsoMu16_eta2p1_MET30_v1",
        "HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1_v1",
        "HLT_IsoMu17_eta2p1_v3",
        "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v4",
        "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v4",
        "HLT_IsoMu17_eta2p1_MediumIsoPFTau35_Trk1_eta2p1_Reg_v3",
        "HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v5",
        "HLT_DoubleIsoMu17_eta2p1_v3",
        "HLT_IsoMu18_v2",
        "HLT_IsoMu18_CentralPFJet30_BTagCSV_p063_v1",
        "HLT_IsoMu18_TriCentralPFJet50_40_30_v2",
        "HLT_IsoMu20_v3",
        "HLT_IsoMu20_eta2p1_LooseIsoPFTau20_v3",
        "HLT_IsoMu22_v2",
        "HLT_IsoMu22_CentralPFJet30_BTagCSV_p063_v1",
        "HLT_IsoMu22_TriCentralPFJet50_40_30_v2",
        "HLT_IsoMu27_v3",
        "HLT_IsoTkMu18_v2",
        "HLT_IsoTkMu20_v4",
        "HLT_IsoTkMu22_v2",
        "HLT_IsoTkMu27_v3",
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_v3",
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v1",
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v1",
        "HLT_Mu17_Mu8_v1",
        "HLT_Mu17_Mu8_DZ_v2",
        "HLT_Mu17_Mu8_SameSign_DZ_v1",
        "HLT_Mu20_Mu10_v1",
        "HLT_Mu20_Mu10_DZ_v1",
        "HLT_Mu20_Mu10_SameSign_DZ_v1",
        "HLT_Mu17_TkMu8_DZ_v2",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2",
        "HLT_Mu25_TkMu0_dEta18_Onia_v2",
        "HLT_Mu27_TkMu8_v2",
        "HLT_Mu30_TkMu11_v2",
        "HLT_Mu30_eta2p1_PFJet150_PFJet50_v1",
        "HLT_Mu40_TkMu11_v2",
        "HLT_Mu40_eta2p1_PFJet200_PFJet50_v3",
        "HLT_Mu20_v2",
        "HLT_TkMu20_v2",
        "HLT_Mu24_eta2p1_v2",
        "HLT_TkMu24_eta2p1_v2",
        "HLT_Mu27_v2",
        "HLT_TkMu27_v2",
        "HLT_Mu50_v2",
        "HLT_Mu45_eta2p1_v2",
        "HLT_PFHT350_PFMET100_v1",
        "HLT_PFHT550_4JetPt50_v1",
        "HLT_PFHT650_4JetPt50_v1",
        "HLT_PFHT750_4JetPt50_v3",
        "HLT_PFHT600_v3",
        "HLT_PFHT650_v3",
        "HLT_PFHT800_v2",
        "HLT_PFJet40_v4",
        "HLT_PFJet60_v4",
        "HLT_PFJet80_v4",
        "HLT_PFJet140_v4",
        "HLT_PFJet200_v4",
        "HLT_PFJet260_v4",
        "HLT_PFJet320_v4",
        "HLT_PFJet400_v4",
        "HLT_PFJet450_v4",
        "HLT_PFJet500_v4",
        "HLT_DiCentralPFJet55_PFMET110_v1",
        "HLT_PFHT200_v2",
        "HLT_PFHT250_v2",
        "HLT_PFHT300_v2",
        "HLT_PFHT350_v3",
        "HLT_PFHT400_v2",
        "HLT_PFHT475_v2",
        "HLT_MET60_IsoTrk35_Loose_v1",
        "HLT_MET75_IsoTrk50_v2",
        "HLT_MET90_IsoTrk50_v2",
        "HLT_PFMET120_BTagCSV_p067_v1",
        "HLT_PFMET120_Mu5_v1",
        "HLT_PFMET170_NoiseCleaned_v3",
        "HLT_PFMET170_HBHECleaned_v2",
        "HLT_PFMET170_JetIdCleaned_v2",
        "HLT_PFMET170_v2",
        "HLT_PFMET90_PFMHT90_IDTight_v2",
        "HLT_PFMET100_PFMHT100_IDTight_v2",
        "HLT_PFMET110_PFMHT110_IDTight_v2",
        "HLT_PFMET120_PFMHT120_IDTight_v2",
        "HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067_v1",
        "HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_v2",
        "HLT_QuadPFJet_BTagCSV_p037_p11_VBF_Mqq200_v1",
        "HLT_QuadPFJet_BTagCSV_p037_VBF_Mqq460_v1",
        "HLT_QuadPFJet_BTagCSV_p037_p11_VBF_Mqq240_v1",
        "HLT_QuadPFJet_BTagCSV_p037_VBF_Mqq500_v1",
        "HLT_QuadPFJet_VBF_v4",
        "HLT_QuadJet45_TripleBTagCSV_p087_v1",
        "HLT_QuadJet45_DoubleBTagCSV_p087_v1",
        "HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v1",
        "HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_v1",
        "HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160_v1",
        "HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v1",
        "HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172_v1",
        "HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6_v1",
        "HLT_Mu8_TrkIsoVVL_v3",
        "HLT_Mu17_TrkIsoVVL_v2",
        "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v3",
        "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v3",
        "HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v3",
        "HLT_BTagMu_DiJet20_Mu5_v2",
        "HLT_BTagMu_DiJet40_Mu5_v2",
        "HLT_BTagMu_DiJet70_Mu5_v2",
        "HLT_BTagMu_DiJet110_Mu5_v2",
        "HLT_BTagMu_Jet300_Mu5_v2",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3",
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3",
        "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v3",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3",
        "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3",
        "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3",
        "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v3",
        "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v3",
        "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v3",
        "HLT_TripleMu_12_10_5_v2",
        "HLT_Mu3er_PFHT140_PFMET125_v1",
        "HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067_v1",
        "HLT_Mu6_PFHT200_PFMET100_v1",
        "HLT_Mu14er_PFMET100_v1",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3",
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3",
        "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3",
        "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2",
        "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3",
        "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v1",
        "HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_v1",
        "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v3",
        "HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v3",
        "HLT_Dimuon0_Jpsi_Muon_v2",
        "HLT_Dimuon0_Upsilon_Muon_v2",
        "HLT_QuadMuon0_Dimuon0_Jpsi_v2",
        "HLT_QuadMuon0_Dimuon0_Upsilon_v2",
        "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v1",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v1",
        "HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v2",
        "HLT_Photon90_CaloIdL_PFHT500_v3",
        "HLT_DoubleMu8_Mass8_PFHT250_v1",
        "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250_v1",
        "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250_v1",
        "HLT_DoubleMu8_Mass8_PFHT300_v4",
        "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v4",
        "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v4",
        "HLT_Mu10_CentralPFJet30_BTagCSV_p13_v1",
        "HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13_v1",
        "HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400_v1",
        "HLT_Ele15_IsoVVVL_PFHT350_PFMET50_v2",
        "HLT_Ele15_IsoVVVL_PFHT600_v3",
        "HLT_Ele15_IsoVVVL_PFHT350_v2",
        "HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_v1",
        "HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400_v1",
        "HLT_Mu15_IsoVVVL_PFHT350_PFMET50_v2",
        "HLT_Mu15_IsoVVVL_PFHT600_v3",
        "HLT_Mu15_IsoVVVL_PFHT350_v2",
        "HLT_Mu16_TkMu0_dEta18_Onia_v2",
        "HLT_Mu16_TkMu0_dEta18_Phi_v2",
        "HLT_Mu8_v3",
        "HLT_Mu17_v2",
        "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v3",
        "HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v3",
        "HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v3",
        "HLT_Ele33_CaloIdM_TrackIdM_PFJet30_v3",
        "HLT_Ele115_CaloIdVT_GsfTrkIdT_v2",
        "HLT_Mu55_v1",
        "HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v2",
        "HLT_Photon90_CaloIdL_PFHT600_v2",
        "HLT_ECALHT800_v2",
        "HLT_MET200_v1",
        "HLT_Physics_v2",
        "HLT_DoubleMu3_Mass10_v1",
        "HLT_Random_v1",
        "HLT_ZeroBias_v2",
        "HLT_AK4CaloJet30_v3",
        "HLT_AK4CaloJet40_v2",
        "HLT_AK4CaloJet50_v2",
        "HLT_AK4CaloJet80_v2",
        "HLT_AK4CaloJet100_v2",
        "HLT_AK4PFJet30_v3",
        "HLT_AK4PFJet50_v3",
        "HLT_AK4PFJet80_v3",
        "HLT_AK4PFJet100_v3",
        "MC_ReducedIterativeTracking_v1",
        "MC_PFMET_v2",
        "MC_AK4PFJets_v2",
        "MC_PFHT_v2",
        "MC_PFMHT_v2",
        "MC_CaloMET_v1",
        "MC_CaloMET_JetIdCleaned_v1",
        "MC_AK4CaloJets_v1",
        "MC_CaloHT_v1",
        "MC_CaloMHT_v1",
        "MC_AK8PFJets_v2",
        "MC_AK8TrimPFJets_v2",
        "MC_AK8PFHT_v2",
        "MC_AK8CaloHT_v1",
        "MC_Diphoton10_10_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass10_v1",
        "MC_DoubleEle5_CaloIdL_GsfTrkIdVL_MW_v2",
        "MC_Ele5_WPLoose_Gsf_v2",
        "MC_Ele15_Ele10_CaloIdL_TrackIdL_IsoVL_DZ_v2",
        "MC_IsoMu_v4",
        "MC_IsoTkMu15_v3",
        "MC_DoubleMu_TrkIsoVVL_DZ_v1",
        "MC_DoubleGlbTrkMu_TrkIsoVVL_DZ_v1",
        "MC_DoubleMuNoFiltersNoVtx_v1",
        "HLT_Photon500_v1",
        "HLT_Photon600_v1",
        "HLT_Mu300_v1",
        "HLT_Mu350_v1",
        "HLT_MET250_v1",
        "HLT_MET300_v1",
        "HLT_PFMET300_v1",
        "HLT_PFMET400_v1",
        "HLT_HT2000_v1",
        "HLT_HT2500_v1",
        "HLT_IsoTrackHE_v1",
        "HLT_IsoTrackHB_v1"
    )

process.p = cms.Path(process.nTupler)
