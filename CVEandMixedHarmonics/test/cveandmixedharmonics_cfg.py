import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options   = cms.untracked.PSet( wantSummary =
cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("Configuration.StandardSequences.Digi_cff")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v15'

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True),   # otherwise it won't filter the events
)

#Reject beam scraping events standard pp configuration
process.NoScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity7/AOD/PromptReco-v1/000/285/480/00000/02BA31E5-08AF-E611-AAA3-FA163ED00180.root'
#'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAMinimumBias1/AOD/PromptReco-v1/000/285/480/00000/32A34AA3-2CAF-E611-9C0D-FA163E8F093D.root'
'root://cmsxrootd.fnal.gov//store/user/davidlw/PAHighMultiplicity7/RecoSkim2016_Pbp_V0Cascade_v1/170301_191341/0000/pPb_HM_1.root'
)
)

process.load("RecoVertex.V0Producer.generalV0Candidates_cff")
process.generalV0CandidatesNew = process.generalV0Candidates.clone()
process.load("TrackingCode.pileUpFilter.pileUpFilter_cff")
process.load("CVEandMixedHarmonics.CVEandMixedHarmonics.cveandmixedharmonics_cfi")

#define the cuts
process.ana.useCentrality = False
process.ana.doEffCorrection = True
#process.ana.generalV0_ksName = 'generalV0Candidates:Kshort'
#process.ana.generalV0_laName = 'generalV0Candidates:Lambda'
#process.ana.Nmin = 320
#process.ana.Nmax = 1000
process.ana.etaHighHF = 5.0
process.ana.etaLowHF = 3.0
process.ana.n1 = +1
process.ana.n2 = +1
process.ana.n3 = -2

process.TFileService = cms.Service("TFileService",fileName = cms.string("test.root"))
process.p = cms.Path(  process.hfCoincFilter *
                       process.PAprimaryVertexFilter *
                       process.NoScraping *
                       process.hltHM *
		       process.olvFilter_pPb8TeV_dz1p0 *
		       #process.generalV0CandidatesNew *
 		       process.ana)
