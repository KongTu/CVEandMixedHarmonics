import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

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
process.GlobalTag.globaltag = '75X_dataRun2_PromptHI_v3'

process.ana = EDAnalyzer("V0AnalyzerSimpleNtuple",

                trackSrc = cms.InputTag("generalTracks"),
                vertexSrc = cms.InputTag("offlinePrimaryVertices"),
                genParticleSrc = cms.InputTag("genParticles"),
                generalV0_ks = cms.InputTag("generalV0CandidatesNew:Kshort"),
                generalV0_la = cms.InputTag("generalV0CandidatesNew:Lambda"),
                doGenParticle = cms.bool(False),

)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'/store/user/davidlw/HIMinimumBias5/RecoSkim2015_pprereco_v5/160727_021048/0000/PbPb_MB_1.root'
)
)

process.TFileService = cms.Service("TFileService",fileName = cms.string("v0Ntuple.root"))
process.p = cms.Path(  process.ana )