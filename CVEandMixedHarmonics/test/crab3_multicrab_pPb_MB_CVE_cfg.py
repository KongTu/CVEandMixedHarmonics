from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

import FWCore.ParameterSet.Config as cms
#load the cfi file and rewrite cross section parameter each time:
process = cms.Process('Demo')
process.load("CVEandMixedHarmonics.CVEandMixedHarmonics.cveandmixedharmonics_cfi")

ntrkRange = [0,35,60,90,120]

hmPaths = ['HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_*_v*']

outputName = "multicrab_CVEandMixedHarmonics_pPb_MB_v1_resubmit"

jsonFile = ['/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/HI/Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt',
'/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/HI/Cert_285952-286009_HI8TeV_PromptReco_Pbp_Collisions16_JSON_NoL1T.txt']

config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.allowUndistributedCMSSW = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cveandmixedharmonics_cfg.py'
config.Data.allowNonValidInputDataset = True
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = jsonFile[0]
config.Data.unitsPerJob = 20
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = outputName

config.Site.storageSite = 'T2_US_MIT'

if __name__ == '__main__':
   from CRABAPI.RawCommand import crabCommand
   from CRABClient.ClientExceptions import ClientException
   from httplib import HTTPException

   config.General.workArea = outputName

   def submit(config):
      try:
           crabCommand('submit', config = config)
      except HTTPException as hte:
           print "Failed submitting task: %s" % (hte.headers)
      except ClientException as cle:
          print "Failed submitting task: %s" % (cle)
 
   print 'This is using the pPb json %r' % jsonFile[0]

   #pPb 8TeV 
   #185-250:
   for num in range(7,8):
	for ntrk in range(2,3):


        	print 'double check that multiplicity range is from %r to %r' % (ntrkRange[ntrk],ntrkRange[ntrk+1])
		process.ana.Nmin = ntrkRange[ntrk]
  		process.ana.Nmax = ntrkRange[ntrk+1]
   		process.hltHM.HLTPaths = [ hmPaths[0] ]
		print 'double check that HM path is %r' % hmPaths[0]
   		RequestName = outputName + "_" + str(ntrkRange[ntrk]) + "_" + str(ntrkRange[ntrk+1]) + "_" + str(num)
   		DataSetName = '/PAMinimumBias' + str(num+1) + '/PARun2016C-PromptReco-v1/AOD'
   		print 'double check that dataset is %r' % DataSetName
		config.General.requestName = RequestName
   		config.Data.inputDataset = DataSetName
   		submit(config)




 


# python crab3_ppTrackingAnalyzer.py to execute 
# ./multicrab -c status -w crab_projects/ to check status 
