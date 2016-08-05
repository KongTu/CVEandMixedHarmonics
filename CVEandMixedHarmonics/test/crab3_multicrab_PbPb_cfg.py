from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

import FWCore.ParameterSet.Config as cms
#load the cfi file and rewrite cross section parameter each time:
process = cms.Process('Demo')
process.load("CVEandMixedHarmonics.CVEandMixedHarmonics.cveandmixedharmonics_cfi")

ntrkRange = [60,80,100,120,140,160,180]

outputName = "multicrab_CVEandMixedHarmonics_PbPb_30_100_v1"

config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.allowUndistributedCMSSW = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cveandmixedharmonics_cfg.py'
config.Data.allowNonValidInputDataset = True
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
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
   
   for num in range(0,6):
		
	print 'double check that centrality range is fram %r to %r' % (ntrkRange[num],ntrkRange[num+1])
      		
	process.ana.Nmin = ntrkRange[num]
	process.ana.Nmax = ntrkRange[num+1]
       	RequestName = outputName + "_" + str(num)
       	DataSetName = '/HIMinimumBias5/HIRun2015-02May2016-v1/AOD'
       	config.General.requestName = RequestName
       	config.Data.inputDataset = DataSetName
       	submit(config)

# python crab3_ppTrackingAnalyzer.py to execute 
# ./multicrab -c status -w crab_projects/ to check status 