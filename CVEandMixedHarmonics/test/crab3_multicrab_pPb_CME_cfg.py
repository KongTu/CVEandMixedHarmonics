from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

import FWCore.ParameterSet.Config as cms
#load the cfi file and rewrite cross section parameter each time:
process = cms.Process('Demo')
process.load("CVEandMixedHarmonics.CVEandMixedHarmonics.cmeandmixedharmonics_cfi")

ntrkRange = [120,150,185,250,300,350,400]

hmPaths = ['HLT_PAFullTracks_Multiplicity120*_v*',
	   'HLT_PAFullTracks_Multiplicity150*_v*',
           'HLT_PAFullTracks_Multiplicity185*_v*',
           'HLT_PAFullTracks_Multiplicity250*_v*']

datasetName = ['/PAHighMultiplicity0/PARun2016C-PromptReco-v1/AOD',
	       '/PAHighMultiplicity1/PARun2016C-PromptReco-v1/AOD',
	       '/PAHighMultiplicity2/PARun2016C-PromptReco-v1/AOD',
	       '/PAHighMultiplicity3/PARun2016C-PromptReco-v1/AOD',
	       '/PAHighMultiplicity4/PARun2016C-PromptReco-v1/AOD',
               '/PAHighMultiplicity5/PARun2016C-PromptReco-v1/AOD',
               '/PAHighMultiplicity6/PARun2016C-PromptReco-v1/AOD',
	       '/PAHighMultiplicity7/PARun2016C-PromptReco-v1/AOD']

outputName = "multicrab_CMEandMixedHarmonics_pPb_HM_v1"

jsonFile = ['/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/HI/Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt',
'/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/HI/Cert_285952-286009_HI8TeV_PromptReco_Pbp_Collisions16_JSON_NoL1T.txt']

config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.allowUndistributedCMSSW = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmeandmixedharmonics_cfg.py'
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

 
   #pPb 8TeV:
   #120-150:
   print 'double check that multiplicity range is from %r to %r' % (ntrkRange[0],ntrkRange[1])
   print 'double check that HM path is %r' % hmPaths[0]
   print 'double check that dataset is %r' % datasetName[0]
   process.ana.Nmin = ntrkRange[0]
   process.ana.Nmax = ntrkRange[1]
   process.hltHM.HLTPaths = [ hmPaths[0] ]
   RequestName = outputName + "_" + str(ntrkRange[0]) + "_" + str(ntrkRange[1])
   DataSetName = datasetName[0]
   config.General.requestName = RequestName
   config.Data.inputDataset = DataSetName
   submit(config)
 
   #150-185:
   print 'double check that multiplicity range is from %r to %r' % (ntrkRange[1],ntrkRange[2])
   print 'double check that HM path is %r and %r' % (hmPaths[0], hmPaths[1])
   print 'double check that dataset is %r' % datasetName[0]
   process.ana.Nmin = ntrkRange[1]
   process.ana.Nmax = ntrkRange[2]
   process.hltHM.HLTPaths = [ hmPaths[0], hmPaths[1] ]
   RequestName = outputName + "_" + str(ntrkRange[1]) + "_" + str(ntrkRange[2])
   DataSetName = datasetName[0]
   config.General.requestName = RequestName
   config.Data.inputDataset = DataSetName
   submit(config)

   #185-250:
   for num in range(1,7):

        print 'double check that multiplicity range is from %r to %r' % (ntrkRange[2],ntrkRange[3])
	process.ana.Nmin = ntrkRange[2]
  	process.ana.Nmax = ntrkRange[3]
   	process.hltHM.HLTPaths = [ hmPaths[2] ]
	print 'double check that HM path is %r' % hmPaths[2]
   	RequestName = outputName + "_" + str(ntrkRange[2]) + "_" + str(ntrkRange[3]) + "_" + str(num)
   	DataSetName = datasetName[num]
   	print 'double check that dataset is %r' % datasetName[num]
	config.General.requestName = RequestName
   	config.Data.inputDataset = DataSetName
   	submit(config)



   #250-400:
   for n in range(0,3):
	print 'double check that multiplicity range is from %r to %r' % (ntrkRange[n+3],ntrkRange[n+4])
        print 'double check that HM path is %r' % hmPaths[3]
        print 'double check that dataset is %r' % datasetName[7]
        process.ana.Nmin = ntrkRange[n+3]
        process.ana.Nmax = ntrkRange[n+4]
        process.hltHM.HLTPaths = [ hmPaths[3] ]
        RequestName = outputName + "_" + str(ntrkRange[n+3]) + "_" + str(ntrkRange[n+4])
        DataSetName = datasetName[7]
        config.General.requestName = RequestName
        config.Data.inputDataset = DataSetName
        submit(config)

 


# python crab3_ppTrackingAnalyzer.py to execute 
# ./multicrab -c status -w crab_projects/ to check status 
