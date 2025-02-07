from WMCore.Configuration import Configuration
name = '16'
steam_dir = 'chench'

config = Configuration()
config.section_("General")
config.General.requestName   = 'herwigM3500_R0-5_off'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.inputFiles = ['Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L1FastJet_AK8PFchs.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK8PFchs.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK8PFchs.txt','Summer16_07Aug2017_V11_MC_L1FastJet_AK8PFPuppi.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK8PFPuppi.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK8PFPuppi.txt','Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFPuppi.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK4PFPuppi.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFPuppi.txt']
#config.JobType.inputFiles = ['PHYS14_25_V2_All_L1FastJet_AK4PFchs.txt','PHYS14_25_V2_All_L2Relative_AK4PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK4PFchs.txt','PHYS14_25_V2_All_L1FastJet_AK8PFchs.txt','PHYS14_25_V2_All_L2Relative_AK8PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK8PFchs.txt']
# Name of the CMSSW configuration file
#config.JobType.psetName    = 'bkg_ana.py'
config.JobType.psetName    = 'analysis_sig.py'
#config.JobType.allowUndistributedCMSSW = True
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
#config.Data.inputDataset = '/WJetsToLNu_13TeV-madgraph-pythia8-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDataset = '/WkkToWRadionToWWW_M3500-R0-5_13TeV-madgraph-herwigpp/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
#config.Data.inputDBS = 'global'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob =2
#config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob =450
config.Data.totalUnits = -1
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/phys_b2g/' + steam_dir + '/' + name + '/'
# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'herwigM3500_R0-5_off'


config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'
