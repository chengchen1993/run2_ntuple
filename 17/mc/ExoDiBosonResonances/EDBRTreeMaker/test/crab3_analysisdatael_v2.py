from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = 'singleEl-16B-N1'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName='Analysis'
config.JobType.sendExternalFolder=True#  = 'Analysis'
config.JobType.inputFiles = ['L1PrefiringMaps_new.root','Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK4PFchs.txt','Summer16_23Sep2016BCDV3_DATA_L2Relative_AK4PFchs.txt','Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK4PFchs.txt','Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK4PFchs.txt','Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK8PFchs.txt','Summer16_23Sep2016BCDV3_DATA_L2Relative_AK8PFchs.txt','Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK8PFchs.txt','Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK8PFchs.txt','Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK8PFPuppi.txt','Summer16_23Sep2016BCDV3_DATA_L2Relative_AK8PFPuppi.txt','Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK8PFPuppi.txt','Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK8PFPuppi.txt']

# Name of the CMSSW configuration file
#config.JobType.psetName    = 'bkg_ana.py'
config.JobType.psetName    = 'analysis.py'
#config.JobType.allowUndistributedCMSSW = True
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/SingleElectron/Run2016B-03Feb2017_ver2-v2/MINIAOD'
#config.Data.inputDataset = '/SingleElectron/Run2016B-PromptReco-v2/MINIAOD'#Run2015D-PromptReco-v3/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 50
config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'#'Cert_246908-254879_13TeV_PromptReco_Collisions15_JSON.txt' #'Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON.txt'#Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt'#'lumiSummary_13_07_2015_JSON.txt'#https://twiki.cern.ch/twiki/pub/CMS/ExoDijet13TeV/lumiSummary_13_07_2015_JetHT.json'#https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
config.Data.runRange = ''#'250843-250932' # '193093-194075'
#config.Data.runRange = '251244-251252'#'250843-250932' # '193093-194075'
#config.Data.outLFNDirBase='/store/user/chench/'# = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'singleEl-16B-N1'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T3_US_FNALLPC'


##config.Data.inputDataset = '/WJetsToLNu_13TeV-madgraph-pythia8-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDBS = 'global'
##config.Data.inputDBS = 'phys03'
#config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob =10 
#config.Data.totalUnits = 279
#config.Data.publication = False
#
## This string is used to construct the output dataset name
#config.Data.outputDatasetTag = 'WJets100To200_weight'
#
#config.section_("Site")
## Where the output files will be transmitted to
