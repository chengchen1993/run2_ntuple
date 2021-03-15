from WMCore.Configuration import Configuration

config = Configuration()
config.section_("General")
config.General.requestName   = 'QCDPT600to800'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName='Analysis'
config.JobType.sendExternalFolder=True#  = 'Analysis'
config.JobType.inputFiles=['L1PrefiringMaps_new.root','Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L1FastJet_AK8PFchs.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK8PFchs.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK8PFchs.txt','Summer16_07Aug2017_V11_MC_L1FastJet_AK8PFPuppi.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK8PFPuppi.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK8PFPuppi.txt','Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFPuppi.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK4PFPuppi.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFPuppi.txt']
#config.JobType.inputFiles=['L1PrefiringMaps_new.root','PHYS14_25_V2_All_L1FastJet_AK4PFchs.txt','PHYS14_25_V2_All_L2Relative_AK4PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK4PFchs.txt','PHYS14_25_V2_All_L1FastJet_AK8PFchs.txt','PHYS14_25_V2_All_L2Relative_AK8PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK8PFchs.txt']
# Name of the CMSSW configuration file
#config.JobType.psetName    = 'bkg_ana.py'
config.JobType.psetName    = 'analysis.py'
#config.JobType.allowUndistributedCMSSW = True
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
#config.Data.inputDataset = '/WJetsToLNu_13TeV-madgraph-pythia8-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDataset = '/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'Automatic'
config.Data.unitsPerJob =180#5
config.Data.totalUnits = -1
# This string is used to construct the output dataset name
name='WWW'
steam_dir='chench'
config.Data.outLFNDirBase='/store/group/phys_b2g/chench/'#='/store/group/phys_b2g/chench/'#+'/'+name+'/'
config.Data.publication = False
config.Data.outputDatasetTag = 'QCDPT600to800'
config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'
