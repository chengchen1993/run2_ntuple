import FWCore.ParameterSet.Config as cms

process = cms.Process( "TEST" )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
#,
#				     SkipEvent = cms.untracked.vstring('ProductNotFound'))
filterMode = False # True                
 
######## Sequence settings ##########
corrJetsOnTheFly = True
runOnMC = False
DOHLTFILTERS = True
#useJSON = not (runOnMC)
#JSONfile = 'Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
#****************************************************************************************************#

#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
if runOnMC:
   process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'#'MCRUN2_74_V9::All'
   #process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'#'MCRUN2_74_V9::All'
elif not(runOnMC):
   process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v4'
#sed -i 's/80X_dataRun2_2016SeptRepro_v4/80X_dataRun2_2016SeptRepro_v4/g'  analysis*.py
hltFiltersProcessName = 'RECO'
if runOnMC:
   hltFiltersProcessName = 'PAT' #'RECO'


######### read JSON file for data ##########					                                                             
'''if not(runOnMC) and useJSON:

  import FWCore.PythonUtilities.LumiList as LumiList
  import FWCore.ParameterSet.Types as CfgTypes
  process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
  myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
  process.source.lumisToProcess.extend(myLumis) 
'''
#if DOHLTFILTERS and not(runOnMC):
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
   reverseDecision = cms.bool(False)
)
process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
   reverseDecision = cms.bool(False)
)
######### read JSON file for data ##########                                                                                                 
'''if not(runOnMC) and useJSON:

  import FWCore.PythonUtilities.LumiList as LumiList
  import FWCore.ParameterSet.Types as CfgTypes
  process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
  myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
  process.source.lumisToProcess.extend(myLumis) 
'''

# ---------------------------------------------------------
# DeepAK8: set up TransientTrackBuilder
process.load('Configuration.StandardSequences.MagneticField_cff')
process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName=cms.string('TransientTrackBuilder')
)

####### Redo Jet clustering sequence ##########

from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS, ak8PFJetsCHS, ak8PFJetsCHSPruned, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSPrunedMass, ak8PFJetsCHSSoftDropMass
from CommonTools.PileupAlgos.Puppi_cff import puppi
process.puppi = puppi.clone()
process.puppi.useExistingWeights = True
process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
process.ak8PFJetsCHS = ak8PFJetsCHS.clone( src = 'puppi', jetPtMin = 100.0 )
process.ak8PFJetsCHSPruned = ak8PFJetsCHSPruned.clone( src = 'puppi', jetPtMin = 100.0 )
process.ak8PFJetsCHSPrunedMass = ak8PFJetsCHSPrunedMass.clone()
process.ak8PFJetsCHSSoftDrop = ak8PFJetsCHSSoftDrop.clone( src = 'puppi', jetPtMin = 100.0 )
process.ak8PFJetsCHSSoftDropMass = ak8PFJetsCHSSoftDropMass.clone()
process.NjettinessAK8 = cms.EDProducer("NjettinessAdder",
                   src = cms.InputTag("ak8PFJetsCHS"),
                   Njets = cms.vuint32(1, 2, 3, 4),
                   # variables for measure definition : 
                   measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
                   beta = cms.double(1.0),          # CMS default is 1
                   R0 = cms.double( 0.8 ),          # CMS default is jet cone size
                   Rcutoff = cms.double( 999.0),       # not used by default
                   # variables for axes definition :
                   axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
                   nPass = cms.int32(0),         # not used by default
                   akAxesR0 = cms.double(-999.0)        # not used by default
                   )
process.substructureSequence = cms.Sequence()
process.substructureSequence+=process.puppi
process.substructureSequence+=process.ak8PFJetsCHS
process.substructureSequence+=process.NjettinessAK8
process.substructureSequence+=process.ak8PFJetsCHSPruned
process.substructureSequence+=process.ak8PFJetsCHSPrunedMass
process.substructureSequence+=process.ak8PFJetsCHSSoftDrop
process.substructureSequence+=process.ak8PFJetsCHSSoftDropMass
####### Redo pat jets sequence ##########
process.redoPatJets = cms.Sequence()
process.redoPrunedPatJets = cms.Sequence()
process.redoSoftDropPatJets = cms.Sequence()
from V3Resonances.EDBRJets.redoPatJets_cff import patJetCorrFactorsAK8, patJetsAK8, selectedPatJetsAK8
# Redo pat jets from ak8PFJetsCHS
process.patJetCorrFactorsAK8 = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsCHS' )
process.patJetsAK8 = patJetsAK8.clone( jetSource = 'ak8PFJetsCHS' )
process.patJetsAK8.userData.userFloats.src = [ cms.InputTag("ak8PFJetsCHSPrunedMass"), cms.InputTag("ak8PFJetsCHSSoftDropMass"), cms.InputTag("NjettinessAK8:tau1"), cms.InputTag("NjettinessAK8:tau2"), cms.InputTag("NjettinessAK8:tau3"),cms.InputTag("NjettinessAK8:tau4")]
process.patJetsAK8.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK8") )
process.selectedPatJetsAK8 = selectedPatJetsAK8.clone( cut = cms.string('pt > 100') )
process.redoPatJets+=process.patJetCorrFactorsAK8
process.redoPatJets+=process.patJetsAK8
process.redoPatJets+=process.selectedPatJetsAK8
# Redo pat jets ak8PFJetsCHSPruned
process.patJetCorrFactorsAK8Pruned = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsCHSPruned' )
process.patJetsAK8Pruned = patJetsAK8.clone( jetSource = 'ak8PFJetsCHSPruned' )
process.patJetsAK8Pruned.userData.userFloats.src = [ "" ]
#process.patJetsAK8Pruned.userData.userFloats =cms.PSet(src = cms.VInputTag(""))
process.patJetsAK8Pruned.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK8Pruned") )
process.selectedPatJetsAK8Pruned = selectedPatJetsAK8.clone(cut = 'pt > 100', src = "patJetsAK8Pruned")
process.redoPrunedPatJets+=process.patJetCorrFactorsAK8Pruned
process.redoPrunedPatJets+=process.patJetsAK8Pruned
process.redoPrunedPatJets+=process.selectedPatJetsAK8Pruned
# Redo pat jets ak8PFJetsCHSSoftDrop
process.patJetCorrFactorsAK8Softdrop = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsCHSSoftDrop' )
process.patJetsAK8Softdrop = patJetsAK8.clone( jetSource = 'ak8PFJetsCHSSoftDrop' )
process.patJetsAK8Softdrop.userData.userFloats.src = [ "" ]
#process.patJetsAK8Softdrop.userData.userFloats =cms.PSet(src = cms.VInputTag(""))
process.patJetsAK8Softdrop.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK8Softdrop") )
process.selectedPatJetsAK8Softdrop = selectedPatJetsAK8.clone(cut = 'pt > 100', src = "patJetsAK8Softdrop")
process.redoSoftDropPatJets+=process.patJetCorrFactorsAK8Softdrop
process.redoSoftDropPatJets+=process.patJetsAK8Softdrop
process.redoSoftDropPatJets+=process.selectedPatJetsAK8Softdrop

option = 'RECO'

process.load("V3Resonances.EDBRCommon.goodMuons_cff")
process.load("V3Resonances.EDBRCommon.goodElectrons_cff")
process.load("V3Resonances.EDBRCommon.goodJets_cff")
process.load("V3Resonances.EDBRCommon.leptonicW_cff")
process.load("V3Resonances.EDBRCommon.hadronicW_cff")
process.load("V3Resonances.EDBRCommon.goodPuppi_cff")

if option == 'RECO':
    process.goodMuons.src = "slimmedMuons"
    process.goodElectrons.src = "slimmedElectrons"
    process.goodJets.src = "slimmedJetsAK8"
    process.Wtoenu.MET  = "slimmedMETs"
    process.Wtomunu.MET = "slimmedMETs"
    process.goodPuppi.src = "selectedPatJetsAK8"

process.goodOfflinePrimaryVertex = cms.EDFilter("VertexSelector",
                                       src = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                       cut = cms.string("chi2!=0 && ndof >= 4.0 && abs(z) <= 24.0 && abs(position.Rho) <= 2.0"),
                                       filter = cms.bool(True)
                                       )

if option == 'RECO':
    process.hadronicV.cut = ' '
if option == 'GEN':
    process.hadronicV.cut = ' '
process.hadronicVFilter = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("hadronicV"),
                                       minNumber = cms.uint32(1),
                                       filter = cms.bool(True)
                                       )
process.leptonSequence = cms.Sequence(process.muSequence +
                                      process.eleSequence +
                                      process.leptonicVSequence)
process.jetSequence = cms.Sequence(process.substructureSequence +
                                   process.redoPatJets + 
                                   process.redoPrunedPatJets+
                                   process.redoSoftDropPatJets+
                                   process.fatJetsSequence +
                                   process.fatPuppiSequence+
                                   process.hadronicV +
                                   process.hadronicVFilter) 
if filterMode == False:
    process.goodOfflinePrimaryVertex.filter = False
    process.Wtomunu.cut = ''
    process.Wtoenu.cut = ''
    process.hadronicV.cut = ''
    process.hadronicVFilter.minNumber = 0

######### JEC ########
METS = "slimmedMETs"
jetsAK8 = "slimmedJetsAK8"
jetsAK8pruned = "slimmedJetsAK8"
jetsAK8softdrop = "slimmedJetsAK8"
jetsAK8puppi = "cleanPuppi"
if runOnMC:
   jecLevelsAK8chs = [
                                   'Summer16_23Sep2016V3_MC_L1FastJet_AK8PFchs.txt',
                                   'Summer16_23Sep2016V3_MC_L2Relative_AK8PFchs.txt',
                                   'Summer16_23Sep2016V3_MC_L3Absolute_AK8PFchs.txt'
     ]
   jecLevelsAK8chsGroomed = [
                                   'Summer16_23Sep2016V3_MC_L2Relative_AK8PFchs.txt',
                                   'Summer16_23Sep2016V3_MC_L3Absolute_AK8PFchs.txt'
     ]
   jecLevelsAK8puppi = [
                                   'Summer16_23Sep2016V3_MC_L1FastJet_AK8PFPuppi.txt',
                                   'Summer16_23Sep2016V3_MC_L2Relative_AK8PFPuppi.txt',
                                   'Summer16_23Sep2016V3_MC_L3Absolute_AK8PFPuppi.txt'
     ]
   jecLevelsAK8puppiGroomed = [
                                   'Summer16_23Sep2016V3_MC_L2Relative_AK8PFPuppi.txt',
                                   'Summer16_23Sep2016V3_MC_L3Absolute_AK8PFPuppi.txt'
     ]
   BjecLevelsAK4chs = [
                                   'Summer16_23Sep2016V3_MC_L1FastJet_AK4PFPuppi.txt',
                                   'Summer16_23Sep2016V3_MC_L2Relative_AK4PFPuppi.txt',
                                   'Summer16_23Sep2016V3_MC_L3Absolute_AK4PFPuppi.txt'
     ]
   jecLevelsAK4chs = [
          'Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt',
          'Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt',
          'Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt'
    ]
else:
   jecLevelsAK8chs = [
                                   'Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK8PFchs.txt',
                                   'Summer16_23Sep2016BCDV3_DATA_L2Relative_AK8PFchs.txt',
                                   'Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK8PFchs.txt',
				   'Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK8PFchs.txt'
     ]
   jecLevelsAK8chsGroomed = [
                                   'Summer16_23Sep2016BCDV3_DATA_L2Relative_AK8PFchs.txt',
                                   'Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK8PFchs.txt',
				   'Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK8PFchs.txt'
     ]
   jecLevelsAK8puppi = [
                                   'Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK8PFPuppi.txt',
                                   'Summer16_23Sep2016BCDV3_DATA_L2Relative_AK8PFPuppi.txt',
                                   'Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK8PFPuppi.txt',
                                   'Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK8PFPuppi.txt'
     ]
   jecLevelsAK8puppiGroomed = [
                                   'Summer16_23Sep2016BCDV3_DATA_L2Relative_AK8PFPuppi.txt',
                                   'Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK8PFPuppi.txt',
                                   'Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK8PFPuppi.txt'
     ]
   BjecLevelsAK4chs = [
                                   'Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK4PFPuppi.txt',
                                   'Summer16_23Sep2016BCDV3_DATA_L2Relative_AK4PFPuppi.txt',
                                   'Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK4PFPuppi.txt',
                                   'Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK4PFPuppi.txt'

     ]
   jecLevelsAK4chs = [
                                   'Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK4PFPuppi.txt',
                                   'Summer16_23Sep2016BCDV3_DATA_L2Relative_AK4PFPuppi.txt',
                                   'Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK4PFPuppi.txt',
				   'Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK4PFPuppi.txt'
     ]

process.treeDumper = cms.EDAnalyzer("EDBRTreeMaker",
                                    originalNEvents = cms.int32(1),
                                    crossSectionPb = cms.double(1),
                                    targetLumiInvPb = cms.double(1.0),
                                    isGen = cms.bool(False),
		                isJEC = cms.bool(corrJetsOnTheFly),
			      RunOnMC = cms.bool(runOnMC), 
			      generator =  cms.InputTag("generator"),
                                    genSrc =  cms.InputTag("prunedGenParticles"),
                                    pileup  =   cms.InputTag("slimmedAddPileupInfo"),
		                looseMuonSrc = cms.InputTag("looseMuons"),
                                    looseElectronSrc = cms.InputTag("looseElectrons"),
                                vetoMuonSrc = cms.InputTag("vetoMuons"),
                                    vetoElectronSrc = cms.InputTag("vetoElectrons"),
                                    t1muSrc = cms.InputTag("slimmedMuons"),
                                    metSrc = cms.InputTag("slimmedMETs"),
                                    mets = cms.InputTag(METS),
                                    ak4jetsSrc = cms.InputTag("cleanPuppiAK4"),
                                    ak8JetSrc = cms.InputTag("jetsAK8"),# 
                                    hadronicVSrc = cms.InputTag("hadronicV"),
                                    jets = cms.InputTag("slimmedJets"),
                                    fatjets = cms.InputTag(jetsAK8),
                                    prunedjets = cms.InputTag(jetsAK8pruned),
                                    softdropjets = cms.InputTag(jetsAK8softdrop),
                                    puppijets = cms.InputTag(jetsAK8puppi),
                                    jecAK8chsPayloadNames = cms.vstring( jecLevelsAK8chs ),
			      jecAK8chsPayloadNamesGroomed = cms.vstring( jecLevelsAK8chsGroomed ),
			      jecAK4chsPayloadNames = cms.vstring( jecLevelsAK4chs ),
                                    BjecAK4chsPayloadNames = cms.vstring( BjecLevelsAK4chs ),
			      jecAK8puppiPayloadNames = cms.vstring( jecLevelsAK8puppi ),
                                    jecAK8puppiPayloadNamesGroomed = cms.vstring( jecLevelsAK8puppiGroomed ),
			      jecpath = cms.string(''),
			      rho = cms.InputTag("fixedGridRhoFastjetAll"),
			      vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                    hltToken    = cms.InputTag("TriggerResults","","HLT"),
                                    muPaths1     = cms.vstring("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v*"),
                                    muPaths2     = cms.vstring("HLT_PFHT800_v*"), 
                                    muPaths3     = cms.vstring("HLT_PFHT900_v*"),
                                    muPaths4     = cms.vstring("HLT_PFJet450_v*"), 
                                    muPaths5     = cms.vstring("HLT_PFJet500_v*"), 
                                    muPaths6     = cms.vstring("HLT_AK8PFJet450_v*"),
                                    muPaths7     = cms.vstring("HLT_AK8PFJet500_v*"),
                                    muPaths8     = cms.vstring("HLT_AK8PFJet360_TrimMass30_v*"),
                                    muPaths9     = cms.vstring("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v*"),
                                    noiseFilter = cms.InputTag('TriggerResults','', hltFiltersProcessName),
                                    noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),
                                    noiseFilterSelection_HBHENoiseIsoFilter = cms.string("Flag_HBHENoiseIsoFilter"),
                                    noiseFilterSelection_GlobalTightHaloFilter = cms.string('Flag_globalTightHalo2016Filter'),
                                    noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),
                                    noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),
                                    noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
                                    )


if option=='GEN':
    process.treeDumper.metSrc = 'genMetTrue'
    process.treeDumper.isGen  = True
 

process.analysis = cms.Path(process.leptonSequence +
                            #process.substructureSequence+
                            #process.redoPatJets+
                            #process.redoPrunedPatJets+
                            #process.redoSoftDropPatJets+
                            process.jetSequence +
                            process.treeDumper)

if option=='RECO':
    process.analysis.replace(process.leptonSequence, process.goodOfflinePrimaryVertex + process.leptonSequence)

process.load("V3Resonances.EDBRCommon.data.RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8")
process.source.inputCommands = ['keep *','drop *_isolatedTracks_*_*']
process.source.fileNames = [
'/store/data/Run2016B/JetHT/RAW/v2/000/274/443/00000/F4AFF614-1E2C-E611-94A9-02163E0141E3.root',
]

process.maxEvents.input = 5000
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.limit = 99999999

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("RStreeEDBR_pickup.root")
                                   )
