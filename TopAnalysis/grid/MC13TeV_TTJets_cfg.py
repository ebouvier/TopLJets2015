from WMCore.Configuration import Configuration
import os
config = Configuration()

config.section_("General")
config.General.requestName = "MC13TeV_TTJets"
config.General.workArea = "grid"
config.General.transferOutputs=True

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "/afs/cern.ch/user/b/byates/CMSSW_7_6_3/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py"
config.JobType.disableAutomaticOutputCollection = False
config.JobType.pyCfgParams = ['runOnData=False']
config.JobType.inputFiles = ['Fall15_25nsV2_MC.db']

config.section_("Data")
config.Data.inputDataset = "/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/MINIAODSIM"
config.Data.inputDBS = "global"
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 2
config.Data.publication = True
config.Data.ignoreLocality = False
config.Data.outLFNDirBase = '/store/user/byates/402de72/'

config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
